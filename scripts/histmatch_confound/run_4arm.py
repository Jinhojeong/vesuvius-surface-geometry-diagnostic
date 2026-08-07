"""Four-arm 2x2: {original, histogram-matched} x {plans CTNormalization, instance z-score}.

Why this exists. Diego-dcv preregistered (villa #191, comment 5213508599) that
histogram-matched copies of the located volumes should recover PART of the located
vs non-located recall gap, and that the recovered share measures a brightness
effect while the surviving share measures intrinsic contrast. That decomposition
assumes the recovery route is the model's response to level. villa's
`vesuvius.predict` normalises nnU-Net checkpoints with instance z-score whatever
the plans say (inference.py:225 CLI default, :388 nnunet branch never reads the
plans' normalization scheme), and instance z-score standardises by the input's OWN
statistics. Histogram-matching a located volume toward the non-located reference
moves it from about +0.99 sigma off the plans intensity fingerprint to about
-0.36 sigma, i.e. much closer to correctly normalised BY CONSTRUCTION. So on a
wrapper-path bench part of any "recovery" is the normalization being accidentally
fixed rather than the model tolerating a level shift.

The 2x2 separates the two. B->D is what a wrapper-path bench measures as recovery
from matching. A->C is what the model does with a level-matched input when it is
normalised the way the plans prescribe. The difference is the size of the confound.

Geometry is identical to run_shells.py (which copies TAUIL-Abd-Elilah's
m7_margin_fp.py @ 9afa412): centred 256^3 block out of the 320^3 uint8 volume,
nnUNetPredictor sliding window at patch 192^3, tile_step_size 0.5, gaussian
weighting, mirroring off, fp16 autocast, TRIM=64 so the scored crop is the inner
128^3, mask sigmoid(logit1 - logit0) > 0.2, class 2 excluded.

Two things are deliberately held fixed across all four arms so only the network
input varies:
  * the labels, and therefore n_sheet / n_scored / the shells;
  * the margin mask, which is computed from the ORIGINAL CT in every arm. The
    histogram LUT is monotone but not affine, so it perturbs the Hessian that
    across_sheet_dirs eigendecomposes. Freezing the margin on the original CT
    keeps the label-side geometry byte-identical to control_his60.jsonl.

Histogram LUTs are READ from /mnt/vesuvius/kaggle892/histmatch_check.json, which
histmatch892.py built against the pooled non-located reference. Nothing here
recomputes the reference.

recall / precision / pred-positive-fraction are derived from the same counts the
existing artifacts carry, so they are comparable line for line:
    TP        = n_pred_scored - n_fp
    recall    = TP / n_sheet
    precision = TP / n_pred_scored
    ppf       = n_pred_scored / n_scored

One JSON line per (volume, arm) so an ssh drop costs at most the arm in flight.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path

import numpy as np
import tifffile
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geom import margin_mask, distance_profile, analyse, descriptives, TRIM, SIZE, THRESH

IM = Path("/mnt/vesuvius/kaggle892/images")
LB = Path("/mnt/vesuvius/kaggle892/labels")
MODEL = "/mnt/vesuvius/models/surface_m7_nnunet"
HISTMATCH = Path("/mnt/vesuvius/kaggle892/histmatch_check.json")

ARMS = ["A_orig_ctplans", "B_orig_instzs", "C_match_ctplans", "D_match_instzs"]
# E and F are the affine control. The histogram LUT is monotone but not affine, so
# arm C mixes a level/scale change with the LUT's rank compression. E rescales the
# original block by the single affine map that reproduces the matched copy's block
# mean and std exactly, so E and C carry the same first two moments and differ only
# in the nonlinear part. F is the same input under instance z-score, which is
# algebraically invariant to any affine map, so F must come back identical to B;
# that is a self-test of the invariance rather than a new measurement.
AFFINE_ARMS = ["E_affine_ctplans", "F_affine_instzs"]
OOM_RETRIES = 40
OOM_SLEEP = 60


class PatchInstanceZ(torch.nn.Module):
    """villa's instance_zscore is applied by Volume.__getitem__ (volume.py:853) to
    whatever slice VCDataset asked for, and VCDataset asks per patch
    (vc_dataset.py:461 __getitem__ -> Volume slice). So the wrapper-path per-patch
    convention normalises each 192^3 patch by its own mean and std. This wrapper
    reproduces that inside nnUNet's own sliding window."""

    def __init__(self, net: torch.nn.Module):
        super().__init__()
        self.net = net

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        dims = tuple(range(1, x.ndim))
        m = x.mean(dim=dims, keepdim=True)
        s = x.std(dim=dims, unbiased=False, keepdim=True).clamp(min=1e-8)
        return self.net((x - m) / s)


def build_predictor(device: str = "cuda"):
    from nnunetv2.inference.predict_from_raw_data import nnUNetPredictor
    p = nnUNetPredictor(tile_step_size=0.5, use_gaussian=True, use_mirroring=False,
                        perform_everything_on_device=True, device=torch.device(device),
                        verbose=False, verbose_preprocessing=False, allow_tqdm=False)
    p.initialize_from_trained_model_folder(MODEL, use_folds=(0,),
                                           checkpoint_name="checkpoint_best.pth")
    p.network.load_state_dict(p.list_of_parameters[0])
    p.network = p.network.to(p.device).eval()
    plans = json.loads((Path(MODEL) / "plans.json").read_text())
    ip = plans["foreground_intensity_properties_per_channel"]["0"]
    return p, ip


def ct_normalize(x: np.ndarray, ip: dict) -> np.ndarray:
    """nnunetv2 CTNormalization, verbatim maths."""
    x = x.astype(np.float32, copy=True)
    np.clip(x, ip["percentile_00_5"], ip["percentile_99_5"], out=x)
    x -= float(ip["mean"])
    x /= max(float(ip["std"]), 1e-8)
    return x


def instance_zscore(x: np.ndarray) -> np.ndarray:
    """villa volume.py:853-858 verbatim maths, applied to the array that is fed."""
    x = x.astype(np.float32, copy=True)
    x -= float(x.mean())
    x /= max(float(x.std()), 1e-8)
    return x


def centre_block(ct: np.ndarray) -> tuple[np.ndarray, int]:
    off = (ct.shape[0] - SIZE) // 2
    return ct[off:off + SIZE, off:off + SIZE, off:off + SIZE], off


def predict_prob(blk: np.ndarray, pred_obj, ip: dict, norm: str, inst_mode: str, off: int):
    """norm in {ct_plans, instance_zscore}. inst_mode in {volume, patch} and only
    bites when norm == instance_zscore."""
    patch_wrapped = False
    if norm == "ct_plans":
        x = ct_normalize(blk, ip)
    elif norm == "instance_zscore":
        if inst_mode == "volume":
            x = instance_zscore(blk)
        elif inst_mode == "patch":
            # feed raw intensities; the wrapper z-scores each 192^3 patch
            x = blk.astype(np.float32, copy=True)
            pred_obj.network = PatchInstanceZ(pred_obj.network)
            patch_wrapped = True
        else:
            raise ValueError(inst_mode)
    else:
        raise ValueError(norm)
    try:
        # the box is shared; another job's allocation can transiently starve this one.
        # Back off and retry rather than dropping the volume, so the cohort stays whole.
        for attempt in range(OOM_RETRIES):
            try:
                with torch.no_grad():
                    lg = pred_obj.predict_sliding_window_return_logits(
                        torch.from_numpy(x)[None])
                break
            except torch.OutOfMemoryError:
                torch.cuda.empty_cache()
                if attempt == OOM_RETRIES - 1:
                    raise
                print(f"MARKER oom retry {attempt+1}/{OOM_RETRIES}, sleeping "
                      f"{OOM_SLEEP}s", flush=True)
                time.sleep(OOM_SLEEP)
    finally:
        if patch_wrapped:
            pred_obj.network = pred_obj.network.net
    lg = lg.float().cpu().numpy()
    lo, hi = TRIM, SIZE - TRIM
    l0 = lg[0, lo:hi, lo:hi, lo:hi]
    l1 = lg[1, lo:hi, lo:hi, lo:hi]
    p = 1.0 / (1.0 + np.exp(-(l1 - l0)))
    return p, off + lo, off + hi


def rates(desc: dict) -> dict:
    tp = desc["n_pred_scored"] - desc["n_fp"]
    return {"recall": tp / desc["n_sheet"] if desc["n_sheet"] else None,
            "precision": tp / desc["n_pred_scored"] if desc["n_pred_scored"] else None,
            "ppf": desc["n_pred_scored"] / desc["n_scored"] if desc["n_scored"] else None,
            "tp": int(tp)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", required=True, help="json file: list of sample ids")
    ap.add_argument("--out", required=True)
    ap.add_argument("--arms", default=",".join(ARMS))
    ap.add_argument("--inst-mode", default="volume", choices=("volume", "patch"))
    ap.add_argument("--device", default="cuda")
    a = ap.parse_args()

    names = json.loads(Path(a.names).read_text())
    arms = a.arms.split(",")
    luts = json.loads(HISTMATCH.read_text())["luts"]
    outp = Path(a.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if outp.exists():
        for ln in outp.read_text().splitlines():
            if ln.strip():
                r = json.loads(ln)
                done.add((r["sample"], r["arm"], r.get("inst_mode")))

    pred_obj, ip = build_predictor(a.device)
    print(f"MARKER init ok intensity={ip} arms={arms} inst_mode={a.inst_mode} "
          f"n={len(names)}", flush=True)

    t0 = time.time()
    fh = outp.open("a")
    for k, nm in enumerate(names):
        todo = [arm for arm in arms
                if (nm, arm, a.inst_mode if "instzs" in arm else None) not in done]
        if not todo:
            continue
        ct = np.asarray(tifffile.imread(str(IM / f"{nm}.tif")))
        lab = np.asarray(tifffile.imread(str(LB / f"{nm}.tif")))
        mg_full = margin_mask(ct, lab)  # frozen on the ORIGINAL ct for every arm

        blk_o, off = centre_block(ct)
        stats = {"orig_blk_mean": float(blk_o.mean()), "orig_blk_std": float(blk_o.std()),
                 "plans_mean": float(ip["mean"]), "plans_std": float(ip["std"])}
        matched = None
        blk_aff = None
        if any(arm[0] in "CDEF" and arm[1] == "_" for arm in todo):
            lut = np.asarray(luts[nm], dtype=np.uint8)
            matched = lut[ct]
            blk_m, _ = centre_block(matched)
            stats["match_blk_mean"] = float(blk_m.mean())
            stats["match_blk_std"] = float(blk_m.std())
            stats["lut_distinct_levels"] = int(np.unique(lut).size)
            occ = np.unique(blk_o)
            stats["lut_distinct_levels_on_block"] = int(np.unique(lut[occ]).size)
            stats["block_distinct_levels"] = int(occ.size)
            # affine copy: same block mean and std as the matched copy, float32,
            # no clipping and no quantisation, so the map stays exactly affine
            a_ = stats["match_blk_std"] / max(stats["orig_blk_std"], 1e-8)
            b_ = stats["match_blk_mean"] - a_ * stats["orig_blk_mean"]
            stats["affine_a"], stats["affine_b"] = a_, b_
            blk_aff = blk_o.astype(np.float32) * a_ + b_
        for key, pre in (("orig", "orig"), ("match", "match")):
            if f"{pre}_blk_mean" in stats:
                stats[f"{pre}_sigma"] = ((stats[f"{pre}_blk_mean"] - stats["plans_mean"])
                                         / stats["plans_std"])

        for arm in todo:
            t1 = time.time()
            if arm.startswith(("E_", "F_")):
                blk = blk_aff
            else:
                blk, _ = centre_block(ct if arm.startswith(("A_", "B_")) else matched)
            norm = "ct_plans" if "ctplans" in arm else "instance_zscore"
            p, lo, hi = predict_prob(blk, pred_obj, ip, norm, a.inst_mode, off)
            sl = (slice(lo, hi),) * 3
            labc, mg = lab[sl], mg_full[sl]
            pred = p > THRESH
            r = analyse(labc, mg, pred)
            r["shells"] = distance_profile(labc, pred)
            r["desc"] = descriptives(labc, pred)
            r.update(rates(r["desc"]))
            r.update({"sample": nm, "arm": arm, "norm": norm,
                      "inst_mode": a.inst_mode if norm == "instance_zscore" else None,
                      "crop": [lo, hi], "seconds": round(time.time() - t1, 1)})
            r.update(stats)
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            # volumes with no labelled sheet in the crop score None; the row is already
            # written, so the print must not be what ends the cohort
            fm = lambda v: "none" if v is None else f"{v:.4f}"
            print(f"MARKER [{k+1}/{len(names)}] {nm} {arm} status={r['status']} "
                  f"recall={fm(r['recall'])} prec={fm(r['precision'])} "
                  f"ppf={fm(r['ppf'])} {time.time()-t0:.0f}s", flush=True)
    fh.close()
    print("MARKER done", flush=True)


if __name__ == "__main__":
    main()
