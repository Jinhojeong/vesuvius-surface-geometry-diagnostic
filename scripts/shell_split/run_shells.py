"""m7 Euclidean-shell false-positive placement, ours.

Geometry is copied from TAUIL-Abd-Elilah's m7_margin_fp.py @ 9afa412 and must stay copied:
centred 256^3 bbox, TRIM=64 so the scored crop is the inner 128^3, prediction mask
sigmoid(logit1 - logit0) > 0.2, class 2 excluded, one-voxel Euclidean shells off the
labelled sheet.

The inference path is NOT his. This runs the same published checkpoint
(hf scrollprize/surface_m7_nnunet, fold_0/checkpoint_best.pth) through nnUNetv2 2.8.1's
own sliding window, which is the path the model card prescribes:

  * crop the centred 256^3 block out of the 320^3 uint8 volume FIRST, so the network sees
    exactly the voxels his --bbox asked for
  * CTNormalization with the constants in plans.json (clip to [p00.5, p99.5], then
    (x - mean) / std). These are global fingerprint constants, not per-image statistics,
    so cropping before normalising cannot change them
  * nnUNetPredictor.predict_sliding_window_return_logits: patch 192^3, tile_step_size 0.5,
    gaussian importance weighting, mirroring OFF (he passed --disable_tta), fp16 autocast
    (that is nnUNet's own default inference path, not a concession to the 11GB card)
  * logits averaged, then sigmoid(l1 - l0). Equivalently p1 / (p0 + p1) of the 3-class
    softmax; the ignore class drops out of the pairwise sigmoid either way.

Writes one JSON line per volume so an ssh drop costs at most the volume in flight.
"""
from __future__ import annotations
import argparse, json, os, sys, time
from pathlib import Path

import numpy as np
import tifffile
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geom import margin_mask, distance_profile, analyse, descriptives, TRIM, SIZE, THRESH

IM = Path("/mnt/vesuvius/kaggle892/images")
LB = Path("/mnt/vesuvius/kaggle892/labels")
MODEL = "/mnt/vesuvius/models/surface_m7_nnunet"


def build_predictor(device: str = "cuda"):
    from nnunetv2.inference.predict_from_raw_data import nnUNetPredictor
    p = nnUNetPredictor(tile_step_size=0.5, use_gaussian=True, use_mirroring=False,
                        perform_everything_on_device=True, device=torch.device(device),
                        verbose=False, verbose_preprocessing=False, allow_tqdm=False)
    p.initialize_from_trained_model_folder(MODEL, use_folds=(0,),
                                           checkpoint_name="checkpoint_best.pth")
    params = p.list_of_parameters[0]
    p.network.load_state_dict(params)
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


def predict_prob(ct: np.ndarray, pred_obj, ip: dict) -> tuple[np.ndarray, int, int]:
    off = (ct.shape[0] - SIZE) // 2
    blk = ct[off:off + SIZE, off:off + SIZE, off:off + SIZE]
    x = torch.from_numpy(ct_normalize(blk, ip))[None]
    with torch.no_grad():
        lg = pred_obj.predict_sliding_window_return_logits(x)
    lg = lg.float().cpu().numpy()
    lo, hi = TRIM, SIZE - TRIM
    l0 = lg[0, lo:hi, lo:hi, lo:hi]
    l1 = lg[1, lo:hi, lo:hi, lo:hi]
    p = 1.0 / (1.0 + np.exp(-(l1 - l0)))
    return p, off + lo, off + hi


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", required=True, help="json file: list of sample ids")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cache", default=None)
    ap.add_argument("--device", default="cuda")
    a = ap.parse_args()

    names = json.loads(Path(a.names).read_text())
    outp = Path(a.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if outp.exists():
        for ln in outp.read_text().splitlines():
            if ln.strip():
                done.add(json.loads(ln)["sample"])
    pred_obj, ip = build_predictor(a.device)
    print(f"MARKER init ok  intensity={ip}  todo={len(names)-len(done)}", flush=True)

    t0 = time.time()
    fh = outp.open("a")
    for k, nm in enumerate(names):
        if nm in done:
            continue
        t1 = time.time()
        ct = np.asarray(tifffile.imread(str(IM / f"{nm}.tif")))
        lab = np.asarray(tifffile.imread(str(LB / f"{nm}.tif")))
        p, lo, hi = predict_prob(ct, pred_obj, ip)
        sl = (slice(lo, hi),) * 3
        labc = lab[sl]
        mg = margin_mask(ct, lab)[sl]
        pred = p > THRESH
        r = analyse(labc, mg, pred)
        r["shells"] = distance_profile(labc, pred)
        r["desc"] = descriptives(labc, pred)
        r["sample"] = nm
        r["crop"] = [lo, hi]
        r["seconds"] = round(time.time() - t1, 1)
        fh.write(json.dumps(r) + "\n")
        fh.flush()
        if a.cache:
            Path(a.cache).mkdir(parents=True, exist_ok=True)
            np.save(Path(a.cache) / f"{nm}.npy", p.astype(np.float16))
        s1 = (r.get("shells") or {}).get("shell_1")
        print(f"MARKER [{k+1}/{len(names)}] {nm} status={r['status']} "
              f"enr={r.get('enrichment')} shell1={s1} {time.time()-t0:.0f}s", flush=True)
    fh.close()
    print("MARKER done", flush=True)


if __name__ == "__main__":
    main()
