"""Leg-2 range check: `surface_recto_3dunet` over the 892 public volumes.

Why this exists. TAUIL's three-gate conjunction gives leg 2 an independent vote,
and I committed to establishing whether it has one before the preregistration is
written. In July this checkpoint returned pooled dice 0.23 on Dataset059 patches,
which is chance at a ground-truth fill near 25 percent, and a dilation probe ruled
out a small-offset explanation. That was a different population, so it does not
settle the question on his 892. This does.

Geometry is leg 1's, unchanged, so the two legs land on the same voxels per volume
and can be compared row by row: centred 256^3 block out of the 320^3 volume,
TRIM=64 so the scored crop is the inner 128^3, label class 2 excluded.

Normalisation is NOT m7's. This checkpoint declares `zscore` with an empty
`intensity_properties` block, so CT normalisation is not defined for it and per
volume z-scoring is the only scheme it supports. Running it under m7's constants
would be measuring the wrong thing. The declared patch is 256^3 and the block is
256^3, so each volume is one forward pass with no tiling and no blending seam.

Nothing is thresholded here that cannot be re-thresholded later. Each row carries a
2000-bin probability histogram over the scored region and over the labelled sheet,
which is enough to recover recall, precision and dice at any operating point, an
oracle-best threshold, and a precision ceiling, without touching the GPU again.
"""
from __future__ import annotations
import argparse, json, os, sys, time
from pathlib import Path

import numpy as np
import tifffile
import torch

MODEL_DIR = "/mnt/vesuvius/models/surface_recto_3dunet"
IM = Path("/mnt/vesuvius/kaggle892/images")
LB = Path("/mnt/vesuvius/kaggle892/labels")
GROUPS = Path("/mnt/vesuvius/kaggle892/groups892.json")

SIZE, TRIM = 256, 64
PS = 256                                   # the checkpoint's own declared patch
THRESH_GRID = (0.05, 0.10, 0.20, 0.30, 0.50, 0.70)   # leg 1 carries these six
NBIN = 2000


def load_model():
    """vesuvius' own inference recipe, as in diag191.py. The state dict arrives with
    `_orig_mod.` prefixes from a compiled training run."""
    from vesuvius.models.build.build_network_from_config import NetworkFromConfig
    ckpt = torch.load(os.path.join(MODEL_DIR, "checkpoint_inference_ready.pth"),
                      map_location="cpu", weights_only=False)
    mc = ckpt.get("model_config") or json.load(
        open(os.path.join(MODEL_DIR, "config.json")))

    class CM:
        def __init__(s, mc):
            s.model_config = mc
            s.targets = mc.get("targets", {})
            s.train_patch_size = mc.get("train_patch_size", mc.get("patch_size"))
            s.train_batch_size = mc.get("train_batch_size", mc.get("batch_size", 2))
            s.in_channels = mc.get("in_channels", 1)
            s.autoconfigure = mc.get("autoconfigure", False)
            s.model_name = mc.get("model_name", "Model")
            s.spacing = [1] * len(s.train_patch_size)

    net = NetworkFromConfig(CM(mc))

    def strip(k):
        for p in ("module.", "_orig_mod."):
            while k.startswith(p):
                k = k[len(p):]
        return k

    sd = {strip(k): v for k, v in (ckpt.get("model", ckpt)).items()}
    net.load_state_dict(sd, strict=True)
    return net.cuda().eval(), mc, ckpt


@torch.no_grad()
def predict_block(net, blk: np.ndarray) -> np.ndarray:
    """One forward pass over the 256^3 block. Per-volume z-score, fp16 autocast,
    softmax over the two output channels."""
    a = blk.astype(np.float32)
    a = (a - a.mean()) / (a.std() + 1e-6)
    t = torch.from_numpy(a[None, None]).cuda()
    with torch.autocast("cuda", torch.float16):
        out = net(t)
    if isinstance(out, dict):
        out = out.get("surface", list(out.values())[0])
    if isinstance(out, (list, tuple)):
        out = out[0]
    return torch.softmax(out.float(), 1)[0, 1].cpu().numpy()


def score(nm: str, p: np.ndarray, labc: np.ndarray, group: str) -> dict:
    scored = labc != 2
    sheet = labc == 1
    n_scored, n_sheet = int(scored.sum()), int(sheet.sum())
    row = {
        "sample": nm,
        "group": group,
        "population": "located" if group in ("located", "intersecting", "iou1")
                      else "nonlocated",
        "n_scored": n_scored,
        "n_sheet": n_sheet,
        "n_nonsheet_scored": n_scored - n_sheet,
        "base_rate": (n_sheet / n_scored) if n_scored else None,
        "crop": [TRIM, SIZE - TRIM],
    }
    if not n_scored:
        row["status"] = "degenerate"
        row["reason"] = "no_scored_voxels_in_crop: the whole inner 128^3 is class 2"
        return row
    ps_ = p[scored]
    psh = p[sheet] if n_sheet else np.empty(0, np.float32)
    row["hist_scored"] = np.histogram(ps_, bins=NBIN, range=(0.0, 1.0))[0].tolist()
    row["hist_sheet"] = np.histogram(psh, bins=NBIN, range=(0.0, 1.0))[0].tolist()
    row["p_mean_scored"] = float(ps_.mean())
    row["p_mean_sheet"] = float(psh.mean()) if n_sheet else None
    sweep = {}
    for t in THRESH_GRID:
        npred = int((ps_ > t).sum())
        ntp = int((psh > t).sum()) if n_sheet else 0
        sweep[f"{t:.2f}"] = {"n_pred_scored": npred, "n_tp": ntp}
    row["threshold_sweep"] = sweep
    if not n_sheet:
        row["status"] = "degenerate"
        row["reason"] = "no_labelled_sheet_in_crop: recall and dice are undefined"
        return row
    row["status"] = "ok"
    row["reason"] = None
    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--device", default="cuda")
    a = ap.parse_args()

    names = json.loads(Path(a.names).read_text())
    groups = json.loads(GROUPS.read_text())
    outp = Path(a.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if outp.exists():
        for ln in outp.read_text().splitlines():
            if ln.strip():
                done.add(json.loads(ln)["sample"])
    todo = [n for n in names if n not in done]

    net, mc, ckpt = load_model()
    print(f"MARKER init ok  todo={len(todo)}  already={len(done)}", flush=True)
    print(f"MARKER model patch={mc.get('patch_size')} ps_used={PS} "
          f"norm={ckpt.get('normalization_scheme')} "
          f"intensity_properties={ckpt.get('intensity_properties')}", flush=True)
    print(f"MARKER geom SIZE={SIZE} TRIM={TRIM} grid={THRESH_GRID} nbin={NBIN}",
          flush=True)

    # his tables carry pred_on_empty_ct; measure it once, as for leg 1
    zero = np.zeros((SIZE, SIZE, SIZE), dtype=np.uint8)
    pz = predict_block(net, zero)
    print(f"MARKER pred_on_empty_ct={float((pz > 0.5).mean()):.8f} "
          f"max_p={float(pz.max()):.6f}", flush=True)
    del zero, pz

    lo, hi = TRIM, SIZE - TRIM
    fh = outp.open("a")
    t_start = time.time()
    for i, nm in enumerate(todo):
        t0 = time.time()
        ct = np.asarray(tifffile.imread(str(IM / f"{nm}.tif")))
        lab = np.asarray(tifffile.imread(str(LB / f"{nm}.tif")))
        off = (ct.shape[0] - SIZE) // 2
        blk = ct[off:off + SIZE, off:off + SIZE, off:off + SIZE]
        p = predict_block(net, blk)
        sl = slice(off + lo, off + hi)
        row = score(nm, p[lo:hi, lo:hi, lo:hi], lab[sl, sl, sl],
                    groups.get(nm, "nonlocated"))
        row["seconds"] = round(time.time() - t0, 2)
        fh.write(json.dumps(row) + "\n")
        fh.flush()
        if (i + 1) % 25 == 0:
            el = time.time() - t_start
            print(f"MARKER progress {i+1}/{len(todo)}  {el/60:.1f}min  "
                  f"eta {(el/(i+1))*(len(todo)-i-1)/60:.1f}min", flush=True)
    fh.close()
    print("MARKER STATUS=DONE", flush=True)


if __name__ == "__main__":
    main()
