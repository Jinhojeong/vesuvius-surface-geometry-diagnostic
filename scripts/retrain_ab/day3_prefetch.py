"""Prereg build, day 3a: prefetch the CT for every frozen crop.

Training reads CT + labels. Labels live locally in the two arm trees; the CT
streams from the public S3 zarr at the label grid's level (L1, 17.28 um). This
fetches each frozen crop's 160^3 CT block once into a local cache so the twelve
seed runs never touch the network. Resumable, one .npy per crop key.

Alignment gate: for the first 10 crops the mean CT inside the label mask must
exceed the mean outside (papyrus is bright). A crop failing that prints ALIGN
FAIL and the run exits nonzero at the end, because a coordinate-frame mistake
here would silently train both arms on wrong voxels.
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path

import numpy as np
import zarr

S3 = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
CT_URL = f"{S3}/PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr"
LEVEL = 1
FRZ = Path("/mnt/vesuvius/experiments/retrain_ab/frozen")
V2 = Path("/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired")
CACHE = Path("/mnt/vesuvius/experiments/retrain_ab/ct_cache")
CACHE.mkdir(parents=True, exist_ok=True)
CROP = 160


def main() -> None:
    crops = json.loads((FRZ / "crops.json").read_text())["crops"]
    ct = zarr.open_array(f"{CT_URL}/{LEVEL}", mode="r")
    print(f"ct zarr L{LEVEL} shape={ct.shape} dtype={ct.dtype}", flush=True)

    origins: dict[str, tuple[int, int, int]] = {}
    align_fail = 0
    t0 = time.time()
    done = 0
    for i, c in enumerate(crops):
        op = CACHE / f"{c['key']}.npy"
        if op.exists():
            done += 1
            continue
        f = c["file"]
        if f not in origins:
            with np.load(V2 / f) as d:
                origins[f] = (int(d["z0"]), int(d["y0"]), int(d["x0"]))
        oz, oy, ox = origins[f]
        z0, y0, x0 = c["corner"]
        gz, gy, gx = oz + z0, oy + y0, ox + x0
        blk = np.asarray(ct[gz:gz + CROP, gy:gy + CROP, gx:gx + CROP])
        if blk.shape != (CROP, CROP, CROP):
            pad = [(0, CROP - s) for s in blk.shape]
            blk = np.pad(blk, pad)
        np.save(op, blk)
        done += 1
        if i < 10:
            with np.load(V2 / f) as d:
                lab = d["labels"][z0:z0 + CROP, y0:y0 + CROP, x0:x0 + CROP]
            m = lab > 0
            if m.any() and (~m).any():
                inside, outside = float(blk[m].mean()), float(blk[~m].mean())
                ok = inside > outside
                align_fail += 0 if ok else 1
                print(f"ALIGN {'ok' if ok else 'FAIL'} {c['tile']} "
                      f"in={inside:.1f} out={outside:.1f}", flush=True)
        if done % 100 == 0:
            el = time.time() - t0
            print(f"[{done}/{len(crops)}] {el/60:.0f}min", flush=True)
    print(f"DONE {done}/{len(crops)} align_fail={align_fail}", flush=True)
    sys.exit(1 if align_fail else 0)


if __name__ == "__main__":
    main()
