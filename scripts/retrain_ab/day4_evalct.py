"""Prereg build, day 4a: cache the CT for every eval tile that carries sites.

The instrument predicts whole 256x512x512 eval tiles rather than per-site
windows, because sites average ~27 per tile and their 160^3 windows overlap
heavily, and a single tile pass also gives the pooled probability histogram the
budget rule needs. This fetches each needed tile's CT once from the S3 zarr at
L1 into /mnt/vesuvius/experiments/retrain_ab/evalct/. Resumable, CPU only.
"""
from __future__ import annotations
import json, time
from pathlib import Path

import numpy as np
import zarr

S3 = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
CT_URL = f"{S3}/PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr"
LEVEL = 1
AB = Path("/mnt/vesuvius/experiments/retrain_ab")
V2 = Path("/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired")
OUT = AB / "evalct"
OUT.mkdir(exist_ok=True)
SHAPE = (256, 512, 512)


def main() -> None:
    tiles = set()
    for f in ("sites_primary.json", "sites_onesided.json", "sites_offset.json",
              "sites_background.json", "sites_trainsplit.json"):
        for s in json.loads((AB / "frozen" / f).read_text())["sites"]:
            tiles.add(s["tile"])
    tiles = sorted(tiles)
    print(f"{len(tiles)} eval tiles carry frozen sites", flush=True)
    ct = zarr.open_array(f"{CT_URL}/{LEVEL}", mode="r")
    t0 = time.time()
    for i, name in enumerate(tiles):
        op = OUT / f"{name}.npy"
        if op.exists():
            continue
        slab, tile = name.split("_", 1)
        with np.load(V2 / slab / f"{tile}.npz") as d:
            oz, oy, ox = int(d["z0"]), int(d["y0"]), int(d["x0"])
        blk = np.asarray(ct[oz:oz + SHAPE[0], oy:oy + SHAPE[1], ox:ox + SHAPE[2]])
        if blk.shape != SHAPE:
            blk = np.pad(blk, [(0, SHAPE[k] - blk.shape[k]) for k in range(3)])
        np.save(op, blk)
        if (i + 1) % 20 == 0:
            print(f"[{i+1}/{len(tiles)}] {(time.time()-t0)/60:.0f}min", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
