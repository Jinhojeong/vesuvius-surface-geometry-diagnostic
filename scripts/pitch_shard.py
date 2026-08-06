"""Sharded crossing-position extraction. ray_metrics is reused verbatim from
iyando's pitch_qa.py; slab_canvases is copied below with ONE change, hoisting
the invariant labels.max() out of the per-key loop guard (the original calls
it once per dict key, ~16k full-array scans per slab, 98% of runtime by
cProfile; outputs are bit-identical, checked by equiv_check.py). Each shard
handles every n-th slab and also emits the per-slice origin."""
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/mnt/vesuvius/vesuvius-sheet-tools/scripts")
import pitch_qa as PQ


def slab_canvases_fast(run: Path, z0: int):
    """pitch_qa.slab_canvases with the labels.max() hoist. Same output."""
    tdir = run / "blocks" / f"z{z0}"
    tpath = run / f"stitch_table_z{z0}.json"
    gpath = run / "global_table.json"
    if not tdir.is_dir() or not tpath.exists():
        return None
    if "t" not in PQ._GLOBAL_TABLE_CACHE:
        PQ._GLOBAL_TABLE_CACHE["t"] = (
            json.load(open(gpath)) if gpath.exists() else None)
    global_table = PQ._GLOBAL_TABLE_CACHE["t"]
    slab_table = json.load(open(tpath))
    canvases = {zl: np.zeros((PQ.GRID, PQ.GRID), dtype=np.int64)
                for zl in PQ.Z_LOCALS}
    tiles = {}
    for p in tdir.glob("tile_*.npz"):
        m = PQ.TILE_RE.fullmatch(p.name)
        tiles[(int(m.group(1)), int(m.group(2)))] = p
    for (y0, x0), p in sorted(tiles.items()):
        tkey = f"y{y0}_x{x0}"
        src = None
        if global_table is not None:
            src = global_table.get(f"z{z0}/{tkey}")
        if src is None:
            src = slab_table.get(tkey)
        if not src:
            continue
        with np.load(p) as d:
            labels = d["labels"]
        lmax = int(labels.max())          # hoisted: was recomputed per key
        lut = np.zeros(lmax + 1, dtype=np.int64)
        for k, v in src.items():
            if int(k) <= lmax:
                lut[int(k)] = v
        own_y = PQ.OVERLAP // 2 if (y0 - PQ.STRIDE, x0) in tiles else 0
        own_x = PQ.OVERLAP // 2 if (y0, x0 - PQ.STRIDE) in tiles else 0
        for zl in PQ.Z_LOCALS:
            sl2d = lut[labels[min(zl, labels.shape[0] - 1), own_y:, own_x:]]
            canvases[zl][y0 + own_y : y0 + own_y + sl2d.shape[0],
                         x0 + own_x : x0 + own_x + sl2d.shape[1]] = sl2d
    return canvases


def main():
    run = Path(sys.argv[1])
    shard, n_shards = int(sys.argv[2]), int(sys.argv[3])
    outdir = Path(sys.argv[4])
    outdir.mkdir(parents=True, exist_ok=True)
    z0s = sorted(int(p.name[1:]) for p in (run / "blocks").iterdir()
                 if p.name.startswith("z") and 896 <= int(p.name[1:]) <= 10976)
    mine = z0s[shard::n_shards]
    print(f"shard {shard}/{n_shards}: slabs {mine}", flush=True)
    pos_f = open(outdir / f"pos_{shard}.csv", "w", newline="")
    org_f = open(outdir / f"org_{shard}.csv", "w", newline="")
    pw, ow = csv.writer(pos_f), csv.writer(org_f)
    for z0 in mine:
        canv = slab_canvases_fast(run, z0)
        if canv is None:
            print(f"z{z0}: no canvases", flush=True)
            continue
        n_pos = 0
        for zl, canvas in sorted(canv.items()):
            res = PQ.ray_metrics(canvas, 6.0)
            if res is None:
                continue
            metrics, cents, (cy, cx) = res
            ow.writerow([z0 + zl, round(float(cy), 3), round(float(cx), 3)])
            for i in range(len(metrics)):
                for k, r in enumerate(cents[i]):
                    pw.writerow([z0 + zl, i * 6.0, k, round(float(r), 2),
                                 round(float(r) * PQ.UM_PER_VOX, 1)])
                    n_pos += 1
        pos_f.flush(); org_f.flush()
        print(f"z{z0} done ({n_pos} positions)", flush=True)
    pos_f.close(); org_f.close()
    print(f"SHARD_{shard}_DONE", flush=True)


if __name__ == "__main__":
    main()
