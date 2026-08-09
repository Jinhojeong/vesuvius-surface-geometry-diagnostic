"""Prereg build, day 3c: precompute the boundary target for every crop and arm.

The smoke run spent most of each step decompressing a 256 MB label npz per
crop. The boundary target is deterministic given (tree, file, corner), so both
arms' targets are computed once here and packed with np.packbits to 512 KB per
crop. Training then reads CT (16 MB) plus a bitmap (0.5 MB) and the step time
returns to GPU-bound. 4,000 crops x 2 arms, CPU only, resumable.

Also emits frozen/flips.json, the per-crop count of differing target voxels
between arms over ALL 4,000 crops, which is the manipulation-check table the
prereg references (ab_check.json covered 20 crops; this is the full list).
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np

AB = Path("/mnt/vesuvius/experiments/retrain_ab")
TREES = {"v1": Path("/mnt/vesuvius/p1218_full/blocks"),
         "v2": Path("/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired")}
PS = 160


def boundary_target(lab: np.ndarray) -> np.ndarray:
    out = np.zeros(lab.shape, bool)
    fg = lab > 0
    for ax in range(3):
        a = [slice(None)] * 3
        b = [slice(None)] * 3
        a[ax] = slice(None, -1)
        b[ax] = slice(1, None)
        a, b = tuple(a), tuple(b)
        la, lb = lab[a], lab[b]
        diff = (la > 0) & (lb > 0) & (la != lb)
        out[a] |= diff
        out[b] |= diff
    return out & fg


def main() -> None:
    crops = json.loads((AB / "frozen/crops.json").read_text())["crops"]
    for arm in TREES:
        (AB / f"target_cache_{arm}").mkdir(exist_ok=True)

    by_file: dict[str, list[dict]] = {}
    for c in crops:
        by_file.setdefault(c["file"], []).append(c)

    flips: dict[str, int] = {}
    n_done = 0
    for f, cs in sorted(by_file.items()):
        need = [c for c in cs if not all(
            (AB / f"target_cache_{arm}" / f"{c['key']}.npy").exists()
            for arm in TREES)]
        if not need and all(c["key"] in flips for c in cs):
            continue
        labs = {}
        for arm, tree in TREES.items():
            with np.load(tree / f) as d:
                labs[arm] = d["labels"].astype(np.int32)
        for c in cs:
            z0, y0, x0 = c["corner"]
            sl = (slice(z0, z0 + PS), slice(y0, y0 + PS), slice(x0, x0 + PS))
            ts = {}
            for arm in TREES:
                t = boundary_target(labs[arm][sl])
                if t.shape != (PS, PS, PS):
                    # edge tiles are smaller than the standard block; outside
                    # the volume is background, matching the CT cache zero-pad
                    t = np.pad(t, [(0, PS - t.shape[k]) for k in range(3)])
                ts[arm] = t
                np.save(AB / f"target_cache_{arm}" / f"{c['key']}.npy",
                        np.packbits(t))
            flips[c["key"]] = int((ts["v1"] != ts["v2"]).sum())
            n_done += 1
        if n_done % 200 < len(cs):
            print(f"{n_done} crops cached", flush=True)

    stats = {k: v for k, v in flips.items()}
    by_kind: dict[str, list[int]] = {}
    for c in crops:
        if c["key"] in flips:
            by_kind.setdefault(c["kind"], []).append(flips[c["key"]])
    summary = {k: {"n": len(v), "median": float(np.median(v)),
                   "zero_flip": int(sum(1 for x in v if x == 0))}
               for k, v in by_kind.items()}
    (AB / "frozen/flips.json").write_text(json.dumps(
        {"summary": summary, "per_crop": stats}, indent=0))
    print(json.dumps(summary, indent=1), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
