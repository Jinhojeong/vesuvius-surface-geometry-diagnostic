"""How much voxel mass does plans CTNormalization clip away in each input arm?

CTNormalization clips to [percentile_00_5, percentile_99_5] = [0, 212] before the
affine step. If the affine control (arm E) pushed a large share of the block below
0, its recall loss would be a clipping artifact rather than a level response, and
the E vs C comparison would not be clean. This measures that directly.

CPU only, no model, no GPU. Reads the same tifs and the same per-volume LUTs and
affine constants that run_4arm.py used (the constants are read back out of the
jsonl the run wrote, so this cannot drift from what was actually fed).
"""
from __future__ import annotations
import json
import statistics as st
from pathlib import Path

import numpy as np
import tifffile

D = Path("/mnt/vesuvius/experiments/histmatch_confound")
IM = Path("/mnt/vesuvius/kaggle892/images")
SIZE = 256
LO, HI = 0.0, 212.0  # plans percentile_00_5 / percentile_99_5


def main() -> None:
    rows = {}
    for f in ("loc60_4arm.jsonl", "loc60_affine.jsonl"):
        for ln in (D / f).read_text().splitlines():
            if ln.strip():
                r = json.loads(ln)
                rows.setdefault(r["sample"], {}).update(r)

    out = []
    for s, r in sorted(rows.items()):
        if "affine_a" not in r:
            continue
        ct = np.asarray(tifffile.imread(str(IM / f"{s}.tif")))
        off = (ct.shape[0] - SIZE) // 2
        blk = ct[off:off + SIZE, off:off + SIZE, off:off + SIZE].astype(np.float32)
        n = blk.size
        aff = blk * r["affine_a"] + r["affine_b"]
        out.append({
            "sample": s,
            "affine_a": round(r["affine_a"], 4), "affine_b": round(r["affine_b"], 3),
            "orig_frac_below_lo": float((blk < LO).sum()) / n,
            "orig_frac_above_hi": float((blk > HI).sum()) / n,
            "affine_frac_below_lo": float((aff < LO).sum()) / n,
            "affine_frac_above_hi": float((aff > HI).sum()) / n,
        })

    # matched copies: apply the shipped LUT, which is uint8 -> uint8 so it cannot
    # leave [0, 255]; only the upper clip at 212 can bite
    luts = json.loads((Path("/mnt/vesuvius/kaggle892/histmatch_check.json")).read_text())["luts"]
    for row in out:
        s = row["sample"]
        ct = np.asarray(tifffile.imread(str(IM / f"{s}.tif")))
        off = (ct.shape[0] - SIZE) // 2
        m = np.asarray(luts[s], dtype=np.uint8)[ct][off:off + SIZE, off:off + SIZE,
                                                    off:off + SIZE]
        row["match_frac_above_hi"] = float((m > HI).sum()) / m.size
        row["match_frac_below_lo"] = 0.0

    summary = {k: round(st.median([r[k] for r in out]), 6)
               for k in out[0] if k not in ("sample",)}
    res = {"clip_bounds": [LO, HI], "n": len(out), "medians": summary, "per_volume": out}
    (D / "clipcheck.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(summary, indent=1))
    print("worst affine_frac_below_lo:",
          round(max(r["affine_frac_below_lo"] for r in out), 6))
    print("worst affine_frac_above_hi:",
          round(max(r["affine_frac_above_hi"] for r in out), 6))
    print("worst orig_frac_above_hi  :",
          round(max(r["orig_frac_above_hi"] for r in out), 6))
    print("wrote", D / "clipcheck.json")


if __name__ == "__main__":
    main()
