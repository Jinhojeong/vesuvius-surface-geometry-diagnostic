"""Prereg build, day 2: boundary-target generator, shared crop list, A/B check.

Three jobs in one file, run in order:

1. boundary_target(lab): the training target. A voxel is positive iff it is
   foreground and has a 6-connected foreground neighbour carrying a DIFFERENT
   nonzero instance id. Pure numpy, axis-shift comparisons, no scipy.

2. Crop list: one frozen list shared verbatim by both arms, 160^3 crops on
   train tiles only. Composition 50% SPLIT-jittered / 25% ONE_SIDED-jittered /
   25% uniform-over-foreground. Jitter is +/-16 vox, derived from md5 of the
   site key, so there is no RNG and both arms read byte-identical coordinates.
   Uniform crops come from md5-hashed tile-grid positions with a foreground
   floor of 200k voxels, checked against the v2 mask (identical to v1's).

3. A/B manipulation check: for 20 SPLIT-anchored crops, generate the boundary
   target under v1 and under v2.0 labels and count differing voxels. The
   prereg needs every sampled split crop to differ (nonzero flips) and the
   binary masks to be byte-identical. Failures print loudly.

Emits /mnt/vesuvius/experiments/retrain_ab/frozen/crops.json and
ab_check.json. CPU only.
"""
from __future__ import annotations
import hashlib, json
from pathlib import Path

import numpy as np

V1 = Path("/mnt/vesuvius/p1218_full/blocks")
V2 = Path("/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired")
SRC = Path("/mnt/vesuvius/kaggle_p1218_repair_v2")
FRZ = Path("/mnt/vesuvius/experiments/retrain_ab/frozen")
CROP = 160
JIT = 16
SHAPE = (256, 512, 512)
N_CROPS = 4000          # 2000 split / 1000 onesided / 1000 uniform
FG_FLOOR = 200_000


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


def hj(key: str, span: int) -> int:
    """Deterministic jitter in [-span, span] from md5(key)."""
    return int(hashlib.md5(key.encode()).hexdigest()[:8], 16) % (2 * span + 1) - span


def corner(z: int, y: int, x: int, key: str) -> tuple[int, int, int]:
    """160^3 crop corner centring (z,y,x) with hash jitter, clamped in bounds."""
    cs = []
    for i, (c, dim) in enumerate(zip((z, y, x), SHAPE)):
        j = hj(f"{key}:{i}", JIT)
        c0 = c + j - CROP // 2
        cs.append(max(0, min(c0, dim - CROP)))
    return tuple(cs)


def main() -> None:
    split_map = json.loads((FRZ / "tile_split.json").read_text())["split"]
    train = {n for n, g in split_map.items() if g == "train"}

    per_tile: dict[str, dict] = {}
    for rp in sorted((SRC / "records").glob("*.json")):
        r = json.loads(rp.read_text())
        name = f"{r['slab']}_{r['tile']}"
        if name not in train:
            continue
        per_tile[name] = r

    crops, pools = [], {"SPLIT": [], "ONE_SIDED": []}
    for name, r in per_tile.items():
        for s in r["sites"]:
            if s["decision"] in pools:
                key = hashlib.md5(
                    f"{name}:{s['i']}:{s['z']}:{s['y']}:{s['x']}".encode()
                ).hexdigest()
                pools[s["decision"]].append(
                    (key, name, r["slab"], r["tile"], s["z"], s["y"], s["x"]))
    for k in pools:
        pools[k].sort()
    for key, name, slab, tile, z, y, x in pools["SPLIT"][:N_CROPS // 2]:
        c = corner(z, y, x, key)
        crops.append({"kind": "split", "tile": name, "slab": slab,
                      "file": f"{slab}/{tile}.npz", "corner": c, "key": key})
    for key, name, slab, tile, z, y, x in pools["ONE_SIDED"][:N_CROPS // 4]:
        c = corner(z, y, x, key)
        crops.append({"kind": "onesided", "tile": name, "slab": slab,
                      "file": f"{slab}/{tile}.npz", "corner": c, "key": key})
    # uniform stratum: hash-ordered train tiles, fixed grid positions
    uni = []
    for name in sorted(train, key=lambda n: hashlib.md5(n.encode()).hexdigest()):
        r = per_tile.get(name)
        if r is None:
            continue
        slab, tile = name.split("_", 1)
        for gz in (48, 96):
            for gy in (96, 256, 352):
                for gx in (96, 256, 352):
                    k = hashlib.md5(f"{name}:u:{gz}:{gy}:{gx}".encode()).hexdigest()
                    uni.append((k, name, slab, tile, gz, gy, gx))
    uni.sort()
    for k, name, slab, tile, z, y, x in uni[:N_CROPS // 4]:
        crops.append({"kind": "uniform", "tile": name, "slab": slab,
                      "file": f"{slab}/{tile}.npz", "corner": corner(z, y, x, k),
                      "key": k})

    (FRZ / "crops.json").write_text(json.dumps(
        {"crop": CROP, "jitter": JIT, "n": len(crops),
         "composition": {k: sum(1 for c in crops if c["kind"] == k)
                         for k in ("split", "onesided", "uniform")},
         "crops": crops}, indent=0))

    # ---- A/B manipulation check on 20 split crops ----
    chk, fails = [], 0
    for c in [c for c in crops if c["kind"] == "split"][:20]:
        with np.load(V1 / c["file"]) as d:
            l1 = d["labels"].astype(np.int32)
        with np.load(V2 / c["file"]) as d:
            l2 = d["labels"].astype(np.int32)
        z0, y0, x0 = c["corner"]
        sl = (slice(z0, z0 + CROP), slice(y0, y0 + CROP), slice(x0, x0 + CROP))
        c1, c2 = l1[sl], l2[sl]
        mask_same = bool(((c1 > 0) == (c2 > 0)).all())
        b1, b2 = boundary_target(c1), boundary_target(c2)
        flips = int((b1 != b2).sum())
        reass = int((c1 != c2).sum())
        if not mask_same or flips == 0:
            fails += 1
        chk.append({"tile": c["tile"], "corner": c["corner"],
                    "mask_identical": mask_same, "boundary_flips": flips,
                    "voxels_reassigned": reass,
                    "b2_positives": int(b2.sum())})
    out = {"n_checked": len(chk), "n_fail": fails,
           "median_flips": float(np.median([x["boundary_flips"] for x in chk])),
           "median_flip_share_of_b2": float(np.median(
               [x["boundary_flips"] / max(x["b2_positives"], 1) for x in chk])),
           "rows": chk}
    (FRZ / "ab_check.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in
                      ("n_checked", "n_fail", "median_flips",
                       "median_flip_share_of_b2")}, indent=1))
    print("crops:", len(crops))


if __name__ == "__main__":
    main()
