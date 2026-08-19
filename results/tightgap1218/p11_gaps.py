"""P11 step 2: measure the normal-direction gap at every ray-validated split site.

Follows PREREGISTRATION.md (frozen 9d5d71dbaf45ab85), rule 2: for each split
site, the gap is measured along the label normal between the two instance ids
the split separated (A and B), at the site voxel, in level-0 voxels.

Reads the repair records only. Writes one row per site with its measured gap
and band, plus the diagnostics the prereg promises to report (slab spread,
instance-pair repetition). No crops here.
"""
import collections
import glob
import json
import os

import numpy as np

REPAIR = "/mnt/vesuvius/kaggle_p1218_repair_v2/repairs"
BLOCKS = "/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired"
OUT = "/mnt/vesuvius/tightgap1218"
os.makedirs(OUT, exist_ok=True)

BANDS = [(0, 2), (2, 4), (4, 6), (6, 10), (10, 10**9)]


def band_of(g):
    for lo, hi in BANDS:
        if lo <= g < hi:
            return "%g-%g" % (lo, hi) if hi < 10**9 else "10+"
    return None


def gap_along_normal(lab, p, a, b, maxr=32):
    """Distance in voxels from the A/B interface, walking both ways along the
    local normal estimated from the A-B separation direction."""
    z, y, x = p
    Z, Y, X = lab.shape
    if not (0 <= z < Z and 0 <= y < Y and 0 <= x < X):
        return None
    # local window
    r = 12
    z0, z1 = max(0, z - r), min(Z, z + r + 1)
    y0, y1 = max(0, y - r), min(Y, y + r + 1)
    x0, x1 = max(0, x - r), min(X, x + r + 1)
    w = lab[z0:z1, y0:y1, x0:x1]
    ma = np.argwhere(w == a)
    mb = np.argwhere(w == b)
    if len(ma) < 3 or len(mb) < 3:
        return None
    ca, cb = ma.mean(0), mb.mean(0)
    d = cb - ca
    n = np.linalg.norm(d)
    if n < 1e-6:
        return None
    d = d / n
    # walk from the site along +/- d, record where each label is first seen
    hit_a = hit_b = None
    for t in np.arange(0.0, maxr, 0.5):
        for sgn in (1, -1):
            q = np.round(np.array([z, y, x]) + sgn * t * d).astype(int)
            if not (0 <= q[0] < Z and 0 <= q[1] < Y and 0 <= q[2] < X):
                continue
            v = int(lab[q[0], q[1], q[2]])
            if v == a and hit_a is None:
                hit_a = t
            elif v == b and hit_b is None:
                hit_b = t
        if hit_a is not None and hit_b is not None:
            break
    if hit_a is None or hit_b is None:
        return None
    return float(hit_a + hit_b)


rows = []
files = sorted(glob.glob(REPAIR + "/*.json"))
print("repair files:", len(files), flush=True)
cache_key, cache_lab = None, None
for k, f in enumerate(files):
    d = json.load(open(f))
    reps = d.get("repairs") or []
    if not reps:
        continue
    slab, tile = d["slab"], d["tile"]
    npz = "%s/%s/%s.npz" % (BLOCKS, slab, tile)
    if not os.path.exists(npz):
        continue
    if cache_key != npz:
        try:
            z = np.load(npz)
            cache_lab = z[z.files[0]]
        except Exception:
            continue
        cache_key = npz
    lab = cache_lab
    for r in reps:
        p = r.get("site")
        a, b = r.get("A"), r.get("B")
        if not p or a is None or b is None:
            continue
        g = gap_along_normal(lab, p, int(a), int(b))
        if g is None:
            continue
        rows.append(dict(slab=slab, tile=tile, z=p[0], y=p[1], x=p[2],
                         A=int(a), B=int(b), gap=round(g, 2), band=band_of(g),
                         conf=r.get("conf"), tier=r.get("tier")))
    if (k + 1) % 300 == 0:
        print("  %d/%d files, %d sites" % (k + 1, len(files), len(rows)), flush=True)

by_band = collections.Counter(r["band"] for r in rows)
by_slab = collections.Counter(r["slab"] for r in rows)
pair_rep = collections.Counter((r["slab"], r["A"], r["B"]) for r in rows)
out = dict(prereg="9d5d71dbaf45ab85", n_sites=len(rows),
           by_band=dict(by_band), n_slabs=len(by_slab),
           top_slab_share=round(max(by_slab.values()) / max(1, len(rows)), 4),
           most_repeated_pair=int(max(pair_rep.values()) if pair_rep else 0))
json.dump(out, open(OUT + "/gaps_summary.json", "w"), indent=1)
json.dump(rows, open(OUT + "/sites_gaps.json", "w"))
print(json.dumps(out, indent=1))
