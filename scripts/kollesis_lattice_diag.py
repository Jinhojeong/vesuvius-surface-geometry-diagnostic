"""Arc-axis sanity + ring-restricted contrast arm."""
import gzip, json
import numpy as np
import sys
sys.path.insert(0, "/mnt/vesuvius")
from kollesis_lattice import fence_fit_fast

rows = []
with gzip.open("/mnt/vesuvius/kollesis/clusters_arc.csv.gz", "rt") as f:
    hdr = f.readline()
    for line in f:
        p = line.split(",")
        rows.append((float(p[0]), int(p[5]), float(p[2])))  # s, ring, ratio
s = np.array([r[0] for r in rows])
ring = np.array([r[1] for r in rows])
ratio = np.array([r[2] for r in rows])
print(f"n={len(s)} s range {s.min():.0f}-{s.max():.0f}mm")
for q in (50, 90, 95, 99):
    print(f"  s p{q}: {np.percentile(s, q):.0f}mm, ring p{q}: {np.percentile(ring, q):.0f}")
print("ring hist:", {f"<=16": int((ring<=16).sum()), "17-38": int(((ring>16)&(ring<=38)).sum()), ">38": int((ring>38).sum())})

# ring perimeter sanity from the merged table
tabs = json.loads("{}")
# recompute per-ring perimeter stats via kollesis_lattice tables
from kollesis_lattice import load_origins, load_rings, ring_tables
origins = load_origins("/mnt/vesuvius/kollesis/origins_merged.csv")
rings_d, zs, ths = load_rings("/mnt/vesuvius/kollesis/positions_merged.csv.gz")
T = ring_tables(rings_d, zs, ths)
pers_by_ring = {}
for z, t in T.items():
    for j, per in enumerate(t["pers"]):
        pers_by_ring.setdefault(j, []).append(per)
for j in (0, 10, 20, 30, 38, 45, 60, 80):
    v = pers_by_ring.get(j)
    if v:
        print(f"ring {j}: n_slices {len(v)} perimeter median {np.median(v):.1f}mm range {min(v):.1f}-{max(v):.1f}")

def contrast(mask_keep, label):
    ss, bb = s[mask_keep], (ratio[mask_keep] >= 1.6) & (ratio[mask_keep] <= 2.6)
    sA, sB = ss[bb], ss[~bb]
    rng = np.random.default_rng(20260809)
    for W, name in ((None, "fitted"), (160.0, "prior160")):
        if W is None:
            W, o, _ = fence_fit_fast(sA)
        else:
            offs = np.arange(0, W, 2.0)
            r = np.abs(((sA[:, None] - offs[None, :] + W/2) % W) - W/2)
            o = float(offs[np.argmax((r <= 0.12*W).mean(axis=0))])
        def onf(x):
            return float((np.abs(((x - o + W/2) % W) - W/2) <= 0.12*W).mean())
        obs = onf(sA) - onf(sB)
        lab = bb.copy()
        null = np.empty(400)
        for i in range(400):
            rng.shuffle(lab)
            null[i] = onf(ss[lab]) - onf(ss[~lab])
        print(f"[{label}/{name}] W={W:.1f} o={o:.1f} n_band={bb.sum()} contrast {100*obs:+.3f}pp p={float((null>=obs).mean()):.4f} sd {100*null.std():.3f}pp")

contrast(ring <= 38, "ring<=38")
contrast(ring <= 16, "ring<=16")
