"""Flag precision/recall on the finite-thickness twin, against true labels.

Same instrument as everywhere else (EDT-SDF normals, SPAN 12 / STEP 0.5,
contiguous run through p0, flag > 1.6x median thickness). Truth per sampled
site = number of distinct turn ids inside the run: >=2 means a real
multi-sheet stack. First ground-truthed precision number for the flag.
"""
import json
import sys

import numpy as np
from scipy import ndimage as ndi

RUN = sys.argv[1] if len(sys.argv) > 1 else "/mnt/vesuvius/twin/run1"
N_SAMPLE = 20000
SPAN, STEP = 12.0, 0.5
K = int(2 * SPAN / STEP) + 1
OFFS = np.linspace(-SPAN, SPAN, K)

lab = np.load(f"{RUN}/labels.npy")
mask = lab > 0
sm = ndi.gaussian_filter(
    ndi.distance_transform_edt(~mask).astype(np.float32)
    - ndi.distance_transform_edt(mask).astype(np.float32), 1.0)
g = np.stack(np.gradient(sm), 0).astype(np.float32)

surf = mask & ~ndi.binary_erosion(mask)
pts = np.argwhere(surf)
rng = np.random.default_rng(1218)
pts = pts[rng.choice(len(pts), size=min(N_SAMPLE, len(pts)), replace=False)]
nv = g[:, pts[:, 0], pts[:, 1], pts[:, 2]].T
nv /= np.linalg.norm(nv, axis=1, keepdims=True) + 1e-6
cc = pts[:, None, :].astype(np.float32) + OFFS[None, :, None] * nv[:, None, :]
onmask = ndi.map_coordinates(mask.astype(np.float32), cc.reshape(-1, 3).T,
                             order=0, mode="constant").reshape(len(pts), K) > 0.5
ctr = K // 2
shp = np.array(mask.shape)

rows = []
for i in range(len(pts)):
    if not onmask[i, ctr]:
        continue
    a = ctr
    while a > 0 and onmask[i, a - 1]:
        a -= 1
    b = ctr
    while b < K - 1 and onmask[i, b + 1]:
        b += 1
    th = (b - a) * STEP
    q = np.round(cc[i, a:b + 1]).astype(int)
    ok = ((q >= 0) & (q < shp)).all(1)
    ids = set(int(v) for v in lab[q[ok, 0], q[ok, 1], q[ok, 2]] if v > 0)
    rows.append((th, len(ids)))

th = np.array([r[0] for r in rows])
nid = np.array([r[1] for r in rows])
med = float(np.median(th))
flag = th > 1.6 * med
ctrl = (th >= 0.8 * med) & (th <= 1.2 * med)
multi = nid >= 2

def rate(x, y):
    return 100.0 * (x & y).sum() / max(x.sum(), 1)

print(f"valid sites {len(rows):,}, median thickness {med:.1f} vox-units")
print(f"flagged {flag.sum():,} ({100*flag.mean():.2f}% of sampled)")
print(f"  PRECISION  P(>=2 true sheets | flagged) = {rate(flag, multi):.1f}%")
print(f"  FP breakdown: flagged & single-sheet = {int((flag & ~multi).sum())}")
print(f"true multi-sheet sites in sample: {multi.sum():,} "
      f"({100*multi.mean():.2f}%)")
print(f"  RECALL     P(flagged | >=2 true sheets) = {rate(multi, flag):.1f}%")
print(f"control sites {ctrl.sum():,}")
print(f"  P(>=2 true sheets | control) = {rate(ctrl, multi):.1f}%")
json.dump({
    "n_valid": int(len(rows)), "median": med,
    "n_flagged": int(flag.sum()), "precision_pct": rate(flag, multi),
    "n_multi": int(multi.sum()), "recall_pct": rate(multi, flag),
    "n_control": int(ctrl.sum()), "control_multi_pct": rate(ctrl, multi),
}, open(f"{RUN}/flag_verdict.json", "w"), indent=1)
