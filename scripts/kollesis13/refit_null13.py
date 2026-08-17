#!/usr/bin/env python3
"""kollesis13, step 3: refit-per-permutation nulls for the two low p's.

Per the running protocol: "If the label-permutation comes out positive, a
full refit-per-permutation null is required before any claim."

1) ring<=16, W=160: null with offset REFIT on each permutation (the step-2
   arm holds the offset fitted to the observed test -> upward bias).
2) all, W free: the fit drifts to the grid edge (120), outside 150-170 -
   fails the pre-registered check; a full refit null (W and offset, 200
   permutations) is run anyway to leave the number.

Also: offset consistency across cuts (a real lattice shares one phase).

Reads results/kollesis13/census13_arc.csv.gz; writes refit_null_report.json.
"""
import gzip, csv, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "scripts"))
import kollesis_lattice as K
from kollesis_lattice2 import fit_offset_only, on_frac

OUTDIR = os.path.join(REPO, "results", "kollesis13")

rows = []
with gzip.open(os.path.join(OUTDIR, "census13_arc.csv.gz"), "rt") as f:
    for r in csv.DictReader(f):
        rows.append((float(r["s_mm"]), float(r["ratio"]), int(r["ring"])))
s = np.array([r[0] for r in rows])
ratio = np.array([r[1] for r in rows])
ring = np.array([r[2] for r in rows])
isT_all = ratio < 1.6
rep = {}

# --- offsets per cut at W=160 (a real lattice shares one phase) ---
offs = {}
for name, keep in (("all", np.ones(len(s), bool)),
                   ("ring_le_38", ring <= 38),
                   ("ring_le_16", ring <= 16)):
    o, f = fit_offset_only(s[keep & isT_all], 160.0)
    offs[name] = o
rep["offsets_W160_by_cut"] = offs
print("[offsets W160]", offs, flush=True)

# --- 1) ring<=16, W=160, offset refit per permutation ---
keep = ring <= 16
sv, isT = s[keep], isT_all[keep]
nT = int(isT.sum())
oT, _ = fit_offset_only(sv[isT], 160.0)
obs = on_frac(sv[isT], 160.0, oT) - on_frac(sv[~isT], 160.0, oT)
rng = np.random.default_rng(20260810)
null = []
idx = np.arange(len(sv))
for i in range(2000):
    rng.shuffle(idx)
    t, b = idx[:nT], idx[nT:]
    o, _ = fit_offset_only(sv[t], 160.0)
    null.append(on_frac(sv[t], 160.0, o) - on_frac(sv[b], 160.0, o))
null = np.array(null)
p = float(np.mean(null >= obs))
rep["ring_le_16_W160_refit"] = {"obs_pp": obs * 100, "perm_p": p,
                                "null_mean_pp": float(null.mean() * 100),
                                "null_sd_pp": float(null.std() * 100),
                                "nperm": 2000}
print(f"[ring<=16 W160 refit-offset] obs {obs*100:+.2f}pp p={p:.3f} "
      f"(null {null.mean()*100:+.2f}+-{null.std()*100:.2f}pp)", flush=True)

# --- 2) all, W free, full refit per permutation (200) ---
sv, isT = s, isT_all
nT = int(isT.sum())
W0, o0, _ = K.fence_fit_fast(sv[isT])
obs = on_frac(sv[isT], W0, o0) - on_frac(sv[~isT], W0, o0)
null = []
idx = np.arange(len(sv))
for i in range(200):
    rng.shuffle(idx)
    t, b = idx[:nT], idx[nT:]
    W, o, _ = K.fence_fit_fast(sv[t])
    null.append(on_frac(sv[t], W, o) - on_frac(sv[b], W, o))
    if (i + 1) % 50 == 0:
        print(f"  [refit-full {i+1}/200]", flush=True)
null = np.array(null)
p = float(np.mean(null >= obs))
rep["all_Wfree_full_refit"] = {"W_mm": W0, "offset_mm": o0,
                               "obs_pp": obs * 100, "perm_p": p,
                               "null_mean_pp": float(null.mean() * 100),
                               "null_sd_pp": float(null.std() * 100),
                               "nperm": 200,
                               "note": "W=120 = grid edge, outside 150-170"}
print(f"[all Wfree refit-full] W={W0:.1f} obs {obs*100:+.2f}pp p={p:.3f} "
      f"(null {null.mean()*100:+.2f}+-{null.std()*100:.2f}pp)", flush=True)

json.dump(rep, open(os.path.join(OUTDIR, "refit_null_report.json"), "w"),
          indent=1)
