#!/usr/bin/env python3
"""kollesis13, step 4: diagnose the W-free edge residue + sensitivity floor.

a) Contrast sweep over W 60-200 (offset fitted on test per W): a real
   periodicity peaks at W* (and harmonics); a broad density non-uniformity
   gives drifting signal toward small W with no stable peak.
b) KS tests vs background over s (do the arc distributions differ at all,
   lattice or not?), and the same over ring and z.
c) Sensitivity at W=160 for 1-4 clusters per joint (step 2 already gives
   power 1.00 at 5).

Reads results/kollesis13/census13_arc.csv.gz; writes diagnose_edge.json.
"""
import gzip, csv, json, os, sys
import numpy as np
from scipy.stats import ks_2samp

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "scripts"))
from kollesis_lattice2 import fit_offset_only, on_frac

OUTDIR = os.path.join(REPO, "results", "kollesis13")

rows = []
with gzip.open(os.path.join(OUTDIR, "census13_arc.csv.gz"), "rt") as f:
    for r in csv.DictReader(f):
        rows.append((float(r["s_mm"]), float(r["ratio"]), int(r["ring"]),
                     int(r["z"])))
s = np.array([x[0] for x in rows])
ratio = np.array([x[1] for x in rows])
ring = np.array([x[2] for x in rows])
zz = np.array([x[3] for x in rows])
isT = ratio < 1.6
sT, sB = s[isT], s[~isT]

# a) sweep
ws, cs = [], []
for W in np.arange(60.0, 200.5, 2.0):
    o, _ = fit_offset_only(sT, W)
    cs.append((on_frac(sT, W, o) - on_frac(sB, W, o)) * 100)
    ws.append(W)
cs = np.array(cs); ws = np.array(ws)
top = ws[np.argsort(cs)[-5:]][::-1]
print(f"[scan] contrast max {cs.max():+.2f}pp at W={ws[np.argmax(cs)]:.0f}; "
      f"top-5 W: {sorted(top.tolist())}", flush=True)

# b) KS
ks_s = ks_2samp(sT, sB)
ks_r = ks_2samp(ring[isT], ring[~isT])
ks_z = ks_2samp(zz[isT], zz[~isT])
print(f"[KS] s: D={ks_s.statistic:.4f} p={ks_s.pvalue:.2e} | "
      f"ring: D={ks_r.statistic:.4f} p={ks_r.pvalue:.2e} | "
      f"z: D={ks_z.statistic:.4f} p={ks_z.pvalue:.2e}", flush=True)

# c) fine sensitivity
joints = np.arange(s.min() + 80.0, s.max(), 160.0)
sens = {}
for per_joint in (1, 2, 3, 4):
    n_inj = int(per_joint * len(joints))
    hits = 0
    for trial in range(20):
        rng = np.random.default_rng(3000 + trial)
        sm = sT.copy()
        pick = rng.choice(len(sm), size=n_inj, replace=False)
        sm[pick] = joints[rng.integers(0, len(joints), size=n_inj)] + \
            rng.normal(0, 0.12 * 160 / 2, size=n_inj)
        o, _ = fit_offset_only(sm, 160.0)
        on = np.abs(((np.concatenate([sm, sB]) - o + 80) % 160) - 80) <= 19.2
        nT = len(sm)
        obs = on[:nT].mean() - on[nT:].mean()
        idx = np.arange(len(on)); null = []
        for _ in range(500):
            rng.shuffle(idx)
            null.append(on[idx[:nT]].mean() - on[idx[nT:]].mean())
        hits += (np.mean(np.array(null) >= obs) < 0.05)
    sens[per_joint] = hits / 20
    print(f"  {per_joint}/joint ({n_inj} of {len(sT):,}): "
          f"power {hits/20:.2f}", flush=True)

json.dump({"w_scan": {"W": ws.tolist(), "contrast_pp": cs.tolist()},
           "ks": {"s": [ks_s.statistic, ks_s.pvalue],
                  "ring": [ks_r.statistic, ks_r.pvalue],
                  "z": [ks_z.statistic, ks_z.pvalue]},
           "sensitivity_fine": sens},
          open(os.path.join(OUTDIR, "diagnose_edge.json"), "w"), indent=1)
print("[done]", flush=True)
