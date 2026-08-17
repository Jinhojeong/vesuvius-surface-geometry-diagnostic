#!/usr/bin/env python3
"""kollesis13, step 1: convert census13 clusters (test 1.3-1.6 and
background >=1.6) to arc position.

Inputs, all from this repo:
  results/census13/census13/          per-tile cluster jsons (the 1.3 pass)
  results/kollesis/origins_merged.csv per-slice mask centroids
  results/kollesis/positions_merged.csv.gz  crossings (rings frame)

Output: results/kollesis13/census13_arc.csv.gz + a small convert report.

Cross-check that pins this conversion: run the same code over the census8k
tables with these same inputs and all 491,337 rows of
results/kollesis/clusters_arc.csv.gz come back identical (zero drops, zero
ring changes, |ds| <= 5 um = the csv's own 2-decimal rounding). iyando's
pitch_qa_ray_positions.csv.gz gives identical results as the crossings
input; the merged file is used here so the chain runs from this repo alone.
"""
import gzip, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "scripts"))
import kollesis_lattice as K

C13 = os.path.join(REPO, "results", "census13", "census13")
ORIGINS = os.path.join(REPO, "results", "kollesis", "origins_merged.csv")
RAYS = os.path.join(REPO, "results", "kollesis", "positions_merged.csv.gz")
OUTDIR = os.path.join(REPO, "results", "kollesis13")
os.makedirs(OUTDIR, exist_ok=True)

origins = K.load_origins(ORIGINS)
rings, zs, ths = K.load_rings(RAYS)
tabs = K.ring_tables(rings, zs, ths)

rows = K.scan_census(C13)
print(f"[scan] {len(rows):,} clusters in census13 tiles", flush=True)
conv, drops = K.convert(rows, origins, rings, tabs, zs, ths)
print(f"[convert] {len(conv):,} converted; drops {drops}", flush=True)

with gzip.open(os.path.join(OUTDIR, "census13_arc.csv.gz"), "wt",
               newline="") as f:
    f.write("s_mm,z,ratio,th,n_sites,ring,dr_um,gz,gy,gx\n")
    for c in conv:
        f.write(f"{c['s']:.2f},{c['z']},{c['ratio']},{c['th']},"
                f"{c['n_sites']},{c['ring']},{c['dr_um']:.1f},"
                f"{c['gz']},{c['gy']},{c['gx']}\n")

r = np.array([c["ratio"] for c in conv])
test = (r >= 1.3) & (r < 1.6)
bg = r >= 1.6
rep = {"scanned": len(rows), "converted": len(conv), "drops": drops,
       "test_1.3_1.6": int(test.sum()), "background_ge_1.6": int(bg.sum()),
       "test_ratio_median": float(np.median(r[test])) if test.any() else None,
       "n_slices_used": len({c["z"] for c in conv})}
print(json.dumps(rep, indent=1))
json.dump(rep, open(os.path.join(OUTDIR, "census13_convert_report.json"),
                    "w"), indent=1)
