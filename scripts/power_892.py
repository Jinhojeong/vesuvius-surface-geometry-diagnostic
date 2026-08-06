"""Sensitivity of the 892 located-vs-nonlocated flag-rate comparison.

Resampling power: under the null-preserving resample of the pooled
per-volume flag rates into groups of 189 and 703, an additive shift
delta is applied to the located draw; power = fraction of simulations
with two-sided Mann-Whitney p < 0.05. Establishes what effect size the
reported p = 0.33 null could have detected."""
import json
import sys

import numpy as np
from scipy.stats import mannwhitneyu

SRC = sys.argv[1] if len(sys.argv) > 1 else "results/flag_split_892.json"
d = json.load(open(SRC))
rates = np.array([r["flag_rate"] for r in d["rows"]])
n_loc, n_non = 189, 703
rng = np.random.default_rng(1218)
out = {"n_sims": 500, "n_located": n_loc, "n_nonlocated": n_non,
       "pooled_median": float(np.median(rates)), "power": {}}
for delta_pp in (0.03, 0.06, 0.1, 0.2):
    delta = delta_pp / 100.0
    hits = 0
    for _ in range(out["n_sims"]):
        perm = rng.permutation(rates)
        a = perm[:n_loc] + delta
        b = perm[n_loc:n_loc + n_non]
        p = mannwhitneyu(a, b, alternative="two-sided").pvalue
        hits += p < 0.05
    out["power"][f"{delta_pp}pp"] = hits / out["n_sims"]
    print(f"delta {delta_pp}pp: power {hits / out['n_sims']:.3f}", flush=True)
json.dump(out, open("results/power_892.json", "w"), indent=1)
print("wrote results/power_892.json")
