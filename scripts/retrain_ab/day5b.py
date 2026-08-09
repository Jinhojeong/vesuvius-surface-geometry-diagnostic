"""Prereg build, day 5: was_mega stratum join + covariate table + 059 note.

1. was_mega: join the frozen primary sites against the v1.4-era agreement table
   on (slab, tile, local coords within 3 vox). The agreement rows are the
   54,377 solved sites of the SAMPLED era; the census primary sites are a
   different, denser enumeration, so the join is expected to be partial. If the
   matched stratum is under 500 the prereg reports it descriptively.

2. Covariate table: primary vs onesided distributions on th and ratio, the two
   covariates both lists carry. Quantiles plus a nearest-neighbour matched
   onesided subset (1:1 on standardized th/ratio, greedy, deterministic order)
   so the control can be read at matched covariates as well as raw.

3. Dataset059: different scroll (Scroll 1 patches vs PHerc1218), intersection
   empty by construction; recorded here so the prereg drop is a written fact.
"""
import csv, json
import numpy as np
from pathlib import Path

AB = Path("/mnt/vesuvius/experiments/retrain_ab")
FRZ = AB / "frozen"

prim = json.loads((FRZ / "sites_primary.json").read_text())["sites"]
ones = json.loads((FRZ / "sites_onesided.json").read_text())["sites"]

# ---- 1. was_mega join ----
agr = {}
with open("/mnt/vesuvius/kaggle_p1218_repair_v2/agreement_54377.csv") as f:
    for row in csv.DictReader(f):
        key = (row["slab"], row["tile"])
        agr.setdefault(key, []).append(
            (int(row["lz"]), int(row["ly"]), int(row["lx"]),
             row["was_mega"], row["recast_ok"]))
matched, mega = 0, 0
for s in prim:
    key = (s["slab"], s["tile"].split("_", 1)[1])
    best = None
    for lz, ly, lx, wm, rc in agr.get(key, []):
        d = max(abs(lz - s["z"]), abs(ly - s["y"]), abs(lx - s["x"]))
        if d <= 3 and (best is None or d < best[0]):
            best = (d, wm)
    if best is not None:
        matched += 1
        if best[1] == "1":
            mega += 1
out1 = {"primary_n": len(prim), "agreement_matched": matched,
        "was_mega_1": mega,
        "verdict": ("inferential" if matched >= 500 else "descriptive_only")}

# ---- 2. covariates ----
def stats(v):
    v = np.array(v, float)
    return {"q10": round(float(np.quantile(v, .1)), 3),
            "median": round(float(np.median(v)), 3),
            "q90": round(float(np.quantile(v, .9)), 3)}

pt = [s["th"] for s in prim]; pr = [s["ratio"] for s in prim]
ot = [s["th"] for s in ones]; orr = [s["ratio"] for s in ones]
mt, st_ = np.mean(pt + ot), np.std(pt + ot) + 1e-9
mr, sr = np.mean(pr + orr), np.std(pr + orr) + 1e-9
P = np.array([[(t - mt) / st_, (r - mr) / sr] for t, r in zip(pt, pr)])
O = np.array([[(t - mt) / st_, (r - mr) / sr] for t, r in zip(ot, orr)])
used = np.zeros(len(O), bool)
pairs = []
for i in range(min(len(P), 2000)):
    d = np.linalg.norm(O - P[i], axis=1)
    d[used] = 1e9
    j = int(np.argmin(d))
    if d[j] < 0.5:
        used[j] = True
        pairs.append((i, j, float(d[j])))
out2 = {"primary": {"th": stats(pt), "ratio": stats(pr)},
        "onesided": {"th": stats(ot), "ratio": stats(orr)},
        "matched_pairs_within_0.5sd": len(pairs),
        "matched_onesided_keys": [ones[j]["key"] for _, j, _ in pairs]}

out = {"was_mega": out1, "covariates": out2,
       "dataset059": "different scroll (Scroll 1 vs PHerc1218); intersection "
                     "empty by construction; concordance endpoint dropped"}
(FRZ / "day5_strata.json").write_text(json.dumps(out, indent=1))
print(json.dumps({"was_mega": out1,
                  "cov_primary": out2["primary"],
                  "cov_onesided": out2["onesided"],
                  "matched_pairs": len(pairs)}, indent=1))
