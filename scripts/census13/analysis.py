"""Census13 readout for Diego's two questions."""
import json, glob
import numpy as np

c13 = {}
for f in glob.glob("/mnt/vesuvius/census13/*/*.json"):
    r = json.load(open(f))
    c13[r["tile"]] = r
c8 = {}
for f in glob.glob("/mnt/vesuvius/census8k/*/*.json"):
    r = json.load(open(f))
    if r["tile"] in c13:
        c8[r["tile"]] = r

ratios13, ratios8 = [], []
n_thin = 0
for t, r in c13.items():
    if r.get("thin"):
        n_thin += 1
        continue
    ratios13 += [c["ratio"] for c in r["clusters"]]
    if t in c8 and not c8[t].get("thin"):
        ratios8 += [c["ratio"] for c in c8[t]["clusters"]]

r13 = np.array(ratios13)
r8 = np.array(ratios8)
band = r13[(r13 >= 1.3) & (r13 < 1.6)]
above13 = r13[r13 >= 1.6]
print(json.dumps({
    "tiles": len(c13), "thin": n_thin,
    "clusters_at_1.3": len(r13),
    "band_1.3_1.6": len(band),
    "at_or_above_1.6_in_13pass": len(above13),
    "census8k_same_tiles_clusters": len(r8),
    "band_share_of_all": round(len(band) / max(len(r13), 1), 4),
    "band_over_above": round(len(band) / max(len(above13), 1), 4),
    "band_ratio_median": round(float(np.median(band)), 3) if len(band) else None,
    "consistency_1.6_recount_vs_census8k": round(len(above13) / max(len(r8), 1), 4),
    "median_seconds": float(np.median([r["seconds"] for r in c13.values() if not r.get("thin")])),
}, indent=1))
# per-tile band density for the reply
per = []
for t, r in c13.items():
    if r.get("thin"):
        continue
    rs = np.array([c["ratio"] for c in r["clusters"]])
    per.append((len(rs[(rs >= 1.3) & (rs < 1.6)]), len(rs[rs >= 1.6])))
b, a = np.array([p[0] for p in per]), np.array([p[1] for p in per])
print("per-tile band median", float(np.median(b)), "above median", float(np.median(a)))
json.dump({"ratios13": r13.tolist()}, open("/mnt/vesuvius/census13_pooled.json", "w"))
