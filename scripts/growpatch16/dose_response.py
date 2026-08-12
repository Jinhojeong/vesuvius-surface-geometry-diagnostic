import json
import numpy as np
from scipy.stats import spearmanr

R = json.load(open("/mnt/vesuvius/vcbuild/demo_out/powered16/RESULTS.json"))
cells = R["cells"]
orig = json.load(open("/mnt/vesuvius/hazard_zarr_smoke/demo_sites.json"))
new = R["new_site_selection"]["sites"]

meta = []
for s in orig:
    meta.append({"frag": s["slab"] + "_" + s["tile"], "n_sites": s["n_sites"],
                 "ratio": s["ratio"], "block": "orig"})
for s in new:
    if "/" in s["tile"]:
        frag = s["tile"].replace("/", "_")
    elif "slab" in s:
        frag = s["slab"] + "_" + s["tile"]
    else:
        frag = s["tile"]
    meta.append({"frag": frag, "n_sites": s["n_sites"], "ratio": s["ratio"],
                 "block": "new"})

rows = []
for ck, cv in cells.items():
    frag = ck.split("_", 1)[1]
    m = next((x for x in meta if x["frag"] == frag), None)
    if m is None:
        print("NO MATCH for", ck, frag)
        continue
    aA = cv["B"]["area"]["mean"] - cv["A"]["area"]["mean"]
    dA = cv["B"]["dbl"]["mean"] - cv["A"]["dbl"]["mean"]
    dC = cv["B"]["dbl"]["mean"] - cv["C"]["dbl"]["mean"]
    rows.append((ck, m["n_sites"], m["ratio"], m["block"], aA, dA, dC))

rows.sort(key=lambda r: -r[1])
hdr = ("site", "n_sites", "ratio", "blk", "dArea(B-A)", "dDbl(B-A)", "dDbl(B-C)")
print("%-38s %7s %6s %5s %12s %10s %10s" % hdr)
for r in rows:
    print("%-38s %7d %6.2f %5s %12.0f %10.2f %10.2f" % r)

ns = np.array([r[1] for r in rows], float)
rr = np.array([r[2] for r in rows], float)
aa = np.array([r[4] for r in rows])
dd = np.array([r[5] for r in rows])
dc = np.array([r[6] for r in rows])
print()
for name, y in (("area B-A", aa), ("dbl B-A", dd), ("dbl B-C", dc)):
    rho, p = spearmanr(ns, y)
    print("Spearman n_sites vs %-9s rho %+.3f p %.4f" % (name, rho, p))
for name, y in (("area B-A", aa), ("dbl B-A", dd)):
    rho, p = spearmanr(rr, y)
    print("Spearman ratio   vs %-9s rho %+.3f p %.4f" % (name, rho, p))

big = ns >= 300
print()
print("n_sites>=300 block: n", int(big.sum()),
      "| mean dArea %.0f | mean dDbl %.2f" % (aa[big].mean(), dd[big].mean()))
print("n_sites<300 block:  n", int((~big).sum()),
      "| mean dArea %.0f | mean dDbl %.2f" % (aa[~big].mean(), dd[~big].mean()))
