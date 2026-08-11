import json
import numpy as np
d = "/mnt/vesuvius/experiments/wrapper_check/"
pr = {}
for ln in open(d + "PR1386_default.jsonl"):
    r = json.loads(ln)
    pr[r["sample"]] = r
rows = json.load(open(d + "wrapper_four_way.json"))["rows"]
diffs = {"C": [], "W": [], "P": []}
meds = {"new": [[], []], "C": [[], []], "W": [[], []]}
for r in rows:
    s = r["sample"]
    if s not in pr:
        continue
    new = pr[s].get("scores") or pr[s].get("result") or None
    if new is None:
        for k in pr[s]:
            if isinstance(pr[s][k], list) and len(pr[s][k]) == 3:
                new = pr[s][k]
                break
    for arm in diffs:
        a = r[arm]
        diffs[arm].append(max(abs(new[0] - a[0]), abs(new[1] - a[1])))
    for i in (0, 1):
        meds["new"][i].append(new[i])
        meds["C"][i].append(r["C"][i])
        meds["W"][i].append(r["W"][i])
print("n =", len(diffs["C"]))
print("medians  new [%.4f %.4f]  C [%.4f %.4f]  W [%.4f %.4f]" % (
    np.median(meds["new"][0]), np.median(meds["new"][1]),
    np.median(meds["C"][0]), np.median(meds["C"][1]),
    np.median(meds["W"][0]), np.median(meds["W"][1])))
for arm in ("C", "P", "W"):
    v = diffs[arm]
    print("max per-vol |delta(recall,prec)| vs %s: median %.5f  max %.5f" % (
        arm, np.median(v), max(v)))
