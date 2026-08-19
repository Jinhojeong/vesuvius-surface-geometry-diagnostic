"""Score arm 1 and report the seed spread that sets the bar for later arms."""
import glob, json, sys
import numpy as np, tifffile
sys.path.insert(0,"/mnt/vesuvius")
from p12_score import score

OUT="/mnt/vesuvius/ink9um_w016"
rows=[]
for f in sorted(glob.glob(OUT+"/preds/*.tif")):
    n=f.split("/")[-1][:-4]
    if n.endswith("_reverse"): continue
    img=tifffile.imread(f).astype(np.float32)/255.0
    r=score(img,n); r["seed"]=42 if "seed42" in n else 43
    r["step"]=int(n.split("step-")[1])
    rows.append(r); print("%-34s f1 %.4f (thr %.2f) auc %.4f call %.3f"%(n,r["f1_best"],r["thr_best"],r["auc"],r["calling_rate_at_best"]))

by=lambda s: {r["step"]:r for r in rows if r["seed"]==s}
a,b=by(42),by(43)
spread=[abs(a[k]["f1_best"]-b[k]["f1_best"]) for k in sorted(set(a)&set(b))]
aspread=[abs(a[k]["auc"]-b[k]["auc"]) for k in sorted(set(a)&set(b))]
best=max(rows,key=lambda r:r["f1_best"])
summ=dict(prereg="8ceda1e7a0536613", n_runs=len(rows),
  pos_rate=rows[0]["pos_rate"],
  seed_spread_f1=dict(median=round(float(np.median(spread)),5), max=round(float(np.max(spread)),5)),
  seed_spread_auc=dict(median=round(float(np.median(aspread)),5), max=round(float(np.max(aspread)),5)),
  best_single=dict(tag=best["tag"], f1=best["f1_best"], thr=best["thr_best"], auc=best["auc"]),
  f1_range=[round(min(r["f1_best"] for r in rows),5), round(max(r["f1_best"] for r in rows),5)],
  auc_range=[round(min(r["auc"] for r in rows),5), round(max(r["auc"] for r in rows),5)],
  rows=rows)
json.dump(summ,open(OUT+"/arm1_scores.json","w"),indent=1)
print()
print("seed spread F1 median %.4f max %.4f  <-- the bar for calling an ensemble a gain"%(np.median(spread),np.max(spread)))
print("F1 range across 14: %.4f to %.4f | AUC range %.4f to %.4f"%(summ["f1_range"][0],summ["f1_range"][1],summ["auc_range"][0],summ["auc_range"][1]))
print("best single:", best["tag"], "F1 %.4f AUC %.4f"%(best["f1_best"],best["auc"]))
print("positive rate (trivial floor reference): %.4f"%rows[0]["pos_rate"])
