"""Arms 2 to 4: seed, step and full ensembles, computed from arm 1 probability maps."""
import glob, json, sys
import numpy as np, tifffile
sys.path.insert(0,"/mnt/vesuvius")
from p12_score import score

OUT="/mnt/vesuvius/ink9um_w016"
a1=json.load(open(OUT+"/arm1_scores.json"))
bar=a1["seed_spread_f1"]["median"]; best=a1["best_single"]
def load(seed,step):
    return tifffile.imread("%s/preds/hybrid_3d2d-seed%d-step-%06d.tif"%(OUT,seed,step)).astype(np.float32)/255.0
steps=[10000,20000,30000,40000,50000,60000,75000]
rows=[]

for st in steps:
    m=(load(42,st)+load(43,st))/2.0
    r=score(m,"seed-ens-step%06d"%st); r["arm"]=2; rows.append(r)
for sd in (42,43):
    m=np.mean([load(sd,s) for s in steps],axis=0)
    r=score(m,"step-ens-seed%d"%sd); r["arm"]=3; rows.append(r)
m=np.mean([load(sd,s) for sd in (42,43) for s in steps],axis=0)
r=score(m,"full-ens-14"); r["arm"]=4; rows.append(r)
# late-steps only, a cheap variant of arm 3
m=np.mean([load(sd,s) for sd in (42,43) for s in (50000,60000,75000)],axis=0)
r=score(m,"ens-late6"); r["arm"]=3; rows.append(r)

for r in sorted(rows,key=lambda x:-x["f1_best"]):
    d=r["f1_best"]-best["f1"]
    verdict="gain" if d>bar else ("indistinguishable" if d>-bar else "worse")
    print("%-22s f1 %.4f (thr %.2f) auc %.4f | vs best single %+0.4f -> %s"
          %(r["tag"],r["f1_best"],r["thr_best"],r["auc"],d,verdict))
best_ens=max(rows,key=lambda x:x["f1_best"])
out=dict(prereg="8ceda1e7a0536613", bar_f1=bar, best_single=best,
         best_ensemble=dict(tag=best_ens["tag"],f1=best_ens["f1_best"],auc=best_ens["auc"],thr=best_ens["thr_best"]),
         delta_f1=round(best_ens["f1_best"]-best["f1"],5),
         verdict=("gain" if best_ens["f1_best"]-best["f1"]>bar else "indistinguishable"),
         rows=rows)
json.dump(out,open(OUT+"/arms234_scores.json","w"),indent=1)
print()
print("PRIMARY: best single %.4f -> best ensemble %.4f, delta %+0.4f against bar %.4f => %s"
      %(best["f1"],best_ens["f1_best"],best_ens["f1_best"]-best["f1"],bar,out["verdict"]))
