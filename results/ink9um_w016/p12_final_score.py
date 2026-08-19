"""Final scoring of arms 5 and 6 with the design errors corrected."""
import glob, json, sys
import numpy as np, tifffile
sys.path.insert(0,"/mnt/vesuvius")
from p12_score import score

OUT="/mnt/vesuvius/ink9um_w016"; P=OUT+"/preds56"
a1=json.load(open(OUT+"/arm1_scores.json"))
bar=a1["seed_spread_f1"]["median"]; abar=a1["seed_spread_auc"]["median"]
ref=[r for r in a1["rows"] if r["tag"]=="hybrid_3d2d-seed42-step-020000"][0]
L=lambda p: tifffile.imread(p).astype(np.float32)/255.0

print("reference (default window 2-18, forward): F1 %.4f AUC %.4f | bars F1 %.4f AUC %.4f\n"
      %(ref["f1_best"],ref["auc"],bar,abar))
rows=[]
def add(tag, img, note=""):
    r=score(img,tag); r["note"]=note; rows.append(r)
    d=r["f1_best"]-ref["f1_best"]; da=r["auc"]-ref["auc"]
    v="gain" if d>bar else ("worse" if d<-bar else "indistinguishable")
    print("%-10s f1 %.4f auc %.4f | dF1 %+0.4f dAUC %+0.4f -> %-18s %s"
          %(tag,r["f1_best"],r["auc"],d,da,v,note))

# arm 5, depth window: narrower (centred) and shifted at the same width
for tag,note in (("w_c13","13 of 21, centred"),("w_c9","9 centred"),("w_c5","5 centred"),
                 ("s_up0","17 shifted up, 0-16"),("s_dn4","17 shifted down, 4-20")):
    add(tag, L("%s/%s.tif"%(P,tag)), note)
# arm 6, direction
add("d_rev", L(P+"/d_rev.tif"), "reverse only")
fused=(L(P+"/d_fwd.tif")+L(P+"/d_rev.tif"))/2.0
add("d_fused", fused, "forward and reverse averaged")
mx=np.maximum(L(P+"/d_fwd.tif"), L(P+"/d_rev.tif"))
add("d_max", mx, "elementwise max of the two directions")

best=max(rows,key=lambda r:r["f1_best"])
out=dict(prereg="8ceda1e7a0536613", reference=ref, bar_f1=bar, bar_auc=abar,
  design_corrections=[
   "w_full and w_c17 were the same window; the checkpoint fixes 17 of 21 (indices 2-18)",
   "--layer-start/--layer-end can only narrow or shift, not widen past 17",
   "--direction both writes forward and reverse as two files rather than fusing them"],
  rows=rows, best=dict(tag=best["tag"],f1=best["f1_best"],auc=best["auc"]),
  verdict=("gain" if best["f1_best"]-ref["f1_best"]>bar else "indistinguishable"))
json.dump(out,open(OUT+"/arms56_scores.json","w"),indent=1)
print("\nbest of arms 5-6: %s F1 %.4f (%+0.4f vs reference) -> %s"
      %(best["tag"],best["f1_best"],best["f1_best"]-ref["f1_best"],out["verdict"]))
