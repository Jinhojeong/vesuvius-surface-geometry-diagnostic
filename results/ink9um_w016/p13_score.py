"""Does the w016 reading replicate on the other two validation segments?

Same endpoints and the same rule: the bar is that segment's own seed spread.
"""
import glob, json, sys
import numpy as np, tifffile
sys.path.insert(0,"/mnt/vesuvius")
import p12_score as S

W="/mnt/vesuvius/ink9um_w016"
STEPS=[10000,20000,30000,40000,50000,60000,75000]

def labels_for(seg):
    d=np.load("%s/%s_val_labels.npz"%(W,seg))
    z=int(np.argmax((d["validation"]>0).reshape(d["validation"].shape[0],-1).sum(1)))
    return d["validation"][z]>0, d["ink"][z]>0

def score_map(img, m, y):
    p=img[m].astype(np.float32); t=y[m]
    best=(-1,None)
    for thr in np.arange(0.05,0.96,0.01):
        pred=p>=thr
        tp=float((pred&t).sum()); fp=float((pred&~t).sum()); fn=float((~pred&t).sum())
        f1=0.0 if (2*tp+fp+fn)==0 else 2*tp/(2*tp+fp+fn)
        if f1>best[0]: best=(f1,float(thr))
    return dict(f1=round(best[0],5), thr=round(best[1],3), auc=round(S.auc(p,t),5),
                pos=round(float(t.mean()),5))

out={}
for seg,pre in (("pherc0139-w016",W+"/preds"),("pherc0814-46527",W+"/preds_pherc0814-46527"),
                ("pherc1667-w029",W+"/preds_pherc1667-w029")):
    if seg=="pherc0139-w016":
        m,y=S.load_labels()
    else:
        m,y=labels_for(seg)
    rows={}
    for sd in (42,43):
        for st in STEPS:
            f="%s/hybrid_3d2d-seed%d-step-%06d.tif"%(pre,sd,st)
            img=tifffile.imread(f).astype(np.float32)/255.0
            if img.shape!=m.shape: raise SystemExit("%s shape %s vs mask %s"%(seg,img.shape,m.shape))
            rows[(sd,st)]=score_map(img,m,y)
    ens=np.mean([tifffile.imread("%s/hybrid_3d2d-seed%d-step-%06d.tif"%(pre,sd,st)).astype(np.float32)/255.0
                 for sd in (42,43) for st in STEPS],axis=0)
    e=score_map(ens,m,y)
    spread=[abs(rows[(42,s)]["f1"]-rows[(43,s)]["f1"]) for s in STEPS]
    aspread=[abs(rows[(42,s)]["auc"]-rows[(43,s)]["auc"]) for s in STEPS]
    bf1=max(r["f1"] for r in rows.values()); bauc=max(r["auc"] for r in rows.values())
    bestk=max(rows,key=lambda k:rows[k]["f1"])
    bar=float(np.median(spread))
    out[seg]=dict(pos_rate=list(rows.values())[0]["pos"],
      best_single=dict(seed=bestk[0],step=bestk[1],f1=bf1,auc=rows[bestk]["auc"]),
      f1_range=[round(min(r["f1"] for r in rows.values()),4),round(bf1,4)],
      bar_f1=round(bar,5), bar_auc=round(float(np.median(aspread)),5),
      ensemble=dict(f1=e["f1"],auc=e["auc"],thr=e["thr"]),
      delta_f1=round(e["f1"]-bf1,5), delta_auc=round(e["auc"]-bauc,5),
      beats_all_f1=bool(e["f1"]>bf1), beats_all_auc=bool(e["auc"]>bauc),
      verdict="gain" if (e["f1"]-bf1)>bar else "indistinguishable")
    o=out[seg]
    print("%-16s ink %.3f | best single F1 %.4f (s%d/%dk) | ens F1 %.4f (%+.4f) | bar %.4f -> %-18s | beats all: F1 %s AUC %s"
      %(seg,o["pos_rate"],bf1,bestk[0],bestk[1]//1000,e["f1"],o["delta_f1"],bar,o["verdict"],
        o["beats_all_f1"],o["beats_all_auc"]))
json.dump(out,open(W+"/three_segment_scores.json","w"),indent=1)
print()
print("beats-all-members on F1: %d of 3 | on AUC: %d of 3"%(
  sum(v["beats_all_f1"] for v in out.values()), sum(v["beats_all_auc"] for v in out.values())))
print("preregistered-bar verdict 'gain': %d of 3"%sum(v["verdict"]=="gain" for v in out.values()))
