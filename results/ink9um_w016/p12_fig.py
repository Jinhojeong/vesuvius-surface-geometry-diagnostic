"""Figure: every arm against the preregistered bar, plus the prediction maps."""
import glob, json
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import tifffile

OUT="/mnt/vesuvius/ink9um_w016"
a1=json.load(open(OUT+"/arm1_scores.json")); a2=json.load(open(OUT+"/arms234_scores.json"))
a5=json.load(open(OUT+"/arms56_scores.json"))
ref=a5["reference"]; bar=a1["seed_spread_f1"]["median"]

fig=plt.figure(figsize=(14,7.6))
ax=fig.add_axes([0.30,0.56,0.66,0.36])
items=[("best single checkpoint",ref["f1_best"],"#4d79ff")]
items+= [(r["tag"].replace("full-ens-14","all 14 averaged").replace("step-ens-","step ens ").replace("ens-late6","6 late averaged"),
          r["f1_best"],"#31a354") for r in sorted(a2["rows"],key=lambda x:-x["f1_best"])[:3]]
items+= [(r["note"],r["f1_best"],"#e6842a") for r in a5["rows"] if r["note"].startswith("17 shifted")]
items+= [(r["note"],r["f1_best"],"#c0392b") for r in a5["rows"] if r["note"] in
         ("forward and reverse averaged","13 of 21, centred","reverse only")]
lbl=[i[0] for i in items]; val=[i[1] for i in items]; col=[i[2] for i in items]
y=np.arange(len(items))[::-1]
ax.barh(y,val,color=col,height=0.62)
ax.axvline(ref["f1_best"],color="#4d79ff",lw=1.2)
ax.axvspan(ref["f1_best"]-bar,ref["f1_best"]+bar,color="#4d79ff",alpha=0.12)
ax.set_yticks(y); ax.set_yticklabels(lbl,fontsize=9)
ax.set_xlim(0.30,0.82); ax.set_xlabel("F1 inside the w016 validation mask",fontsize=9)
ax.set_title("Inference-time arms against the preregistered bar (shaded = seed spread, %.3f F1)"%bar,fontsize=10)
for yy,v in zip(y,val): ax.text(v+0.004,yy,"%.3f"%v,va="center",fontsize=8)

d=np.load(OUT+"/w016_val_labels.npz"); zc=10
m=d["validation"][zc]>0
nzy,nzx=np.where(m); y0,y1,x0,x1=nzy.min(),nzy.max()+1,nzx.min(),nzx.max()+1
panels=[("labels",(d["ink"][zc][y0:y1,x0:x1]>0).astype(float),"gray"),
        ("best single",tifffile.imread(OUT+"/preds/hybrid_3d2d-seed42-step-020000.tif")[y0:y1,x0:x1]/255.0,"inferno"),
        ("all 14 averaged",np.mean([tifffile.imread(f)[y0:y1,x0:x1]/255.0 for f in sorted(glob.glob(OUT+"/preds/*.tif")) if not f.endswith("_reverse.tif")],axis=0),"inferno"),
        ("reverse only",tifffile.imread(OUT+"/preds56/d_rev.tif")[y0:y1,x0:x1]/255.0,"inferno")]
for i,(t,img,cm) in enumerate(panels):
    a=fig.add_axes([0.135,0.40-i*0.105,0.845,0.09])
    a.imshow(img,cmap=cm,aspect="auto",vmin=0,vmax=1)
    a.set_xticks([]); a.set_yticks([])
    a.set_ylabel(t,fontsize=8,rotation=0,ha="right",va="center",labelpad=8)
fig.text(0.5,0.955,"Released 9 um ink checkpoints on the w016 held-out mask",ha="center",fontsize=13,weight="bold")
fig.text(0.557,0.02,"maps cropped to the validation mask bounding box, 481 by 2152 at 9.6 um",ha="center",fontsize=8,color="#555")
fig.savefig(OUT+"/ink9um_arms.png",dpi=130,facecolor="white")
print("saved")
