"""Figure: what the tight-contact bands actually look like in CT."""
import glob, json
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT="/mnt/vesuvius/tightgap1218"
bands=["0-2","2-4","4-6","6-10","10+"]
rows=json.load(open(OUT+"/sites_gaps.json"))
gapof={}
for f in glob.glob(OUT+"/crops/*.npz"):
    d=np.load(f,allow_pickle=True)
    gapof.setdefault(str(d["band"]),[]).append((float(d["gap"]),f))
fig,axes=plt.subplots(2,6,figsize=(16.5,6.0))
for j,b in enumerate(bands):
    items=sorted(gapof.get(b,[]))
    for i,ax in enumerate(axes[:,j]):
        if i<len(items):
            g,f=items[i*max(1,len(items)//2)] if len(items)>1 else items[0]
            v=np.load(f,allow_pickle=True)["intensity"]
            ax.imshow(v[64],cmap="gray"); ax.set_title("%s vox  (gap %.1f)"%(b,g),fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
items=[(0,f) for f in sorted(glob.glob(OUT+"/control/*.npz"))]
for i,ax in enumerate(axes[:,5]):
    if i<len(items):
        v=np.load(items[i*20][1],allow_pickle=True)["intensity"]
        ax.imshow(v[64],cmap="gray"); ax.set_title("control (single sheet)",fontsize=8)
    ax.set_xticks([]); ax.set_yticks([])
counts=json.load(open(OUT+"/crops_summary.json"))["realised"]
fig.suptitle("PHerc1218 tight-contact validation set: central slice of one crop per gap band\n"
  "realised crops  0-2:%d  2-4:%d  4-6:%d  6-10:%d  10+:%d  control:60   (gap measured along the label normal, level 0)"
  %(counts.get("0-2",0),counts.get("2-4",0),counts.get("4-6",0),counts.get("6-10",0),counts.get("10+",0)),fontsize=10)
fig.tight_layout(rect=[0,0,1,0.90])
fig.savefig(OUT+"/tightgap_bands.png",dpi=130)
print("saved")
