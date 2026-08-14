"""Figure: what a support-aware stopping rule would remove, and where."""
import json, sys
import numpy as np, zarr, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy import ndimage as ndi
sys.path.insert(0,"/mnt/vesuvius/vcbuild"); import score_ab
ROOT="/mnt/vesuvius/p8_sprint"
rep=json.load(open(ROOT+"/r1_repair.json")); lo,hi=rep["bounds_l0"]
sup=zarr.open(ROOT+"/r1_support.zarr/0",mode="r")
rows=json.load(open(ROOT+"/r1_trim.json"))["rows"]
cmap=ListedColormap([(0.06,0.06,0.08),(0.85,0.25,0.2),(0.55,0.75,0.95)])
fig,axes=plt.subplots(1,3,figsize=(13.5,5.0))
for k,i in enumerate([0,1,2]):
    run=score_ab.load_tifxyz("%s/trace/r1_s%d_before"%(ROOT,i))
    v=run["valid"]; pts=np.stack([run["z"],run["y"],run["x"]],-1)
    inb=((pts>=np.array(lo))&(pts<np.array(hi))).all(-1)&v
    supm=np.zeros(v.shape,bool)
    for a,b in np.argwhere(inb):
        supm[a,b]= int(sup[int(round(pts[a,b,0])),int(round(pts[a,b,1])),int(round(pts[a,b,2]))])>0
    keep=ndi.binary_opening(supm,structure=np.ones((3,3)))
    img=np.zeros(v.shape,int); img[inb & ~keep]=1; img[inb & keep]=2
    r=[x for x in rows if x["seed"]==i and x["arm"]=="before"][0]
    ax=axes[k]; ax.imshow(img,cmap=cmap,vmin=0,vmax=2,interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("seed %d\nrule keeps %.0f%% of in-region quads\nremoved at edge-distance %.1f, kept %.1f"
                 %(i,100*r["retained_frac"],r["removed_mean_edge_dist"],r["kept_mean_edge_dist"]),fontsize=9)
fig.suptitle("PHerc1218: what a support-aware stopping rule would drop\n"
             "blue kept (CT under the quad), red dropped (no CT), patch grid space",fontsize=11)
fig.tight_layout(rect=[0,0,1,0.90])
fig.savefig(ROOT+"/support_rule.png",dpi=130); print("saved")
