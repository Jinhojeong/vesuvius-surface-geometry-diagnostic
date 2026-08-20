"""Figure: the offset distribution against the sign-even null."""
import json
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
D=json.load(open("/mnt/vesuvius/ridgecentre1218/ridge_offsets.json"))
rows=D["rows"]; s=D["summary"]
fig,axes=plt.subplots(1,3,figsize=(13.5,4.0))
for ax,h in zip(axes,(3,4,8)):
    v=np.array([r["off_%d"%h] for r in rows if "off_%d"%h in r])
    c=s["corridors"][str(h)]
    ax.hist(v,bins=np.arange(-h-0.125,h+0.3,0.25),color="#7f9fd6",edgecolor="white",linewidth=0.4)
    ax.axvline(0,color="#333",lw=1.2)
    ax.axvline(c["median"],color="#c0392b",lw=1.6,label="median %.2f"%c["median"])
    ax.axvspan(c["ci95"][0],c["ci95"][1],color="#c0392b",alpha=0.15,label="95%% CI [%.2f, %.2f]"%tuple(c["ci95"]))
    ax.set_title("corridor +/- %d vox\n%.1f%% negative, n=%d crops"%(h,100*c["frac_negative"],c["n"]),fontsize=9)
    ax.set_xlabel("ridge minus label-run centre, voxels\n(negative = ridge inward)",fontsize=8)
    ax.legend(fontsize=7,loc="upper right")
    ax.set_yticks([])
fig.suptitle("Do the repaired PHerc1218 labels sit on the CT ridge? Per-sheet offsets on the 209 published contact crops",
             fontsize=12,weight="bold")
fig.tight_layout(rect=[0,0,1,0.90])
fig.savefig("/mnt/vesuvius/ridgecentre1218/ridge_offsets.png",dpi=130,facecolor="white")
print("saved")
