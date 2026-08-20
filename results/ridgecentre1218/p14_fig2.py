"""Figure: offsets on the CORRECT CT (level 1, the label's own grid), against the
shipped-crop version that measured label against a displaced CT."""
import json
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
A=json.load(open("/mnt/vesuvius/ridgecentre1218/ridge_offsets_L1.json"))
B=json.load(open("/mnt/vesuvius/ridgecentre1218/ridge_offsets.json"))
fig,axes=plt.subplots(1,3,figsize=(13.5,4.2))
for ax,h in zip(axes,(3,4,8)):
    va=np.array([r["off_%d"%h] for r in A["rows"] if "off_%d"%h in r])
    vb=np.array([r["off_%d"%h] for r in B["rows"] if "off_%d"%h in r])
    ca=A["summary"]["corridors"][str(h)]
    bins=np.arange(-h-0.125,h+0.3,0.25)
    ax.hist(vb,bins=bins,color="#cccccc",edgecolor="white",linewidth=0.3,
            label="withdrawn, displaced CT")
    ax.hist(va,bins=bins,color="#7f9fd6",edgecolor="white",linewidth=0.4,alpha=0.9,
            label="CT level 1, the label grid")
    ax.axvline(0,color="#333",lw=1.2)
    ax.axvline(ca["median"],color="#c0392b",lw=1.6)
    ax.axvspan(ca["ci95"][0]-0.02,ca["ci95"][1]+0.02,color="#c0392b",alpha=0.18)
    ax.set_title("corridor +/- %d vox\nmedian %.2f, 95%% CI [%.2f, %.2f], %.1f%% negative"
                 %(h,ca["median"],ca["ci95"][0],ca["ci95"][1],100*ca["frac_negative"]),fontsize=9)
    ax.set_xlabel("ridge minus label-run centre, voxels\n(negative = ridge inward)",fontsize=8)
    ax.legend(fontsize=7,loc="upper right")
    ax.set_yticks([])
fig.suptitle("Do the repaired PHerc1218 labels sit on the CT ridge? 209 contact crops, measured on the label's own grid",
             fontsize=12,weight="bold")
fig.tight_layout(rect=[0,0,1,0.90])
fig.savefig("/mnt/vesuvius/ridgecentre1218/ridge_offsets_L1.png",dpi=130,facecolor="white")
print("saved")
