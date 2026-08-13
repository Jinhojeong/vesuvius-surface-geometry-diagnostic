"""Assemble the before/after flattened-render comparison figure."""
import json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
P8="/mnt/vesuvius/p8_sprint"
q=json.load(open(P8+"/r1_quadscore_v2.json"))
fig,axes=plt.subplots(3,2,figsize=(9.0,13.2))
for i in range(3):
    for j,arm in enumerate(("before","after")):
        im=np.array(Image.open("%s/render/s%d_%s/00.tif"%(P8,i,arm)))
        ax=axes[i,j]; ax.imshow(im,cmap="gray"); ax.set_xticks([]); ax.set_yticks([])
        d=q[i][arm]
        ax.set_title("seed %d  %s\n%d quads in region, %.1f%% over unsupported CT"
                     %(i,"original m7" if arm=="before" else "support-masked m7",
                       d["quads_in_region"],100*d["unsupported_frac"]),fontsize=9)
        nz=float((im>0).mean())
        ax.set_xlabel("rendered texture coverage %.1f%%"%(100*nz),fontsize=8)
fig.suptitle("PHerc1218: surfaces traced from the m7 prediction before and after\n"
             "CT-support masking (villa PR #1156), rendered through vc_render_tifxyz",
             fontsize=11)
fig.tight_layout(rect=[0,0,1,0.965])
fig.savefig(P8+"/render_before_after.png",dpi=130)
print("saved", P8+"/render_before_after.png")
for i in range(3):
    for arm in ("before","after"):
        im=np.array(Image.open("%s/render/s%d_%s/00.tif"%(P8,i,arm)))
        print("s%d %-6s shape %s dtype %s nonzero %.3f mean %.1f"%(i,arm,im.shape,im.dtype,float((im>0).mean()),float(im.mean())))
