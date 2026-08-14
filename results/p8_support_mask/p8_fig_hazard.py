"""Hazard A/B: same seed, no weight vs hazard weight, flattened renders."""
import json
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
P8="/mnt/vesuvius/p8_sprint"
pw=json.load(open("/mnt/vesuvius/vcbuild/demo_out/powered/powered_scores.json"))
runs=pw["runs"] if isinstance(pw,dict) and "runs" in pw else pw
def cell(site,arm):
    v=[r for r in runs if r.get("site","").startswith(site+"_") and r.get("arm")==arm]
    a=np.mean([r["area_vx2"] for r in v]); 
    d=[r["ray"].get("frac_double") for r in v if r["ray"].get("frac_double") is not None]
    return a, (100*np.mean(d) if d else float("nan")), len(v)
sites=["s0","s1","s3"]
fig,axes=plt.subplots(3,2,figsize=(9.0,13.2))
for i,s in enumerate(sites):
    for j,arm in enumerate(("A","B")):
        im=np.array(Image.open("%s/render_hazard/%s_%s/00.tif"%(P8,s,arm)))
        ax=axes[i,j]; ax.imshow(im,cmap="gray"); ax.set_xticks([]); ax.set_yticks([])
        area,dbl,n=cell(s,arm)
        ax.set_title("%s  %s\n%.0fk vx2 traced, %.1f%% double-thickness (mean of %d reps)"
                     %(s,"direction field only" if arm=="A" else "+ hazard weight",
                       area/1000.0,dbl,n),fontsize=9)
fig.suptitle("PHerc1218: the fused-site census as a GrowPatch weight volume\n"
             "same seed and parameters, one replicate rendered per arm through vc_render_tifxyz",
             fontsize=11)
fig.tight_layout(rect=[0,0,1,0.965])
fig.savefig(P8+"/render_hazard_ab.png",dpi=130)
print("saved")
for s in sites:
    for arm in ("A","B"):
        a,d,n=cell(s,arm); print(s,arm,"area %.0f dbl %.2f n %d"%(a,d,n))
