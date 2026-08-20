"""Measure the validation windows of the two additional segments."""
import json
import numpy as np, zarr
OUT={}
for seg in ("pherc0814-46527","pherc1667-w029"):
    d="/mnt/vesuvius/ink9um_labels/"+seg
    g=lambda k: (lambda z: z["0"] if "0" in z else z)(zarr.open("%s/%s_%s.zarr"%(d,seg,k),mode="r"))
    V=np.asarray(g("validation_mask")); I=g("inklabels")
    nz=np.argwhere(V>0)
    if not len(nz):
        print(seg,"no validation voxels"); continue
    lo,hi=nz.min(0),nz.max(0)+1
    zc=int(np.bincount(nz[:,0]).argmax())
    m=V[zc]>0
    ink=np.asarray(I[zc])[m]
    r=dict(shape=list(V.shape), voxels=int((V>0).sum()), annotated_z=zc,
           bbox_lo=lo.tolist(), bbox_hi=hi.tolist(),
           extent=[int(hi[1]-lo[1]),int(hi[2]-lo[2])],
           masked_at_z=int(m.sum()), ink_frac=round(float((ink>0).mean()),5))
    OUT[seg]=r
    print("%-18s val voxels %7d | z %2d | window %dx%d | ink %.3f"%(
        seg,r["voxels"],zc,r["extent"][0],r["extent"][1],r["ink_frac"]))
json.dump(OUT,open("/mnt/vesuvius/ink9um_labels/probe2.json","w"),indent=1)
