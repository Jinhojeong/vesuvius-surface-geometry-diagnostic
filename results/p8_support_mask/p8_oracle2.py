"""Diagnose the 953: is it L0/L1 sampling or real over-masking?
For labeled L1 voxels whose L0 corner lost its positive, inspect the full
2x2x2 L0 CT block: if any sub-voxel has CT>0 the label keeps support."""
import glob, json, re, sys
import numpy as np, zarr

ROOT="/mnt/vesuvius/p8_sprint"; tag=sys.argv[1]
rep=json.load(open("%s/%s_repair.json"%(ROOT,tag))); lo,hi=rep["bounds_l0"]
b=zarr.open("%s/%s_before.zarr/0"%(ROOT,tag),mode="r")
a=zarr.open("%s/%s_after.zarr"%(ROOT,tag),mode="r")
if hasattr(a,"keys") and "0" in a: a=a["0"]
s=zarr.open("%s/%s_support.zarr/0"%(ROOT,tag),mode="r")
l1lo=[v//2 for v in lo]; l1hi=[v//2 for v in hi]
cases=[]; n_block_any_ct=0; n_block_no_ct=0; n_all8_removed=0
for p in glob.glob("/mnt/vesuvius/p1218_repair_v3/blocks_repaired/z*/tile_*.npz"):
    slab=p.split("/")[-2]; m=re.search(r"tile_y(\d+)_x(\d+)",p)
    z0=int(slab[1:]); y0=int(m.group(1)); x0=int(m.group(2))
    if not (l1lo[0]<z0+256 and z0<l1hi[0] and l1lo[1]<y0+512 and y0<l1hi[1]
            and l1lo[2]<x0+512 and x0<l1hi[2]): continue
    with np.load(p) as d: lab=d["labels"]
    zs,ze=max(z0,l1lo[0]),min(z0+lab.shape[0],l1hi[0])
    ys,ye=max(y0,l1lo[1]),min(y0+lab.shape[1],l1hi[1])
    xs,xe=max(x0,l1lo[2]),min(x0+lab.shape[2],l1hi[2])
    if zs>=ze or ys>=ye or xs>=xe: continue
    sub=lab[zs-z0:ze-z0,ys-y0:ye-y0,xs-x0:xe-x0]>0
    if not sub.any(): continue
    bb=b[2*zs:2*ze:2,2*ys:2*ye:2,2*xs:2*xe:2]; aa=a[2*zs:2*ze:2,2*ys:2*ye:2,2*xs:2*xe:2]
    n=[min(sub.shape[i],bb.shape[i]) for i in range(3)]
    sub=sub[:n[0],:n[1],:n[2]]; bb=bb[:n[0],:n[1],:n[2]]; aa=aa[:n[0],:n[1],:n[2]]
    idx=np.argwhere(sub & (bb>=128) & (aa<128))
    for i in idx[:400]:
        Z,Y,X=2*(zs+int(i[0])),2*(ys+int(i[1])),2*(xs+int(i[2]))
        ctblk=np.asarray(s[Z:Z+2,Y:Y+2,X:X+2]); ablk=np.asarray(a[Z:Z+2,Y:Y+2,X:X+2])
        bblk=np.asarray(b[Z:Z+2,Y:Y+2,X:X+2])
        anyct=bool((ctblk>0).any()); n_block_any_ct+=anyct; n_block_no_ct+= (not anyct)
        if (bblk>=128).any() and not (ablk>=128).any(): n_all8_removed+=1
        if len(cases)<6: cases.append({"l0":[Z,Y,X],"ct_block":ctblk.ravel().tolist(),
                                       "before":bblk.ravel().tolist(),"after":ablk.ravel().tolist()})
out={"inspected":n_block_any_ct+n_block_no_ct,"block_has_some_CT":n_block_any_ct,
     "block_entirely_CT0":n_block_no_ct,"label_fully_lost_all8":n_all8_removed,"examples":cases}
print(json.dumps(out,indent=1)[:2500])
json.dump(out,open("%s/%s_oracle2.json"%(ROOT,tag),"w"),indent=1)
