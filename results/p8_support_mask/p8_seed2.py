"""Seeds on supported surface CLOSE to the phantom halo (fixed rule)."""
import hashlib, json, sys
import numpy as np, zarr
from scipy import ndimage as ndi
ROOT="/mnt/vesuvius/p8_sprint"; tag=sys.argv[1]
rep=json.load(open("%s/%s_repair.json"%(ROOT,tag))); lo,hi=rep["bounds_l0"]
b=zarr.open("%s/%s_before.zarr/0"%(ROOT,tag),mode="r")
s=zarr.open("%s/%s_support.zarr/0"%(ROOT,tag),mode="r")
SS=2  # subsample within plane
cands=[]
for z in range(lo[0]+128, hi[0]-128, 128):
    zb=np.asarray(b[z, lo[1]:hi[1]:SS, lo[2]:hi[2]:SS])>=128
    zs=np.asarray(s[z, lo[1]:hi[1]:SS, lo[2]:hi[2]:SS])>0
    sup=zb&zs; ph=zb&~zs
    if sup.sum()<500 or ph.sum()<500: continue
    d=ndi.distance_transform_edt(~ph)          # samples to nearest phantom
    m=sup & (d<=16) & (d>=2)                   # 4..32 L0 vox from phantom
    if not m.any(): continue
    ys,xs=np.nonzero(m)
    # local supported mass so the seed sits on real sheet, not a speck
    for k in range(0,len(ys),max(1,len(ys)//25)):
        yy,xx=int(ys[k]),int(xs[k])
        loc=sup[max(0,yy-12):yy+12, max(0,xx-12):xx+12].sum()
        if loc<120: continue
        Y=lo[1]+SS*yy; X=lo[2]+SS*xx
        cands.append({"l0":[z,Y,X],"dist_to_phantom_vox":int(d[yy,xx]*SS),
                      "local_supported":int(loc),
                      "key":hashlib.md5(("%d:%d:%d"%(z,Y,X)).encode()).hexdigest()})
cands.sort(key=lambda c:(-c["local_supported"], c["key"]))
picked=[]
for c in cands:
    if all(abs(c["l0"][0]-p["l0"][0])>256 or abs(c["l0"][1]-p["l0"][1])>256
           or abs(c["l0"][2]-p["l0"][2])>256 for p in picked): picked.append(c)
    if len(picked)>=3: break
print(json.dumps({"n_cand":len(cands),"picked":picked},indent=1))
json.dump({"rule":"supported m7 positive 4-32 L0 vox from the phantom halo, local supported mass >=120 samples, rank by local mass then md5, min sep 256 L0","n_cand":len(cands),"picked":picked},
          open("%s/%s_seeds.json"%(ROOT,tag),"w"),indent=1)
