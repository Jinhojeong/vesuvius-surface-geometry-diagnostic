"""Pick deterministic seeds: supported real surface adjacent to phantom mass."""
import hashlib, json, sys
import numpy as np, zarr
ROOT="/mnt/vesuvius/p8_sprint"; tag=sys.argv[1]
rep=json.load(open("%s/%s_repair.json"%(ROOT,tag))); lo,hi=rep["bounds_l0"]
b=zarr.open("%s/%s_before.zarr/0"%(ROOT,tag),mode="r")
s=zarr.open("%s/%s_support.zarr/0"%(ROOT,tag),mode="r")
STEP=32; best=[]
for z in range(lo[0]+128, hi[0]-128, STEP*4):
    zb=np.asarray(b[z:z+1, lo[1]:hi[1]:2, lo[2]:hi[2]:2])[0]>=128
    zs=np.asarray(s[z:z+1, lo[1]:hi[1]:2, lo[2]:hi[2]:2])[0]>0
    sup=zb&zs; ph=zb&~zs
    if sup.sum()<200 or ph.sum()<2000: continue
    # for each supported point, count phantom within a 32-sample (=64 vox) box
    ys,xs=np.nonzero(sup)
    for k in range(0, len(ys), max(1,len(ys)//40)):
        yy,xx=int(ys[k]),int(xs[k])
        y0,y1=max(0,yy-32),yy+32; x0,x1=max(0,xx-32),xx+32
        npha=int(ph[y0:y1,x0:x1].sum())
        if npha<300: continue
        Y=lo[1]+2*yy; X=lo[2]+2*xx
        best.append({"l0":[z,Y,X],"phantom_near":npha,
                     "key":hashlib.md5(("%d:%d:%d"%(z,Y,X)).encode()).hexdigest()})
best.sort(key=lambda d:(-d["phantom_near"], d["key"]))
picked=[]; 
for c in best:
    if all(abs(c["l0"][0]-p["l0"][0])>256 or abs(c["l0"][1]-p["l0"][1])>256
           or abs(c["l0"][2]-p["l0"][2])>256 for p in picked):
        picked.append(c)
    if len(picked)>=3: break
print(json.dumps({"n_cand":len(best),"picked":picked},indent=1))
json.dump({"rule":"supported m7 positive with >=300 phantom samples within 64 L0 vox; rank by phantom_near then md5; min separation 256 L0","n_cand":len(best),"picked":picked},
          open("%s/%s_seeds.json"%(ROOT,tag),"w"),indent=1)
