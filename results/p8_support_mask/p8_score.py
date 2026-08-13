"""Phantom-quad counts before vs after, per seed."""
import json, sys
import numpy as np, zarr
sys.path.insert(0,"/mnt/vesuvius/vcbuild")
import score_ab
ROOT="/mnt/vesuvius/p8_sprint"; tag=sys.argv[1]
s=zarr.open("%s/%s_support.zarr/0"%(ROOT,tag),mode="r")
rows=[]
for i in range(3):
    r={"seed":i}
    for arm in ("before","after"):
        d="%s/trace/%s_s%d_%s"%(ROOT,tag,i,arm)
        run=score_ab.load_tifxyz(d)
        if run is None: r[arm]=None; continue
        v=run["valid"]; pts=np.stack([run["z"][v],run["y"][v],run["x"][v]],1)
        n=len(pts); sub=pts[::max(1,n//4000)][:4000]
        sup=np.array([int(s[int(round(p[0])),int(round(p[1])),int(round(p[2]))])
                      for p in sub])
        meta=json.load(open(d+"/meta.json"))
        r[arm]={"quads":int(v.sum()),"area_vx2":meta.get("area_vx2"),
                "sampled":len(sub),"phantom_quads_frac":float((sup==0).mean()),
                "max_gen":meta.get("max_gen")}
    rows.append(r)
print(json.dumps(rows,indent=1))
json.dump(rows,open("%s/%s_quadscore.json"%(ROOT,tag),"w"),indent=1)
