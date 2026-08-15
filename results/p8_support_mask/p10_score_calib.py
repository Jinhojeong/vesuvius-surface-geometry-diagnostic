"""Score the dilation calibration: at-source rule (v2) vs baseline and post-hoc."""
import json, sys, glob
import numpy as np, zarr
sys.path.insert(0,"/mnt/vesuvius/vcbuild"); import score_ab
P8="/mnt/vesuvius/p8_sprint"
lo,hi=json.load(open(P8+"/r1_repair.json"))["bounds_l0"]
sup=zarr.open(P8+"/r1_support.zarr/0",mode="r")
trim={r["seed"]:r for r in json.load(open(P8+"/r1_trim.json"))["rows"] if r["arm"]=="after"}
def stats(path):
    run=score_ab.load_tifxyz(path)
    if run is None: return None
    v=run["valid"]; pts=np.stack([run["z"],run["y"],run["x"]],-1)
    inb=((pts>=np.array(lo))&(pts<np.array(hi))).all(-1)&v
    idx=np.argwhere(inb)
    if not len(idx): return dict(quads=int(v.sum()),inreg=0,supfrac=None)
    s=np.array([int(sup[int(round(pts[a,b,0])),int(round(pts[a,b,1])),int(round(pts[a,b,2]))]) for a,b in idx])
    return dict(quads=int(v.sum()),inreg=int(inb.sum()),supfrac=round(float((s>0).mean()),4))
out={"baseline":{},"rule":{}}
for i in range(3):
    b=stats("%s/trace/r1_s%d_after"%(P8,i)); out["baseline"][i]=b
for f in sorted(glob.glob(P8+"/calib1451/d*_s*")):
    tag=f.split("/")[-1]; d,s=tag.split("_"); s=int(s[1:]); d=int(d[1:])
    st=stats(f)
    if st is None: continue
    base=out["baseline"][s]
    st["retained_vs_base"]=round(st["inreg"]/base["inreg"],4) if base["inreg"] else None
    st["posthoc_target"]=trim[s]["retained_frac"]
    out["rule"][tag]=st
    print(tag, json.dumps(st))
json.dump(out,open(P8+"/pr1451_calib.json","w"),indent=1)
print("saved")
