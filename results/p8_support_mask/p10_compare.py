"""PR #1451 at-source rule vs our post-hoc trim, same seeds."""
import json, sys
import numpy as np, zarr
sys.path.insert(0,"/mnt/vesuvius/vcbuild"); import score_ab
P8="/mnt/vesuvius/p8_sprint"
lo,hi=json.load(open(P8+"/r1_repair.json"))["bounds_l0"]
sup=zarr.open(P8+"/r1_support.zarr/0",mode="r")
trim={r["seed"]:r for r in json.load(open(P8+"/r1_trim.json"))["rows"] if r["arm"]=="after"}
rows=[]
for i in range(3):
    base=score_ab.load_tifxyz("%s/trace/r1_s%d_after"%(P8,i))
    rule=score_ab.load_tifxyz("%s/trace1451/s%d"%(P8,i))
    if rule is None: print("no rule run for s%d"%i); continue
    def stats(run):
        v=run["valid"]; pts=np.stack([run["z"],run["y"],run["x"]],-1)
        inb=((pts>=np.array(lo))&(pts<np.array(hi))).all(-1)&v
        idx=np.argwhere(inb)
        s=np.array([int(sup[int(round(pts[a,b,0])),int(round(pts[a,b,1])),int(round(pts[a,b,2]))]) for a,b in idx])
        return int(v.sum()), int(inb.sum()), float((s>0).mean()) if len(s) else float("nan")
    qb,ib,fb=stats(base); qr,ir,fr=stats(rule)
    mb=json.load(open("%s/trace/r1_s%d_after/meta.json"%(P8,i)))
    mr=json.load(open("%s/trace1451/s%d/meta.json"%(P8,i)))
    rej=0
    for line in open("%s/logs/t1451_s%d.log"%(P8,i)):
        if "support test rejected" in line:
            rej+=int(line.split("rejected")[1].split("of")[0])
    r=dict(seed=i,
        base_quads=qb, rule_quads=qr,
        base_in=ib, rule_in=ir,
        base_supported_frac=round(fb,4), rule_supported_frac=round(fr,4),
        base_area=round(mb["area_vx2"]), rule_area=round(mr["area_vx2"]),
        posthoc_retained=trim[i]["retained_frac"],
        atsource_retained_vs_base=round(ir/ib,4) if ib else None,
        rejects_logged=rej)
    rows.append(r); print(json.dumps(r))
json.dump(rows,open(P8+"/pr1451_compare.json","w"),indent=1)
