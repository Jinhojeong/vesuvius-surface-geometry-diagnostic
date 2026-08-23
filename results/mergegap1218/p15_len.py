"""Is the drop at the loosest band a length bias?

"Merged" requires m7 positive at EVERY step of the gap interval, and the interval
grows with the gap. A longer interval has more chances to contain a zero, so a
lower merge rate at wide gaps could be geometry rather than prediction behaviour.
Test: per-step positive rate, which is length-free, alongside the all-steps rate.
"""
import glob, json, os, sys
sys.path.insert(0,"/mnt/vesuvius/overlap_step2")
import numpy as np
from zarr_http import RemoteZarrLevel
M7=("https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc1218/"
    "representations/predictions/surfaces/"
    "20250521120456-surface-20260413222639-surface-m7-L0-th0.2.zarr")
CACHE="/mnt/vesuvius/hazard_zarr_smoke/m7L1_cache"
V5="/mnt/vesuvius/tightgap1218_v5"; N=128; NPTS=12; STEP=0.5; RADIUS=20.0; UM=17.28
m7=RemoteZarrLevel(M7,1,cache_dir=CACHE); rng=np.random.default_rng(0)
per_band={}
for f in sorted(glob.glob(V5+"/crops/*.npz")):
    d=np.load(f,allow_pickle=True); lab=d["instance"]
    A,B=int(d["A_id"]),int(d["B_id"]); band=str(d["band"])
    sz,sy,sx=[int(v) for v in d["site"]]; z0,y0,x0=sz-N//2,sy-N//2,sx-N//2
    try: pred=m7.read_crop((z0,y0,x0),(z0+N,y0+N,x0+N))
    except Exception: continue
    ca=np.argwhere(lab==A).mean(0); cb=np.argwhere(lab==B).mean(0)
    nv=cb-ca; nn=np.linalg.norm(nv)
    if nn<1e-6: continue
    nv=nv/nn
    pts=np.argwhere(lab==A); c=np.array([64.,64.,64.])
    pts=pts[np.linalg.norm(pts-c,axis=1)<=RADIUS]
    if len(pts)<NPTS: continue
    for p in pts[rng.choice(len(pts),size=NPTS,replace=False)].astype(float):
        t=0.0; eA=None; eB=None
        while t<40:
            t+=STEP
            q=np.round(p+t*nv).astype(int)
            if np.any(q<0) or np.any(q>=np.array(lab.shape)): break
            v=int(lab[q[0],q[1],q[2]])
            if eA is None and v!=A: eA=t
            if eA is not None and v==B: eB=t; break
        if eA is None or eB is None or eB-eA<STEP: continue
        ts=np.arange(eA,eB+1e-9,STEP); vals=[]; bad=False
        for tt in ts:
            q=np.round(p+tt*nv).astype(int)
            if np.any(q<0) or np.any(q>=np.array(pred.shape)): bad=True; break
            vals.append(float(pred[q[0],q[1],q[2]]))
        if bad or not vals: continue
        e=per_band.setdefault(band,{"steps":0,"pos":0,"all":0,"n":0,"len":[]})
        e["steps"]+=len(vals); e["pos"]+=sum(1 for v in vals if v>0)
        e["all"]+= 1 if all(v>0 for v in vals) else 0
        e["n"]+=1; e["len"].append(len(vals))
BANDS=[("0-2","under 34.6"),("2-4","34.6-69.1"),("4-6","69.1-103.7"),("6-10","103.7-172.8"),("10+","above 172.8")]
print("%-8s %-14s %7s %9s %14s %14s"%("band","microns","points","med len","all-steps rate","per-step rate"))
for b,um in BANDS:
    e=per_band.get(b)
    if not e: continue
    print("%-8s %-14s %7d %9.1f %14.3f %14.3f"%(b,um,e["n"],np.median(e["len"]),e["all"]/e["n"],e["pos"]/e["steps"]))
print("\nif the all-steps rate falls with gap but the per-step rate does not,")
print("the fall is a length bias rather than prediction behaviour.")
