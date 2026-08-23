"""Does the published m7 surface prediction merge two touching sheets?

Per PREREGISTRATION.md (a8ae5f7740a5079a). Measured on the 300 v5 contact crops
plus the 60 controls as a negative control. m7 read at level 1, the same grid the
labels and crops sit on.
"""
import glob, json, os, sys
sys.path.insert(0,"/mnt/vesuvius/overlap_step2")
import numpy as np
from zarr_http import RemoteZarrLevel

M7=("https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc1218/"
    "representations/predictions/surfaces/"
    "20250521120456-surface-20260413222639-surface-m7-L0-th0.2.zarr")
CACHE="/mnt/vesuvius/hazard_zarr_smoke/m7L1_cache"
V5="/mnt/vesuvius/tightgap1218_v5"; OUT="/mnt/vesuvius/mergegap1218"
N=128; NPTS=12; STEP=0.5; UM=17.28
os.makedirs(OUT,exist_ok=True)
m7=RemoteZarrLevel(M7,1,cache_dir=CACHE)
print("m7 level1",m7.shape,flush=True)

def walk(lab,p,n,A,B,maxr=40):
    """from a point on A, walk toward B. return (exit_of_A, entry_of_B) in steps."""
    exitA=None; entryB=None
    t=0.0
    while t<maxr:
        t+=STEP
        q=np.round(p+t*n).astype(int)
        if np.any(q<0) or np.any(q>=np.array(lab.shape)): return None,None
        v=int(lab[q[0],q[1],q[2]])
        if exitA is None and v!=A: exitA=t
        if exitA is not None and v==B: entryB=t; break
    return exitA,entryB

rows=[]; disc={"empty_interval":0,"left_crop":0,"never_reached_B":0,"few_points":0}
ctrl_scored=0
rng=np.random.default_rng(0)
for arm in ("crops","control"):
    for f in sorted(glob.glob("%s/%s/*.npz"%(V5,arm))):
        d=np.load(f,allow_pickle=True)
        lab=d["instance"]; sz,sy,sx=[int(v) for v in d["site"]]
        z0,y0,x0=sz-N//2,sy-N//2,sx-N//2
        try: pred=m7.read_crop((z0,y0,x0),(z0+N,y0+N,x0+N))
        except Exception: continue
        if arm=="control":
            A=int(np.unique(lab[lab>0])[0]) if (lab>0).any() else None
            B=-1
        else:
            A,B=int(d["A_id"]),int(d["B_id"])
        if A is None: continue
        pts=np.argwhere(lab==A)
        if len(pts)<NPTS: disc["few_points"]+=1; continue
        if arm=="crops":
            ca=np.argwhere(lab==A).mean(0); cb=np.argwhere(lab==B).mean(0)
            nvec=cb-ca; nn=np.linalg.norm(nvec)
            if nn<1e-6: continue
            nvec=nvec/nn
        else:
            nvec=np.array([0.,0.,1.])
        merged=0; tot=0
        for p in pts[rng.choice(len(pts),size=NPTS,replace=False)].astype(float):
            eA,eB=walk(lab,p,nvec,A,B)
            if eA is None: disc["left_crop"]+=1; continue
            if eB is None: disc["never_reached_B"]+=1; continue
            if eB-eA<STEP: disc["empty_interval"]+=1; continue
            ts=np.arange(eA,eB+1e-9,STEP)
            vals=[]
            bad=False
            for t in ts:
                q=np.round(p+t*nvec).astype(int)
                if np.any(q<0) or np.any(q>=np.array(pred.shape)): bad=True; break
                vals.append(float(pred[q[0],q[1],q[2]]))
            if bad: disc["left_crop"]+=1; continue
            tot+=1
            if all(v>0 for v in vals): merged+=1
        if arm=="control":
            if tot>0: ctrl_scored+=1
            continue
        if tot==0: continue
        rows.append({"file":os.path.basename(f),"band":str(d["band"]),
                     "gap_vox":float(d["gap"]),"gap_um":round(float(d["gap"])*UM,1),
                     "n_points":tot,"merged":merged,"merge_frac":merged/tot})
    print("done",arm,flush=True)

print("\nCONTROL crops that produced any merge score:",ctrl_scored,"(prereg: must be 0)",flush=True)
print("discards:",disc,flush=True)
kept=sum(r["n_points"] for r in rows); dtot=sum(disc.values())
print("sampled points kept %d, discarded %d, rate %.1f%%"%(kept,dtot,100*dtot/(kept+dtot)),flush=True)

BANDS=[("0-2","under 34.6 um"),("2-4","34.6-69.1"),("4-6","69.1-103.7"),
       ("6-10","103.7-172.8"),("10+","above 172.8")]
print("\n%-8s %-16s %6s %10s %12s"%("band","microns","crops","median","pooled"),flush=True)
summ={}
for b,lab_um in BANDS:
    v=[r["merge_frac"] for r in rows if r["band"]==b]
    pm=sum(r["merged"] for r in rows if r["band"]==b)
    pt=sum(r["n_points"] for r in rows if r["band"]==b)
    if not v: continue
    summ[b]=dict(n=len(v),median=round(float(np.median(v)),4),
                 pooled=round(pm/pt,4) if pt else None,points=pt,microns=lab_um)
    print("%-8s %-16s %6d %10.3f %12.3f"%(b,lab_um,len(v),np.median(v),pm/pt),flush=True)

json.dump({"prereg":"a8ae5f7740a5079a","m7":M7,"level":1,
           "control_scored":ctrl_scored,"discards":disc,
           "points_kept":kept,"bands":summ,"rows":rows},
          open(OUT+"/merge_by_gap.json","w"),indent=1)
print("\nwritten",flush=True)
