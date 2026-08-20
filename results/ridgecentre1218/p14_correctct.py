"""p14_ridge2.py with one change: the CT comes from level 1, the label's own grid.

The shipped crop intensity was read out of CT level 0 at level-1 label indices,
so it is a different physical region. Everything else here is byte-identical to
the published estimator.
"""
import glob, json, os
import numpy as np, zarr

CT=("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/"
    "PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr")
SRC="/mnt/vesuvius/kaggle_tightgap1218/crops"; N=128
AXIS=np.array([7593/2.0,7593/2.0]); STEP=0.25; NPTS=12
rng=np.random.default_rng(0)
L1=zarr.open(CT,mode="r")["1"]
print("CT level1",L1.shape,flush=True)

def sample_profile(vol,p,n,half):
    ts=np.arange(-half,half+1e-9,STEP); out=[]
    for t in ts:
        q=p+t*n; i=np.floor(q).astype(int); f=q-i
        if np.any(i<0) or np.any(i+1>=np.array(vol.shape)): return None,None
        c=vol[i[0]:i[0]+2,i[1]:i[1]+2,i[2]:i[2]+2].astype(np.float32)
        if c.shape!=(2,2,2): return None,None
        w=np.array([[[(1-f[0])*(1-f[1])*(1-f[2]),(1-f[0])*(1-f[1])*f[2]],
                     [(1-f[0])*f[1]*(1-f[2]),(1-f[0])*f[1]*f[2]]],
                    [[f[0]*(1-f[1])*(1-f[2]),f[0]*(1-f[1])*f[2]],
                     [f[0]*f[1]*(1-f[2]),f[0]*f[1]*f[2]]]])
        out.append(float((c*w).sum()))
    return ts,np.array(out)

def run_centre(lab,p,n,iid,half=20):
    ext={}
    for sgn,key in ((1,"hi"),(-1,"lo")):
        t=0.0
        while t<half:
            t+=0.5
            q=np.round(p+sgn*t*n).astype(int)
            if np.any(q<0) or np.any(q>=np.array(lab.shape)): break
            if int(lab[q[0],q[1],q[2]])!=iid: break
        ext[key]=t-0.5
    if ext["hi"]+ext["lo"]<1.0: return None
    return (ext["hi"]-ext["lo"])/2.0

rows=[]; disc={"no_ids":0,"few_points":0,"no_run":0,"corridor":0,"flat":0}; done=0
for f in sorted(glob.glob(SRC+"/*.npz")):
    d=np.load(f,allow_pickle=True)
    A,B=int(d["A_id"]),int(d["B_id"]); lab=d["instance"]
    ids=set(np.unique(lab).tolist())
    if A not in ids or B not in ids: disc["no_ids"]+=1; continue
    sz,sy,sx=[int(v) for v in d["site"]]
    z0,y0,x0=sz-N//2,sy-N//2,sx-N//2
    vol=np.asarray(L1[z0:z0+N,y0:y0+N,x0:x0+N])
    ca=np.argwhere(lab==A).mean(0); cb=np.argwhere(lab==B).mean(0)
    n=cb-ca; nn=np.linalg.norm(n)
    if nn<1e-6: disc["no_run"]+=1; continue
    n=n/nn
    site=np.array([float(v) for v in d["site"]])
    if np.dot(n[1:],site[1:]-AXIS)<0: n=-n
    offs={3:[],4:[],8:[]}; pts_signed=[]
    for iid in (A,B):
        pts=np.argwhere(lab==iid)
        if len(pts)<NPTS: disc["few_points"]+=1; continue
        for p in pts[rng.choice(len(pts),size=NPTS,replace=False)].astype(float):
            c=run_centre(lab,p,n,iid)
            if c is None: disc["no_run"]+=1; continue
            centre=p+c*n
            for half in (3,4,8):
                ts,prof=sample_profile(vol,centre,n,half)
                if ts is None: disc["corridor"]+=1; continue
                if float(prof.max()-prof.min())<1.0: disc["flat"]+=1; continue
                v=float(ts[int(np.argmax(prof))]); offs[half].append(v)
                if half==4: pts_signed.append(v)
    if not offs[4]: continue
    rows.append({"file":os.path.basename(f),"band":str(d["band"]),"gap":float(d["gap"]),
                 **{"off_%d"%h: float(np.median(offs[h])) for h in (3,4,8) if offs[h]},
                 "n_points":len(offs[4]),"frac_neg_within":float(np.mean(np.array(pts_signed)<0))})
    done+=1
    if done%25==0: print("  %d crops"%done,flush=True)

def summ(key):
    v=np.array([r[key] for r in rows if key in r])
    if not len(v): return None
    r2=np.random.default_rng(0)
    boot=[float(np.median(r2.choice(v,len(v)))) for _ in range(2000)]
    return dict(n=len(v),median=round(float(np.median(v)),4),
                median_abs=round(float(np.median(np.abs(v))),4),
                frac_negative=round(float((v<0).mean()),4),
                ci95=[round(float(np.percentile(boot,2.5)),4),round(float(np.percentile(boot,97.5)),4)])

fw=np.array([r["frac_neg_within"] for r in rows])
hist,_=np.histogram(fw,bins=[0,.15,.35,.65,.85,1.001])
out=dict(ct="level 1, the label grid", crops_used=len(rows), discarded=disc,
         corridors={k:summ("off_%d"%k) for k in (3,4,8)},
         within_crop_sign=dict(bins="0-.15|.15-.35|.35-.65|.65-.85|.85-1",counts=hist.tolist(),
                               median=round(float(np.median(fw)),4)))
json.dump({"summary":out,"rows":rows},open("/mnt/vesuvius/ridgecentre1218/ridge_offsets_L1.json","w"),indent=1)
print(json.dumps(out,indent=1))
