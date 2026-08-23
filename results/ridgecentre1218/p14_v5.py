"""Ridge versus label-run centre on the version-5 membership.

Per PREREGISTRATION.md (ad186ffbf127216b) and its AMENDMENT_1 (953017184b93ba59).
Estimator unchanged from p14_ridge2.py. Data is v5, all 300 contact crops.
CT is the crop's own shipped intensity, verified byte-identical to CT level 1.
Primary axis is fitted from the crops' own normals; the level-1 volume centre and
the axis-free within-crop test are reported alongside.
"""
import glob, json, os
import numpy as np

SRC="/mnt/vesuvius/tightgap1218_v5/crops"; OUT="/mnt/vesuvius/ridgecentre1218"
STEP=0.25; NPTS=12
os.makedirs(OUT,exist_ok=True)
rng=np.random.default_rng(0)

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

# pass 1: collect normals and sites so the axis can be fitted from the data
recs=[]
for f in sorted(glob.glob(SRC+"/*.npz")):
    d=np.load(f,allow_pickle=True)
    A,B=int(d["A_id"]),int(d["B_id"]); lab=d["instance"]
    ids=set(np.unique(lab).tolist())
    if A not in ids or B not in ids: continue
    ca=np.argwhere(lab==A).mean(0); cb=np.argwhere(lab==B).mean(0)
    n=cb-ca; nn=np.linalg.norm(n)
    if nn<1e-6: continue
    recs.append((f,n/nn,np.array([float(v) for v in d["site"]])))
print("crops with both ids:",len(recs),flush=True)
U=np.array([r[1][1:]/max(np.linalg.norm(r[1][1:]),1e-9) for r in recs])
S=np.array([r[2][1:] for r in recs])
P=np.stack([-U[:,1],U[:,0]],1)
AX_FIT,*_=np.linalg.lstsq(P,(P*S).sum(1),rcond=None)
AX_VOL=np.array([3797/2.0,3797/2.0])
print("fitted axis y=%.0f x=%.0f | level-1 volume centre %.1f"%(AX_FIT[0],AX_FIT[1],AX_VOL[0]),flush=True)

def run(axis,label):
    rng2=np.random.default_rng(0)
    rows=[]; disc={"no_run":0,"corridor":0,"flat":0,"few_points":0}
    for f,n0,site in recs:
        d=np.load(f,allow_pickle=True); A,B=int(d["A_id"]),int(d["B_id"])
        lab=d["instance"]; vol=d["intensity"]
        n=n0.copy()
        if axis is not None and np.dot(n[1:],site[1:]-axis)<0: n=-n
        offs={3:[],4:[],8:[]}; pts=[]
        for iid in (A,B):
            q=np.argwhere(lab==iid)
            if len(q)<NPTS: disc["few_points"]+=1; continue
            for p in q[rng2.choice(len(q),size=NPTS,replace=False)].astype(float):
                c=run_centre(lab,p,n,iid)
                if c is None: disc["no_run"]+=1; continue
                centre=p+c*n
                for half in (3,4,8):
                    ts,prof=sample_profile(vol,centre,n,half)
                    if ts is None: disc["corridor"]+=1; continue
                    if float(prof.max()-prof.min())<1.0: disc["flat"]+=1; continue
                    v=float(ts[int(np.argmax(prof))]); offs[half].append(v)
                    if half==4: pts.append(v)
        if not offs[4]: continue
        rows.append({"file":os.path.basename(f),"band":str(d["band"]),"gap":float(d["gap"]),
                     **{"off_%d"%h: float(np.median(offs[h])) for h in (3,4,8) if offs[h]},
                     "n_points":len(offs[4]),"frac_neg_within":float(np.mean(np.array(pts)<0))})
    def summ(key):
        v=np.array([r[key] for r in rows if key in r])
        if not len(v): return None
        r2=np.random.default_rng(0)
        boot=[float(np.median(r2.choice(v,len(v)))) for _ in range(4000)]
        return dict(n=len(v),median=round(float(np.median(v)),4),
                    median_abs=round(float(np.median(np.abs(v))),4),
                    frac_negative=round(float((v<0).mean()),4),
                    ci95=[round(float(np.percentile(boot,2.5)),4),round(float(np.percentile(boot,97.5)),4)])
    fw=np.array([r["frac_neg_within"] for r in rows])
    hist,_=np.histogram(fw,bins=[0,.15,.35,.65,.85,1.001])
    out=dict(axis=label,crops_used=len(rows),discarded=disc,
             corridors={k:summ("off_%d"%k) for k in (3,4,8)},
             within_crop_sign=dict(bins="0-.15|.15-.35|.35-.65|.65-.85|.85-1",
                                   counts=hist.tolist(),median=round(float(np.median(fw)),4)))
    print("\n=== %s ==="%label,flush=True)
    print(json.dumps(out,indent=1),flush=True)
    return out,rows

a,rows=run(AX_FIT,"axis fitted from the crops' own normals")
b,_=run(AX_VOL,"level-1 volume centre")
json.dump({"prereg":"ad186ffbf127216b","amendment":"953017184b93ba59",
           "membership":"tight-contact v5, 300 contact crops",
           "fitted_axis":[round(float(AX_FIT[0]),1),round(float(AX_FIT[1]),1)],
           "primary":a,"sensitivity_volume_centre":b,"rows":rows},
          open(OUT+"/ridge_offsets_v5.json","w"),indent=1)
print("\nwritten",flush=True)
