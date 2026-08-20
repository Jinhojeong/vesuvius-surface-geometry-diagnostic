"""Is the null a property of the data, or of the axis constant I picked?

(a) axis-free: within a crop, do the ~24 point offsets share a sign?
(b) fit the winding axis from the crops' own normals, redo the sign split under it.
Functions below are copied verbatim from p14_ridge2.py.
"""
import glob, os
import numpy as np

SRC="/mnt/vesuvius/kaggle_tightgap1218/crops"; STEP=0.25; NPTS=12
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

recs=[]
for f in sorted(glob.glob(SRC+"/*.npz")):
    d=np.load(f,allow_pickle=True)
    A,B=int(d["A_id"]),int(d["B_id"]); lab,vol=d["instance"],d["intensity"]
    ids=set(np.unique(lab).tolist())
    if A not in ids or B not in ids: continue
    ca=np.argwhere(lab==A).mean(0); cb=np.argwhere(lab==B).mean(0)
    n=cb-ca; nn=np.linalg.norm(n)
    if nn<1e-6: continue
    n=n/nn
    site=np.array([float(v) for v in d["site"]])
    pt=[]
    for iid in (A,B):
        pts=np.argwhere(lab==iid)
        if len(pts)<NPTS: continue
        sel=pts[rng.choice(len(pts),size=NPTS,replace=False)]
        for p in sel.astype(float):
            c=run_centre(lab,p,n,iid)
            if c is None: continue
            ts,prof=sample_profile(vol,p+c*n,n,4)
            if ts is None or float(prof.max()-prof.min())<1.0: continue
            pt.append(float(ts[int(np.argmax(prof))]))
    if pt: recs.append(dict(file=os.path.basename(f),n=n,site=site,pts=np.array(pt)))
print("crops:",len(recs),"points:",sum(len(r["pts"]) for r in recs))

fr=np.array([(r["pts"]<0).mean() for r in recs])
hist,_=np.histogram(fr,bins=[0,.15,.35,.65,.85,1.001])
print("\n(a) per-crop fraction of its points that are negative (unoriented normal)")
print("    0-.15 | .15-.35 | .35-.65 | .65-.85 | .85-1   ->",hist.tolist())
print("    median %.3f, share in the middle bin %.1f%%"%(np.median(fr),100*hist[2]/len(fr)))

U=np.array([r["n"][1:]/max(np.linalg.norm(r["n"][1:]),1e-9) for r in recs])
S=np.array([r["site"][1:] for r in recs])
P=np.stack([-U[:,1],U[:,0]],1)
axis,*_=np.linalg.lstsq(P,(P*S).sum(1),rcond=None)
print("\n(b) axis fitted from the crops' own normals: y=%.0f x=%.0f  (script assumed 3796, 3796)"%(axis[0],axis[1]))
rr=np.linalg.norm(S-axis,axis=1)
print("    site radius under fitted axis: min %.0f med %.0f max %.0f"%(rr.min(),np.median(rr),rr.max()))
print("\n    orientation rule                 median   share negative")
for name,ax in (("script constant (3796,3796)",np.array([3796.5,3796.5])),
                ("fitted axis",axis),("site centroid",S.mean(0))):
    med=np.array([float(np.median(r["pts"]*(1.0 if np.dot(r["n"][1:],r["site"][1:]-ax)>=0 else -1.0))) for r in recs])
    print("    %-30s %+.3f    %.1f%%"%(name,np.median(med),100*(med<0).mean()))

# (c) does dropping the runs that hit run_centre's 20-voxel walk cap move anything?
kept=[]; drop=0
for f in sorted(glob.glob(SRC+"/*.npz")):
    d=np.load(f,allow_pickle=True)
    A,B=int(d["A_id"]),int(d["B_id"]); lab,vol=d["instance"],d["intensity"]
    ids=set(np.unique(lab).tolist())
    if A not in ids or B not in ids: continue
    ca=np.argwhere(lab==A).mean(0); cb=np.argwhere(lab==B).mean(0)
    n=cb-ca; nn=np.linalg.norm(n)
    if nn<1e-6: continue
    n=n/nn
    site=np.array([float(v) for v in d["site"]])
    if np.dot(n[1:],site[1:]-np.array([3796.5,3796.5]))<0: n=-n
    pt=[]
    for iid in (A,B):
        pts=np.argwhere(lab==iid)
        if len(pts)<NPTS: continue
        for p in pts[rng.choice(len(pts),size=NPTS,replace=False)].astype(float):
            ext={}; capped=False
            for sgn,key in ((1,"hi"),(-1,"lo")):
                t=0.0
                while t<20:
                    t+=0.5
                    q=np.round(p+sgn*t*n).astype(int)
                    if np.any(q<0) or np.any(q>=np.array(lab.shape)): break
                    if int(lab[q[0],q[1],q[2]])!=iid: break
                else: capped=True
                ext[key]=t-0.5
            if ext["hi"]+ext["lo"]<1.0: continue
            if capped: drop+=1; continue
            ts,prof=sample_profile(vol,p+((ext["hi"]-ext["lo"])/2.0)*n,n,4)
            if ts is None or float(prof.max()-prof.min())<1.0: continue
            pt.append(float(ts[int(np.argmax(prof))]))
    if pt: kept.append(float(np.median(pt)))
k=np.array(kept)
print("\n(c) excluding the %d runs that hit the 20-voxel walk cap: crops %d, median %+.3f, neg %.1f%%, median|off| %.3f"
      %(drop,len(k),np.median(k),100*(k<0).mean(),np.median(np.abs(k))))
