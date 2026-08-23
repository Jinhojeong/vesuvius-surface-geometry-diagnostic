"""Why does the walk fail to reach B for 84 percent of discards?

Hypothesis: I sampled points anywhere on instance A. The crop is centred on the
split SITE, so A and B face each other near the centre and diverge away from it.
A global centroid-to-centroid normal then misses B from distant points. That
would be my sampling geometry, not a property of the data.
Test: reach rate as a function of the sampled point's distance from the crop centre.
"""
import glob, json
import numpy as np
V5="/mnt/vesuvius/tightgap1218_v5"; N=128; NPTS=40; STEP=0.5
rng=np.random.default_rng(0)
C=np.array([64.,64.,64.])

def reaches(lab,p,n,A,B,maxr=40):
    t=0.0; exitA=None
    while t<maxr:
        t+=STEP
        q=np.round(p+t*n).astype(int)
        if np.any(q<0) or np.any(q>=np.array(lab.shape)): return "left"
        v=int(lab[q[0],q[1],q[2]])
        if exitA is None and v!=A: exitA=t
        if exitA is not None and v==B: return "reached"
    return "never"

bins={}; overall={"reached":0,"never":0,"left":0}
for f in sorted(glob.glob(V5+"/crops/*.npz"))[:120]:
    d=np.load(f,allow_pickle=True); lab=d["instance"]
    A,B=int(d["A_id"]),int(d["B_id"])
    ca=np.argwhere(lab==A).mean(0); cb=np.argwhere(lab==B).mean(0)
    n=cb-ca; nn=np.linalg.norm(n)
    if nn<1e-6: continue
    n=n/nn
    pts=np.argwhere(lab==A)
    if len(pts)<NPTS: continue
    for p in pts[rng.choice(len(pts),size=NPTS,replace=False)].astype(float):
        r=reaches(lab,p,n,A,B)
        overall[r]+=1
        dist=float(np.linalg.norm(p-C))
        b=int(dist//10)*10
        bins.setdefault(b,{"reached":0,"never":0,"left":0})[r]+=1

print("overall over %d points: %s"%(sum(overall.values()),overall))
print("\n%-14s %7s %8s"%("dist from centre","n","reach %"))
for b in sorted(bins):
    v=bins[b]; tot=sum(v.values())
    if tot<20: continue
    print("%-14s %7d %7.1f%%"%("%d-%d vox"%(b,b+10),tot,100*v["reached"]/tot))
