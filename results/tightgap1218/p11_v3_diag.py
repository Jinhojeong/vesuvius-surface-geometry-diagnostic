"""Why is corr(CT, label) low on the v3 crops when the CT is byte-verified correct?

Two hypotheses.
  H1 the labels are misaligned in v3, in which case a deliberately SHIFTED CT
     scores as well as the true one.
  H2 the v3 sites sit in dense regions where CT barely distinguishes sheet from
     gap, in which case the true CT beats every shift and the contrast is real
     but small.
The shift test separates them without depending on the absolute correlation.
"""
import glob, json
import numpy as np, zarr

CT=("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/"
    "PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr")
V3="/mnt/vesuvius/tightgap1218_v3"; V2="/mnt/vesuvius/kaggle_tightgap1218"; N=128
L1=zarr.open(CT,mode="r")["1"]
def corr(a,b):
    a=a.astype(np.float64).ravel(); b=b.astype(np.float64).ravel()
    return float(np.corrcoef(a,b)[0,1]) if a.std()>1e-9 and b.std()>1e-9 else float("nan")

rng=np.random.default_rng(1)
for tag,src in (("v3",V3+"/crops"),("v2",V2+"/crops")):
    fs=sorted(glob.glob(src+"/*.npz"))
    pick=[fs[i] for i in rng.choice(len(fs),6,replace=False)]
    print("\n===== %s ====="%tag,flush=True)
    print("%-34s %6s %7s %7s %7s %7s %8s %8s"%("crop","r_true","r_sh20","r_sh40","labfrac","empty","mu_lab","mu_bg"))
    for f in pick:
        d=np.load(f,allow_pickle=True); m=(d["instance"]>0)
        sz,sy,sx=[int(q) for q in d["site"]]; z0,y0,x0=sz-N//2,sy-N//2,sx-N//2
        v=np.asarray(L1[z0:z0+N,y0:y0+N,x0:x0+N])
        rs=[]
        for sh in (20,40):
            w=np.asarray(L1[z0+sh:z0+sh+N,y0:y0+N,x0:x0+N])
            rs.append(corr(w,m) if w.shape==(N,N,N) else float("nan"))
        vf=v.astype(np.float64)
        print("%-34s %6.3f %7.3f %7.3f %7.3f %7.3f %8.1f %8.1f"%(
            f.split("/")[-1][:34],corr(v,m),rs[0],rs[1],
            float(m.mean()),float((v==0).mean()),
            float(vf[m].mean()) if m.any() else -1, float(vf[~m].mean())),flush=True)

# aggregate over more crops: does the true CT beat shifts systematically?
print("\n=== aggregate shift test on 30 v3 crops ===",flush=True)
fs=sorted(glob.glob(V3+"/crops/*.npz"))
pick=[fs[i] for i in rng.choice(len(fs),30,replace=False)]
wins=0; n=0; dmu=[]
for f in pick:
    d=np.load(f,allow_pickle=True); m=(d["instance"]>0)
    if not m.any() or m.all(): continue
    sz,sy,sx=[int(q) for q in d["site"]]; z0,y0,x0=sz-N//2,sy-N//2,sx-N//2
    v=np.asarray(L1[z0:z0+N,y0:y0+N,x0:x0+N]); rt=corr(v,m)
    best=-9
    for sh in (-40,-20,20,40):
        if z0+sh<0 or z0+sh+N>L1.shape[0]: continue
        w=np.asarray(L1[z0+sh:z0+sh+N,y0:y0+N,x0:x0+N])
        best=max(best,corr(w,m))
    n+=1
    if rt>best: wins+=1
    vf=v.astype(np.float64); dmu.append(float(vf[m].mean()-vf[~m].mean()))
print("true CT beats every shift on %d of %d crops"%(wins,n))
print("mean CT on label minus mean CT on background: median %.1f, min %.1f, share positive %.0f%%"
      %(np.median(dmu),np.min(dmu),100*np.mean(np.array(dmu)>0)))
