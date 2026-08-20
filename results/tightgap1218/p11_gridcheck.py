"""Independent check: are the shipped crop intensities the CT at the label's own location?

If the crops are sound, intensity should be high exactly where the surface label is.
If the CT was read at level-1 label indices out of a level-0 array, the two arrays
describe different physical boxes and the association collapses.
"""
import glob, json
import numpy as np, zarr

CT=("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/"
    "PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr")
BLOCKS="/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired"
BASE="/mnt/vesuvius/kaggle_tightgap1218"; N=128

z=zarr.open(CT,mode="r")
print("CT pyramid:")
for k in sorted(z.array_keys()): print("   level",k,z[k].shape)

# what grid are the repair labels on?
zs=[]
for f in sorted(glob.glob(BLOCKS+"/*/*.npz"))[:4000]:
    n=f.split("/")[-2]
    if n.startswith("z"): zs.append(int(n[1:]))
print("\nrepair block z0: min %d max %d  -> label grid z extent ~ %d"%(min(zs),max(zs),max(zs)+256))
print("CT level0 z=%d  level1 z=%d"%(z["0"].shape[0],z["1"].shape[0]))

def corr(a,b):
    a=a.astype(np.float64).ravel(); b=b.astype(np.float64).ravel()
    if a.std()<1e-9 or b.std()<1e-9: return float("nan")
    return float(np.corrcoef(a,b)[0,1])

print("\n%-42s %8s %8s %8s %9s %9s"%("crop","ship","L1@idx","L0@2idx","ship0frac","L1 0frac"))
for f in sorted(glob.glob(BASE+"/crops/*.npz"))[:6]:
    d=np.load(f,allow_pickle=True)
    lab=(d["instance"]>0); ship=d["intensity"]
    sz,sy,sx=[int(v) for v in d["site"]]
    z0,y0,x0=sz-N//2,sy-N//2,sx-N//2
    l1=np.asarray(z["1"][z0:z0+N,y0:y0+N,x0:x0+N])
    # same physical box read from level 0: double the indices, take 256^3, mean-pool 2x
    a=np.asarray(z["0"][2*z0:2*z0+2*N,2*y0:2*y0+2*N,2*x0:2*x0+2*N]).astype(np.float32)
    l0=a.reshape(N,2,N,2,N,2).mean((1,3,5)) if a.shape==(2*N,2*N,2*N) else None
    print("%-42s %8.3f %8.3f %8s %9.4f %9.4f"%(
        f.split("/")[-1][:42], corr(ship,lab), corr(l1,lab),
        ("%.3f"%corr(l0,lab)) if l0 is not None else "n/a",
        float((ship==0).mean()), float((l1==0).mean())),flush=True)
