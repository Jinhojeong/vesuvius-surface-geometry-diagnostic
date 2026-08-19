"""P11 control arm, per AMENDMENT_1 (a07c2f86f39b3bd8).

Control site: a labelled voxel whose local run thickness along the label
normal is at or below the tile median, drawn from the repaired label blocks
in block order, seed 0 for the within-block pick. All other frozen rules
unchanged.
"""
import glob, json, os
import numpy as np, zarr

OUT="/mnt/vesuvius/tightgap1218"
CT=("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/"
    "PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr")
BLOCKS="/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired"
N=128; TARGET=60
rng=np.random.default_rng(0)
ct=zarr.open(CT,mode="r"); ct0=ct["0"]; Z,Y,X=ct0.shape

def origin(slab,tile):
    p=tile.split("_"); return int(slab[1:]), int(p[1][1:]), int(p[2][1:])

def run_thickness(lab, p, inst):
    """thickness of the labelled run through p along the coarse normal."""
    z,y,x=p; best=None
    for ax in range(3):
        n=0
        for sgn in (1,-1):
            q=[z,y,x]
            for t in range(1,40):
                q[ax]=[z,y,x][ax]+sgn*t
                if not (0<=q[0]<lab.shape[0] and 0<=q[1]<lab.shape[1] and 0<=q[2]<lab.shape[2]): break
                if int(lab[q[0],q[1],q[2]])!=inst: break
                n+=1
        t_ax=n+1
        if best is None or t_ax<best: best=t_ax
    return best

acc=[]; skip={"oob":0,"octant":0,"overlap":0,"thick":0}
os.makedirs(OUT+"/control",exist_ok=True)
files=sorted(glob.glob(BLOCKS+"/*/*.npz"))
print("label blocks:",len(files),flush=True)
for f in files:
    if len(acc)>=TARGET: break
    slab=f.split("/")[-2]; tile=f.split("/")[-1][:-4]
    try:
        z=np.load(f); lab=z[z.files[0]]
    except Exception:
        continue
    nz=np.argwhere(lab>0)
    if len(nz)<500: continue
    med=None
    pick=nz[rng.choice(len(nz),size=min(40,len(nz)),replace=False)]
    ths=[run_thickness(lab,tuple(p),int(lab[p[0],p[1],p[2]])) for p in pick]
    med=float(np.median(ths))
    for p,th in zip(pick,ths):
        if len(acc)>=TARGET: break
        if th>med: skip["thick"]+=1; continue
        oz,oy,ox=origin(slab,tile)
        cz,cy,cx=oz+int(p[0]),oy+int(p[1]),ox+int(p[2])
        z0,y0,x0=cz-N//2,cy-N//2,cx-N//2
        if z0<0 or y0<0 or x0<0 or z0+N>Z or y0+N>Y or x0+N>X: skip["oob"]+=1; continue
        ov=False
        for a in acc:
            inter=1.0
            for k,dd in enumerate((abs(z0-a[0]),abs(y0-a[1]),abs(x0-a[2]))): inter*=max(0,N-dd)
            if inter>0.25*N**3: ov=True; break
        if ov: skip["overlap"]+=1; continue
        v=np.asarray(ct0[z0:z0+N,y0:y0+N,x0:x0+N]); h=N//2
        octs=[v[:h,:h,:h],v[:h,:h,h:],v[:h,h:,:h],v[:h,h:,h:],v[h:,:h,:h],v[h:,:h,h:],v[h:,h:,:h],v[h:,h:,h:]]
        if any((o>0).mean()<0.01 for o in octs): skip["octant"]+=1; continue
        tag="ctl_%s_%s_%d_%d_%d"%(slab,tile,p[0],p[1],p[2])
        np.savez_compressed("%s/control/%s.npz"%(OUT,tag),intensity=v,
                            thickness=np.float32(th),tile_median=np.float32(med),
                            site=np.array([cz,cy,cx]))
        acc.append((z0,y0,x0))
        if len(acc)%10==0: print("  control",len(acc),flush=True)

out=dict(amendment="a07c2f86f39b3bd8",accepted=len(acc),target=TARGET,skipped=skip)
json.dump(out,open(OUT+"/control_summary.json","w"),indent=1)
print(json.dumps(out,indent=1))
