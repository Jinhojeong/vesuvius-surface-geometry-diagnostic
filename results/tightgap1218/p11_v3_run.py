"""Corrected extraction, per AMENDMENT_2 (f9d57773c715cb65).

One change against the published run: the CT is read from level 1, the grid the
sites and labels are given in. Every other rule is byte-identical, including the
per-band target, the octant threshold, the overlap fraction and census order.

The output directory is required to be absent or empty, which is the structural
fix for the pilot-orphan defect aviad12g found in version 1.
"""
import glob, json, os, sys
import numpy as np, zarr

CT=("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/"
    "PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr")
SRC="/mnt/vesuvius/tightgap1218"          # census lives here, read only
OUT="/mnt/vesuvius/tightgap1218_v3"
BLOCKS="/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired"
N=128; TARGET=60; BASEID=1000000
rng=np.random.default_rng(0)

for sub in ("crops","control"):
    d="%s/%s"%(OUT,sub)
    if os.path.isdir(d) and os.listdir(d):
        sys.exit("REFUSING: %s is not empty. AMENDMENT_2 requires a clean output dir."%d)
    os.makedirs(d,exist_ok=True)

ct1=zarr.open(CT,mode="r")["1"]; Z,Y,X=ct1.shape
print("CT level 1",ct1.shape,"(the published run used level 0, 23247x7593x7593)",flush=True)

def origin(slab,tile):
    p=tile.split("_"); return int(slab[1:]),int(p[1][1:]),int(p[2][1:])
def octant_ok(v):
    h=N//2
    o=[v[:h,:h,:h],v[:h,:h,h:],v[:h,h:,:h],v[:h,h:,h:],
       v[h:,:h,:h],v[h:,:h,h:],v[h:,h:,:h],v[h:,h:,h:]]
    return not any((q>0).mean()<0.01 for q in o)
def overlaps(z0,y0,x0,acc):
    for a in acc:
        i=1.0
        for k,dd in enumerate((abs(z0-a[0]),abs(y0-a[1]),abs(x0-a[2]))): i*=max(0,N-dd)
        if i>0.25*N**3: return True
    return False

# ---------- stage 1: contact arm ----------
rows=json.load(open(SRC+"/sites_gaps.json"))
acc=[]; per={}; skip={"oob":0,"octant":0,"overlap":0,"read":0}
for r in rows:
    b=r["band"]
    if per.get(b,0)>=TARGET: continue
    oz,oy,ox=origin(r["slab"],r["tile"])
    cz,cy,cx=oz+r["z"],oy+r["y"],ox+r["x"]
    z0,y0,x0=cz-N//2,cy-N//2,cx-N//2
    if z0<0 or y0<0 or x0<0 or z0+N>Z or y0+N>Y or x0+N>X: skip["oob"]+=1; continue
    if overlaps(z0,y0,x0,acc): skip["overlap"]+=1; continue
    try: vol=np.asarray(ct1[z0:z0+N,y0:y0+N,x0:x0+N])
    except Exception: skip["read"]+=1; continue
    if not octant_ok(vol): skip["octant"]+=1; continue
    tag="%s_%s_%d_%d_%d"%(r["slab"],r["tile"],r["z"],r["y"],r["x"])
    np.savez_compressed("%s/crops/%s.npz"%(OUT,tag),intensity=vol,
                        gap=np.float32(r["gap"]),band=b,
                        site=np.array([cz,cy,cx]),A=r["A"],B=r["B"])
    acc.append((z0,y0,x0)); per[b]=per.get(b,0)+1
    if len(acc)%25==0: print("  accepted %d | %s"%(len(acc),per),flush=True)
json.dump(dict(amendment="f9d57773c715cb65",ct_level=1,target_per_band=TARGET,
               realised=per,accepted=len(acc),skipped=skip),
          open(OUT+"/crops_summary.json","w"),indent=1)
print("CONTACT",json.dumps(per),"accepted",len(acc),"skipped",skip,flush=True)

# ---------- stage 2: control arm ----------
def run_thickness(lab,p,inst):
    z,y,x=p; best=None
    for ax in range(3):
        n=0
        for sgn in (1,-1):
            q=[z,y,x]
            for t in range(1,40):
                q[ax]=[z,y,x][ax]+sgn*t
                if not(0<=q[0]<lab.shape[0] and 0<=q[1]<lab.shape[1] and 0<=q[2]<lab.shape[2]): break
                if int(lab[q[0],q[1],q[2]])!=inst: break
                n+=1
        t_ax=n+1
        if best is None or t_ax<best: best=t_ax
    return best

cacc=[]; cskip={"oob":0,"octant":0,"overlap":0,"thick":0}
for f in sorted(glob.glob(BLOCKS+"/*/*.npz")):
    if len(cacc)>=TARGET: break
    slab=f.split("/")[-2]; tile=f.split("/")[-1][:-4]
    try:
        z=np.load(f); lab=z[z.files[0]]
    except Exception: continue
    nz=np.argwhere(lab>0)
    if len(nz)<500: continue
    pick=nz[rng.choice(len(nz),size=min(40,len(nz)),replace=False)]
    ths=[run_thickness(lab,tuple(p),int(lab[p[0],p[1],p[2]])) for p in pick]
    med=float(np.median(ths))
    for p,th in zip(pick,ths):
        if len(cacc)>=TARGET: break
        if th>med: cskip["thick"]+=1; continue
        oz,oy,ox=origin(slab,tile)
        cz,cy,cx=oz+int(p[0]),oy+int(p[1]),ox+int(p[2])
        z0,y0,x0=cz-N//2,cy-N//2,cx-N//2
        if z0<0 or y0<0 or x0<0 or z0+N>Z or y0+N>Y or x0+N>X: cskip["oob"]+=1; continue
        if overlaps(z0,y0,x0,cacc): cskip["overlap"]+=1; continue
        v=np.asarray(ct1[z0:z0+N,y0:y0+N,x0:x0+N])
        if not octant_ok(v): cskip["octant"]+=1; continue
        np.savez_compressed("%s/control/ctl_%s_%s_%d_%d_%d.npz"%(OUT,slab,tile,p[0],p[1],p[2]),
                            intensity=v,thickness=np.float32(th),tile_median=np.float32(med),
                            site=np.array([cz,cy,cx]))
        cacc.append((z0,y0,x0))
        if len(cacc)%10==0: print("  control",len(cacc),flush=True)
json.dump(dict(amendment="f9d57773c715cb65",ct_level=1,accepted=len(cacc),
               target=TARGET,skipped=cskip),open(OUT+"/control_summary.json","w"),indent=1)
print("CONTROL accepted",len(cacc),cskip,flush=True)
print("STAGE 1-2 DONE",flush=True)
