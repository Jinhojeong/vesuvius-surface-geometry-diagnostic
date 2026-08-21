"""v4 labels, per AMENDMENT_3 (db61e452e4cce9ca).

The repair blocks overlap (256 arrays on 224/448 strides), and the old stitch let
a later block overwrite an earlier one, which could relabel the split pair out of
existence. Here the site block is written FIRST and every other block fills only
voxels nothing has written yet.

Crops are copied from the v3 tree, so selection is untouched. Gap values are
carried over verbatim and never recomputed.
"""
import glob, hashlib, json, os, shutil
import numpy as np

SRC="/mnt/vesuvius/tightgap1218_v3"; OUT="/mnt/vesuvius/tightgap1218_v4"
BLOCKS="/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired"
N=128; BASE=1000000

for sub in ("crops","control"):
    d="%s/%s"%(OUT,sub)
    if os.path.isdir(d) and os.listdir(d): raise SystemExit("REFUSING: %s not empty"%d)
    os.makedirs(d,exist_ok=True)

def origin_from_path(f):
    slab=f.split("/")[-2]; tile=f.split("/")[-1][:-4]
    p=tile.split("_"); return int(slab[1:]),int(p[1][1:]),int(p[2][1:])
idx=[(*origin_from_path(f),(256,512,512),f) for f in sorted(glob.glob(BLOCKS+"/*/*.npz"))]
print("blocks indexed from filenames:",len(idx),flush=True)

def stitch(z0,y0,x0,site):
    """Site block first, then others only where nothing is written. AMENDMENT_3."""
    hits=[]
    for bz,by,bx,s,f in idx:
        if bz>=z0+N or bz+s[0]<=z0: continue
        if by>=y0+N or by+s[1]<=y0: continue
        if bx>=x0+N or bx+s[2]<=x0: continue
        sz,sy,sx=site
        is_site=(bz<=sz<bz+s[0] and by<=sy<by+s[1] and bx<=sx<bx+s[2])
        hits.append((0 if is_site else 1,bz,by,bx,s,f,is_site))
    hits.sort(key=lambda t:(t[0],t[1],t[2],t[3]))     # site block first, then deterministic order
    # the site can fall inside more than one block because blocks overlap. Exactly
    # one of them is authoritative: the first in this deterministic order, which is
    # also the one written first, so its base is the one A_id and B_id are built on.
    out=np.zeros((N,N,N),np.int32); prov=[]; site_base=None; rank=0
    for i,(_,bz,by,bx,s,f,is_site) in enumerate(hits):
        lab=np.load(f)["labels"]
        az0,az1=max(z0,bz),min(z0+N,bz+s[0]); ay0,ay1=max(y0,by),min(y0+N,by+s[1]); ax0,ax1=max(x0,bx),min(x0+N,bx+s[2])
        sub=lab[az0-bz:az1-bz,ay0-by:ay1-by,ax0-bx:ax1-bx].astype(np.int32)
        rank+=1; base=BASE*rank
        w=out[az0-z0:az1-z0,ay0-y0:ay1-y0,ax0-x0:ax1-x0]
        m=(sub>0)&(w==0)                                   # <- the fix
        w[m]=sub[m]+base
        prov.append({"block":"/".join(f.split("/")[-2:]),"base":base,
                     "site_block":bool(is_site),"authoritative":bool(i==0 and is_site)})
        if i==0 and is_site: site_base=base
    return out,prov,site_base

stats={"contact":{"n":0,"both":0,"one":0,"neither":0,"no_site_block":0,"multi_block":0},
       "control":{"n":0,"multi_block":0}}
lf={"contact":[],"control":[]}; empt={"contact":[],"control":[]}; gapdiff=0
for arm in ("crops","control"):
    key="contact" if arm=="crops" else "control"
    for f in sorted(glob.glob("%s/%s/*.npz"%(SRC,arm))):
        d=dict(np.load(f,allow_pickle=True))
        site=[int(v) for v in d["site"]]; z0,y0,x0=[v-N//2 for v in site]
        lab,prov,sb=stitch(z0,y0,x0,site)
        d["instance"]=lab; d["surface"]=(lab>0).astype(np.uint8); d["id_bases"]=json.dumps(prov)
        e=float((d["intensity"]==0).mean()); d["ct_empty_frac"]=np.float32(e)
        empt[key].append(e); stats[key]["n"]+=1; lf[key].append(float((lab>0).mean()))
        if len(prov)>1: stats[key]["multi_block"]+=1
        if key=="contact":
            if sb is None:
                stats["contact"]["no_site_block"]+=1; d["A_id"]=-1; d["B_id"]=-1
            else:
                A=sb+int(d["A"]); B=sb+int(d["B"]); d["A_id"]=A; d["B_id"]=B
                ids=set(np.unique(lab).tolist())
                if A in ids and B in ids: stats["contact"]["both"]+=1
                elif A in ids or B in ids: stats["contact"]["one"]+=1
                else: stats["contact"]["neither"]+=1
        np.savez_compressed("%s/%s/%s"%(OUT,arm,os.path.basename(f)),**d)
    print("done",arm,flush=True)

for k in lf:
    stats[k]["labelled_frac_median"]=round(float(np.median(lf[k])),4)
    a=np.array(empt[k])
    stats[k]["emptiness"]=dict(median=round(float(np.median(a)),4),p90=round(float(np.percentile(a,90)),4),
                               frac_over_10pct=round(float((a>0.10).mean()),4),
                               frac_over_25pct=round(float((a>0.25).mean()),4))
stats["amendment"]="db61e452e4cce9ca"; stats["ct_level"]=1
json.dump(stats,open(OUT+"/labels_summary.json","w"),indent=1)
print(json.dumps(stats,indent=1),flush=True)

for j in ("crops_summary.json","control_summary.json"): shutil.copy(SRC+"/"+j,OUT+"/"+j)
with open(OUT+"/MANIFEST.jsonl","w") as mf:
    for arm in ("crops","control"):
        for f in sorted(glob.glob("%s/%s/*.npz"%(OUT,arm))):
            d=np.load(f,allow_pickle=True)
            h=hashlib.sha256(open(f,"rb").read()).hexdigest()
            row=dict(arm=arm,file="%s/%s"%(arm,os.path.basename(f)),sha256=h,
                     site=[int(v) for v in d["site"]],ct_empty_frac=round(float(d["ct_empty_frac"]),5))
            if arm=="crops":
                ids=set(np.unique(d["instance"]).tolist())
                row.update(band=str(d["band"]),gap=float(d["gap"]),A_id=int(d["A_id"]),B_id=int(d["B_id"]),
                           both_instances_present=bool(int(d["A_id"])>0 and int(d["A_id"]) in ids and int(d["B_id"]) in ids))
            else:
                row.update(thickness=float(d["thickness"]),tile_median=float(d["tile_median"]))
            mf.write(json.dumps(row)+"\n")
print("manifest written",flush=True)
