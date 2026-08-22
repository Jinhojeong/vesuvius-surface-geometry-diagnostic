"""v5 labels, per AMENDMENT_4 (0a5d38b6e8bba1e9).

The census block, the one the crop filename names, is authoritative. A and B are
block-local ids in THAT block's numbering, so A_id = its base + A. It is written
first; every other block fills only voxels nothing has written yet.
Block shapes are read from the stored headers rather than assumed.
Crops are copied from v3 (selection untouched); gaps carried over verbatim.
"""
import glob, hashlib, json, os, zipfile
import numpy as np, numpy.lib.format as fmt

SRC="/mnt/vesuvius/tightgap1218_v3"; OUT="/mnt/vesuvius/tightgap1218_v5"
BLOCKS="/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired"
N=128; BASE=1000000

for s in ("crops","control"):
    d="%s/%s"%(OUT,s)
    if os.path.isdir(d) and os.listdir(d): raise SystemExit("REFUSING: %s not empty"%d)
    os.makedirs(d,exist_ok=True)

def hdr(f):
    with zipfile.ZipFile(f) as z:
        with z.open("labels.npy") as g:
            v=fmt.read_magic(g)
            s,_,_=fmt.read_array_header_1_0(g) if v==(1,0) else fmt.read_array_header_2_0(g)
    return s
def og(f):
    s=f.split("/")[-2]; t=f.split("/")[-1][:-4]; p=t.split("_")
    return int(s[1:]),int(p[1][1:]),int(p[2][1:])
idx=[(*og(f),hdr(f),f) for f in sorted(glob.glob(BLOCKS+"/*/*.npz"))]
from collections import Counter
print("blocks indexed with REAL shapes:",len(idx),dict(Counter(t[3] for t in idx)),flush=True)

def census_block(name):
    """the block the census read, encoded in the crop filename"""
    base=os.path.basename(name)[:-4]
    if base.startswith("ctl_"): base=base[4:]
    slab,rest=base.split("_tile_",1)
    tile="tile_"+rest.rsplit("_",3)[0]
    return "%s/%s/%s.npz"%(BLOCKS,slab,tile)

def stitch(z0,y0,x0,cpath):
    hits=[]
    for bz,by,bx,s,f in idx:
        if bz>=z0+N or bz+s[0]<=z0: continue
        if by>=y0+N or by+s[1]<=y0: continue
        if bx>=x0+N or bx+s[2]<=x0: continue
        hits.append((0 if f==cpath else 1,bz,by,bx,s,f))
    hits.sort(key=lambda t:(t[0],t[1],t[2],t[3]))
    out=np.zeros((N,N,N),np.int32); prov=[]; cbase=None; rank=0
    for pri,bz,by,bx,s,f in hits:
        lab=np.load(f)["labels"]
        az0,az1=max(z0,bz),min(z0+N,bz+s[0]); ay0,ay1=max(y0,by),min(y0+N,by+s[1]); ax0,ax1=max(x0,bx),min(x0+N,bx+s[2])
        sub=lab[az0-bz:az1-bz,ay0-by:ay1-by,ax0-bx:ax1-bx].astype(np.int32)
        rank+=1; base=BASE*rank
        w=out[az0-z0:az1-z0,ay0-y0:ay1-y0,ax0-x0:ax1-x0]
        assert w.shape==sub.shape,(w.shape,sub.shape,f)
        m=(sub>0)&(w==0); w[m]=sub[m]+base
        prov.append({"block":"/".join(f.split("/")[-2:]),"base":base,"census_block":bool(pri==0)})
        if pri==0: cbase=base
    return out,prov,cbase

stats={"contact":{"n":0,"both":0,"one":0,"neither":0,"no_census_block":0,"multi_block":0},
       "control":{"n":0,"multi_block":0}}
lf={"contact":[],"control":[]}; empt={"contact":[],"control":[]}
for arm in ("crops","control"):
    key="contact" if arm=="crops" else "control"
    for f in sorted(glob.glob("%s/%s/*.npz"%(SRC,arm))):
        d=dict(np.load(f,allow_pickle=True))
        site=[int(v) for v in d["site"]]; z0,y0,x0=[v-N//2 for v in site]
        lab,prov,cb=stitch(z0,y0,x0,census_block(f))
        d["instance"]=lab; d["surface"]=(lab>0).astype(np.uint8); d["id_bases"]=json.dumps(prov)
        e=float((d["intensity"]==0).mean()); d["ct_empty_frac"]=np.float32(e)
        empt[key].append(e); stats[key]["n"]+=1; lf[key].append(float((lab>0).mean()))
        if len(prov)>1: stats[key]["multi_block"]+=1
        if key=="contact":
            if cb is None:
                stats["contact"]["no_census_block"]+=1; d["A_id"]=-1; d["B_id"]=-1
            else:
                A=cb+int(d["A"]); B=cb+int(d["B"]); d["A_id"]=A; d["B_id"]=B
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
stats["amendment"]="0a5d38b6e8bba1e9"; stats["ct_level"]=1
json.dump(stats,open(OUT+"/labels_summary.json","w"),indent=1)
print(json.dumps(stats,indent=1),flush=True)

import shutil
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
