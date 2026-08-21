"""v3 stage 3: stitch labels with per-block id bases, add emptiness, write the manifest.
Identical to p11_labels2.py and p11_final.py except for the output root."""
import glob, hashlib, json, os
import numpy as np

OUT="/mnt/vesuvius/tightgap1218_v3"
BLOCKS="/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired"
N=128; BASE=1000000

idx=[]
for f in sorted(glob.glob(BLOCKS+"/*/*.npz")):
    z=np.load(f); idx.append((int(z["z0"]),int(z["y0"]),int(z["x0"]),z["labels"].shape,f))
print("blocks indexed:",len(idx),flush=True)

def stitch(z0,y0,x0,site):
    out=np.zeros((N,N,N),np.int32); prov=[]; site_base=None; rank=0
    for bz,by,bx,sh,f in idx:
        if bz>=z0+N or bz+sh[0]<=z0: continue
        if by>=y0+N or by+sh[1]<=y0: continue
        if bx>=x0+N or bx+sh[2]<=x0: continue
        lab=np.load(f)["labels"]
        az0,az1=max(z0,bz),min(z0+N,bz+sh[0])
        ay0,ay1=max(y0,by),min(y0+N,by+sh[1])
        ax0,ax1=max(x0,bx),min(x0+N,bx+sh[2])
        sub=lab[az0-bz:az1-bz,ay0-by:ay1-by,ax0-bx:ax1-bx].astype(np.int32)
        rank+=1; base=BASE*rank; m=sub>0
        out[az0-z0:az1-z0,ay0-y0:ay1-y0,ax0-x0:ax1-x0][m]=sub[m]+base
        prov.append({"block":"/".join(f.split("/")[-2:]),"base":base})
        sz,sy,sx=site
        if bz<=sz<bz+sh[0] and by<=sy<by+sh[1] and bx<=sx<bx+sh[2]: site_base=base
    return out,prov,site_base

stats={"contact":{"n":0,"both":0,"one":0,"neither":0,"no_site_block":0,"multi_block":0},
       "control":{"n":0,"multi_block":0}}
lf={"contact":[],"control":[]}; empt={"contact":[],"control":[]}
for arm in ("crops","control"):
    key="contact" if arm=="crops" else "control"
    for f in sorted(glob.glob("%s/%s/*.npz"%(OUT,arm))):
        d=dict(np.load(f,allow_pickle=True))
        site=[int(v) for v in d["site"]]; z0,y0,x0=[v-N//2 for v in site]
        lab,prov,sb=stitch(z0,y0,x0,site)
        d["instance"]=lab; d["surface"]=(lab>0).astype(np.uint8)
        d["id_bases"]=json.dumps(prov)
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
        np.savez_compressed(f,**d)
    print("done",arm,flush=True)

for k in lf:
    stats[k]["labelled_frac_median"]=round(float(np.median(lf[k])),4)
    a=np.array(empt[k])
    stats[k]["emptiness"]=dict(median=round(float(np.median(a)),4),p90=round(float(np.percentile(a,90)),4),
                               frac_over_10pct=round(float((a>0.10).mean()),4),
                               frac_over_25pct=round(float((a>0.25).mean()),4))
stats["amendment"]="f9d57773c715cb65"; stats["ct_level"]=1
json.dump(stats,open(OUT+"/labels_summary.json","w"),indent=1)
print(json.dumps(stats,indent=1),flush=True)

with open(OUT+"/MANIFEST.jsonl","w") as mf:
    for arm in ("crops","control"):
        for f in sorted(glob.glob("%s/%s/*.npz"%(OUT,arm))):
            d=np.load(f,allow_pickle=True)
            h=hashlib.sha256(open(f,"rb").read()).hexdigest()
            row=dict(arm=arm,file="%s/%s"%(arm,os.path.basename(f)),sha256=h,
                     site=[int(v) for v in d["site"]],ct_empty_frac=round(float(d["ct_empty_frac"]),5))
            if arm=="crops":
                row.update(band=str(d["band"]),gap=float(d["gap"]),A_id=int(d["A_id"]),B_id=int(d["B_id"]),
                           both_instances_present=bool(int(d["A_id"])>0 and
                               set([int(d["A_id"]),int(d["B_id"])]).issubset(set(np.unique(d["instance"]).tolist()))))
            else:
                row.update(thickness=float(d["thickness"]),tile_median=float(d["tile_median"]))
            mf.write(json.dumps(row)+"\n")
print("manifest written",flush=True)
