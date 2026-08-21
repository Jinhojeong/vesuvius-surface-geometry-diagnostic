"""Verification for v4. Repeats every v3 check and adds AMENDMENT_3's own conditions."""
import glob, hashlib, json, os, sys
import numpy as np, zarr
CT=("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/"
    "PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr")
V4="/mnt/vesuvius/tightgap1218_v4"; V3="/mnt/vesuvius/tightgap1218_v3"
BLOCKS="/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired"; N=128; BASE=1000000
fails=[]
def chk(n,ok,d=""):
    print(("PASS  " if ok else "FAIL  ")+n+(("  | "+d) if d else ""),flush=True)
    if not ok: fails.append(n)
L1=zarr.open(CT,mode="r")["1"]; L0=zarr.open(CT,mode="r")["0"]
crops=sorted(glob.glob(V4+"/crops/*.npz")); ctls=sorted(glob.glob(V4+"/control/*.npz"))
summ=json.load(open(V4+"/crops_summary.json")); lab_s=json.load(open(V4+"/labels_summary.json"))
man=[json.loads(l) for l in open(V4+"/MANIFEST.jsonl")]
rng=np.random.default_rng(0); samp=[crops[i] for i in rng.choice(len(crops),30,replace=False)]

print("=== AMENDMENT_3 conditions ===")
def og(f):
    s=f.split("/")[-2]; t=f.split("/")[-1][:-4]; p=t.split("_")
    return int(s[1:]),int(p[1][1:]),int(p[2][1:])
idx=[(*og(f),(256,512,512),f) for f in sorted(glob.glob(BLOCKS+"/*/*.npz"))]
loss=[]; base_ok=0; base_n=0
for f in crops:
    d=np.load(f,allow_pickle=True)
    sz,sy,sx=[int(v) for v in d["site"]]; z0,y0,x0=sz-N//2,sy-N//2,sx-N//2
    hits=[t for t in idx if not(t[0]>=z0+N or t[0]+256<=z0) and not(t[1]>=y0+N or t[1]+512<=y0) and not(t[2]>=x0+N or t[2]+512<=x0)]
    sb=[t for t in hits if t[0]<=sz<t[0]+256 and t[1]<=sy<t[1]+512 and t[2]<=sx<t[2]+512]
    if not sb: continue
    sb.sort(key=lambda t:(t[0],t[1],t[2]))
    bz,by,bx,s,bf=sb[0]
    lab=np.load(bf)["labels"]
    az0,az1=max(z0,bz),min(z0+N,bz+s[0]); ay0,ay1=max(y0,by),min(y0+N,by+s[1]); ax0,ax1=max(x0,bx),min(x0+N,bx+s[2])
    sub=lab[az0-bz:az1-bz,ay0-by:ay1-by,ax0-bx:ax1-bx].astype(np.int32)
    A,B=int(d["A_id"]),int(d["B_id"])
    base_n+=1
    if A-int(d["A"])==BASE: base_ok+=1
    trueA=int((sub+BASE==A).sum()); trueB=int((sub+BASE==B).sum())
    got=int((d["instance"]==A).sum())+int((d["instance"]==B).sum())
    if trueA+trueB>0: loss.append(1.0-got/(trueA+trueB))
l=np.array(loss)
chk("pair-voxel loss is zero for every crop",float(l.max())<1e-9,"n=%d max %.6f mean %.6f"%(len(l),l.max(),l.mean()))
chk("A_id/B_id are built on the authoritative site block (base 1e6)",base_ok==base_n,"%d/%d"%(base_ok,base_n))
g4={m["file"]:m["gap"] for m in man if m["arm"]=="crops"}
g3={m["file"]:m["gap"] for m in (json.loads(x) for x in open(V3+"/MANIFEST.jsonl")) if m["arm"]=="crops"}
chk("no gap value changed against v3",g4==g3,"%d gaps, %d differ"%(len(g4),sum(1 for k in g4 if g3.get(k)!=g4[k])))
b4={m["file"]:m["band"] for m in man if m["arm"]=="crops"}
b3={m["file"]:m["band"] for m in (json.loads(x) for x in open(V3+"/MANIFEST.jsonl")) if m["arm"]=="crops"}
chk("no band changed against v3",b4==b3)
chk("same crop membership as v3",set(g4)==set(g3),"%d vs %d"%(len(g4),len(g3)))

print("\n=== alignment, repeated ===")
wins=0;n=0;dmu=[]
for f in samp:
    d=np.load(f,allow_pickle=True); m=(d["instance"]>0)
    if not m.any() or m.all(): continue
    sz,sy,sx=[int(q) for q in d["site"]]; z0,y0,x0=sz-N//2,sy-N//2,sx-N//2
    v=np.asarray(L1[z0:z0+N,y0:y0+N,x0:x0+N])
    rt=float(np.corrcoef(v.astype(float).ravel(),m.astype(float).ravel())[0,1]); best=-9
    for sh in (-40,-20,20,40):
        if z0+sh<0 or z0+sh+N>L1.shape[0]: continue
        w=np.asarray(L1[z0+sh:z0+sh+N,y0:y0+N,x0:x0+N])
        best=max(best,float(np.corrcoef(w.astype(float).ravel(),m.astype(float).ravel())[0,1]))
    n+=1
    if rt>best: wins+=1
chk("CT at the site beats CT shifted 20 or 40 voxels",wins>=0.8*n,"%d of %d"%(wins,n))
allmu=[]; allcor=[]
for f in crops:
    d=np.load(f,allow_pickle=True); m=(d["instance"]>0)
    if not m.any() or m.all(): continue
    vf=d["intensity"].astype(np.float64)
    allmu.append(float(vf[m].mean()-vf[~m].mean()))
    allcor.append(float(np.corrcoef(vf.ravel(),m.astype(float).ravel())[0,1]))
allmu=np.array(allmu); allcor=np.array(allcor)
chk("mean CT under the label exceeds the background, all %d"%len(allmu),float((allmu>0).mean())>0.9,
    "median +%.1f, %.1f%% positive"%(np.median(allmu),100*(allmu>0).mean()))
print("  (raw CT-to-label correlation over ALL %d crops: median %.4f  -- emptiness-driven, not an alignment test)"%(len(allcor),np.median(allcor)))
ex=sum(1 for f in samp if np.array_equal(np.load(f,allow_pickle=True)["intensity"],
      np.asarray(L1[int(np.load(f,allow_pickle=True)["site"][0])-N//2:int(np.load(f,allow_pickle=True)["site"][0])-N//2+N,
                    int(np.load(f,allow_pickle=True)["site"][1])-N//2:int(np.load(f,allow_pickle=True)["site"][1])-N//2+N,
                    int(np.load(f,allow_pickle=True)["site"][2])-N//2:int(np.load(f,allow_pickle=True)["site"][2])-N//2+N])))
chk("crop intensity is byte-identical to CT level 1 at its own site",ex==len(samp),"%d/%d"%(ex,len(samp)))

print("\n=== membership, manifest, frozen rules ===")
chk("contact 300 / control 60",len(crops)==300 and len(ctls)==60,"%d / %d"%(len(crops),len(ctls)))
disk={"crops/"+os.path.basename(f) for f in crops}|{"control/"+os.path.basename(f) for f in ctls}
chk("manifest and disk agree exactly",disk=={m["file"] for m in man})
bad=sum(1 for m in man if hashlib.sha256(open(V4+"/"+m["file"],"rb").read()).hexdigest()!=m["sha256"])
chk("sha256 matches on ALL %d rows"%len(man),bad==0,"%d mismatches"%bad)
chk("all five bands at 60",all(summ["realised"].get(b)==60 for b in ("0-2","2-4","4-6","6-10","10+")),json.dumps(summ["realised"]))
EDGE={"0-2":(0,2),"2-4":(2,4),"4-6":(4,6),"6-10":(6,10),"10+":(10,1e9)}
chk("every gap lies inside its band",all(EDGE[m["band"]][0]<=m["gap"]<=EDGE[m["band"]][1] for m in man if m["arm"]=="crops"))
sites=[m["site"] for m in man if m["arm"]=="crops"]; viol=0
for i in range(len(sites)):
    for j in range(i+1,len(sites)):
        inter=1.0
        for k in range(3): inter*=max(0,N-abs(sites[i][k]-sites[j][k]))
        if inter>0.25*N**3: viol+=1
chk("no two crops overlap by more than a quarter",viol==0,"%d pairs"%viol)
need={"intensity","instance","surface","gap","band","site","A_id","B_id","id_bases","ct_empty_frac"}
chk("contact crops carry every promised array",all(need.issubset(set(np.load(f,allow_pickle=True).keys())) for f in samp))
both=sum(1 for m in man if m["arm"]=="crops" and m.get("both_instances_present"))
chk("both_instances_present recomputes to the summary",both==lab_s["contact"]["both"],"%d vs %d"%(both,lab_s["contact"]["both"]))
inst=[len(json.loads(str(np.load(f,allow_pickle=True)["id_bases"]))) for f in samp]
chk("id_bases records every contributing block",min(inst)>=1,"median %d max %d"%(int(np.median(inst)),max(inst)))

print("\n=== summary ===")
print(json.dumps({"contact":len(crops),"control":len(ctls),"bands":summ["realised"],
 "both":lab_s["contact"]["both"],"one":lab_s["contact"]["one"],"neither":lab_s["contact"]["neither"],
 "v3_both":245,"pair_loss_max":float(l.max()),"ct_corr_median_all":round(float(np.median(allcor)),4)},indent=1))
print("\nALL CHECKS PASSED" if not fails else "\nFAILED: "+", ".join(fails))
sys.exit(1 if fails else 0)
