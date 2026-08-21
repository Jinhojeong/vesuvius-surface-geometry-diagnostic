"""Verification for the corrected v3 set, before anything is uploaded or posted.

Every check either PASSes or FAILs out loud. The first three exist because the
defect that caused this rebuild would have been caught by them.
"""
import glob, hashlib, json, os, sys
import numpy as np, zarr

CT=("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/"
    "PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr")
OUT="/mnt/vesuvius/tightgap1218_v3"; N=128
fails=[]
def chk(name,ok,detail=""):
    print(("PASS  " if ok else "FAIL  ")+name+(("  | "+detail) if detail else ""),flush=True)
    if not ok: fails.append(name)

L1=zarr.open(CT,mode="r")["1"]; L0=zarr.open(CT,mode="r")["0"]
crops=sorted(glob.glob(OUT+"/crops/*.npz")); ctls=sorted(glob.glob(OUT+"/control/*.npz"))
summ=json.load(open(OUT+"/crops_summary.json"))
csum=json.load(open(OUT+"/control_summary.json"))
lab_s=json.load(open(OUT+"/labels_summary.json"))
man=[json.loads(l) for l in open(OUT+"/MANIFEST.jsonl")]

print("=== 1. the defect that caused this rebuild ===")
rng=np.random.default_rng(0); samp=[crops[i] for i in rng.choice(len(crops),30,replace=False)]
cors=[]; exact=0
for f in samp:
    d=np.load(f,allow_pickle=True); v=d["intensity"]; m=(d["instance"]>0)
    a=v.astype(np.float64).ravel(); b=m.astype(np.float64).ravel()
    cors.append(float(np.corrcoef(a,b)[0,1]) if a.std()>0 and b.std()>0 else float("nan"))
    sz,sy,sx=[int(q) for q in d["site"]]; z0,y0,x0=sz-N//2,sy-N//2,sx-N//2
    if np.array_equal(v,np.asarray(L1[z0:z0+N,y0:y0+N,x0:x0+N])): exact+=1
cors=np.array(cors)
# alignment is tested by shifting, not by the absolute correlation, which tracks
# how much empty space a crop holds rather than how well CT and labels line up
wins=0;n=0;dmu=[]
for f in samp:
    d=np.load(f,allow_pickle=True); m=(d["instance"]>0)
    if not m.any() or m.all(): continue
    sz,sy,sx=[int(q) for q in d["site"]]; z0,y0,x0=sz-N//2,sy-N//2,sx-N//2
    v=np.asarray(L1[z0:z0+N,y0:y0+N,x0:x0+N])
    rt=float(np.corrcoef(v.astype(np.float64).ravel(),m.astype(np.float64).ravel())[0,1])
    best=-9
    for sh in (-40,-20,20,40):
        if z0+sh<0 or z0+sh+N>L1.shape[0]: continue
        w=np.asarray(L1[z0+sh:z0+sh+N,y0:y0+N,x0:x0+N])
        best=max(best,float(np.corrcoef(w.astype(np.float64).ravel(),m.astype(np.float64).ravel())[0,1]))
    n+=1
    if rt>best: wins+=1
    vf=v.astype(np.float64); dmu.append(float(vf[m].mean()-vf[~m].mean()))
chk("CT at the site beats CT shifted 20 or 40 voxels",wins>=0.8*n,"%d of %d"%(wins,n))
# run this one over EVERY crop rather than a sample; it needs no network
allmu=[]
for f in crops:
    d=np.load(f,allow_pickle=True); m=(d["instance"]>0)
    if not m.any() or m.all(): continue
    vf=d["intensity"].astype(np.float64)
    allmu.append(float(vf[m].mean()-vf[~m].mean()))
allmu=np.array(allmu)
chk("mean CT under the label exceeds the background, all %d crops"%len(allmu),
    float((allmu>0).mean())>0.9,
    "median +%.1f grey, %.1f%% positive, min %.1f"%(np.median(allmu),100*(allmu>0).mean(),allmu.min()))
chk("crop intensity is byte-identical to CT level 1 at its own site",exact==len(samp),"%d/%d"%(exact,len(samp)))
d=np.load(samp[0],allow_pickle=True); sz,sy,sx=[int(q) for q in d["site"]]
z0,y0,x0=sz-N//2,sy-N//2,sx-N//2
wrong=np.asarray(L0[z0:z0+N,y0:y0+N,x0:x0+N])
chk("crop intensity is NOT the old level-0 read",not np.array_equal(d["intensity"],wrong))

print("\n=== 2. membership and manifest ===")
chk("contact file count matches crops_summary",len(crops)==summ["accepted"],"%d vs %d"%(len(crops),summ["accepted"]))
chk("control file count matches control_summary",len(ctls)==csum["accepted"],"%d vs %d"%(len(ctls),csum["accepted"]))
disk={"crops/"+os.path.basename(f) for f in crops}|{"control/"+os.path.basename(f) for f in ctls}
mset={m["file"] for m in man}
chk("manifest and disk agree exactly",disk==mset,"disk %d manifest %d, sym-diff %d"%(len(disk),len(mset),len(disk^mset)))
bad=0
for m in rng.choice(man,25,replace=False):
    h=hashlib.sha256(open(OUT+"/"+m["file"],"rb").read()).hexdigest()
    if h!=m["sha256"]: bad+=1
chk("sha256 matches on a 25-file sample",bad==0,"%d mismatches"%bad)
tags=[os.path.basename(f) for f in crops]
chk("no duplicate crop names",len(set(tags))==len(tags))

print("\n=== 3. the frozen rules actually held ===")
chk("all five bands at the 60 target",all(summ["realised"].get(b)==60 for b in ("0-2","2-4","4-6","6-10","10+")),json.dumps(summ["realised"]))
EDGE={"0-2":(0,2),"2-4":(2,4),"4-6":(4,6),"6-10":(6,10),"10+":(10,1e9)}
off=[m for m in man if m["arm"]=="crops" and not (EDGE[m["band"]][0]<=m["gap"]<=EDGE[m["band"]][1])]
chk("every crop's gap lies inside its stated band",len(off)==0,"%d outside"%len(off))
sites=[[int(v) for v in m["site"]] for m in man if m["arm"]=="crops"]
viol=0
for i in range(len(sites)):
    for j in range(i+1,len(sites)):
        inter=1.0
        for k in range(3): inter*=max(0,N-abs(sites[i][k]-sites[j][k]))
        if inter>0.25*N**3: viol+=1
chk("no two accepted crops overlap by more than a quarter",viol==0,"%d violating pairs"%viol)
oob=[s for s in sites if any(v-N//2<0 for v in s) or s[0]+N//2>L1.shape[0] or s[1]+N//2>L1.shape[1] or s[2]+N//2>L1.shape[2]]
chk("every crop lies wholly inside the level-1 volume",len(oob)==0,"%d out of bounds"%len(oob))

print("\n=== 4. AMENDMENT_2 failure conditions ===")
chk("accepts MORE than the displaced run (else the premise was wrong)",len(crops)>254,"%d vs 254"%len(crops))
chk("alignment holds under the shift test (the amendment's raw-correlation form of this is superseded, see RESULTS)",wins>=0.8*n,
    "%d of %d; raw corr median %.3f is emptiness-driven, not alignment"%(wins,n,float(np.nanmedian(cors))))
chk("band counts are NOT identical to v2 (else the change did not take)",
    summ["realised"]!={"0-2":14,"2-4":60,"4-6":60,"6-10":60,"10+":60})

print("\n=== 5. contents every consumer needs ===")
need={"intensity","instance","surface","gap","band","site","A_id","B_id","id_bases","ct_empty_frac"}
missing=[os.path.basename(f) for f in samp if not need.issubset(set(np.load(f,allow_pickle=True).keys()))]
chk("contact crops carry every promised array",len(missing)==0,str(missing[:3]))
cneed={"intensity","instance","surface","thickness","tile_median","site","ct_empty_frac"}
cmiss=[os.path.basename(f) for f in ctls[:10] if not cneed.issubset(set(np.load(f,allow_pickle=True).keys()))]
chk("control crops carry every promised array",len(cmiss)==0,str(cmiss[:3]))
both=sum(1 for m in man if m["arm"]=="crops" and m.get("both_instances_present"))
chk("both_instances_present recomputes to the summary",both==lab_s["contact"]["both"],"%d vs %d"%(both,lab_s["contact"]["both"]))

print("\n=== summary ===")
print(json.dumps({"contact":len(crops),"control":len(ctls),"bands":summ["realised"],
                  "both_instances":both,"skipped":summ["skipped"],
                  "emptiness_contact":lab_s["contact"]["emptiness"],
                  "ct_corr_median":round(float(np.nanmedian(cors)),3)},indent=1))
print(("\nALL CHECKS PASSED" if not fails else "\nFAILED: "+", ".join(fails)))
sys.exit(1 if fails else 0)
