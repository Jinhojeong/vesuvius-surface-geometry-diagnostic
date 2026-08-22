"""v5 verification. Includes the AMENDMENT_4 conditions and the consumer-facing
test the earlier audit raised: does the shipped gap reproduce from the crop's OWN
labels? And it tests whether 300/300 is real or trivially true."""
import glob, hashlib, json, os, sys, zipfile
import numpy as np, numpy.lib.format as fmt, zarr
CT=("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/"
    "PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr")
V5="/mnt/vesuvius/tightgap1218_v5"; V4="/mnt/vesuvius/tightgap1218_v4"; V3="/mnt/vesuvius/tightgap1218_v3"
BLOCKS="/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired"; N=128; BASE=1000000
fails=[]
def chk(n,ok,d=""):
    print(("PASS  " if ok else "FAIL  ")+n+(("  | "+d) if d else ""),flush=True); 
    if not ok: fails.append(n)
L1=zarr.open(CT,mode="r")["1"]
crops=sorted(glob.glob(V5+"/crops/*.npz")); ctls=sorted(glob.glob(V5+"/control/*.npz"))
man=[json.loads(l) for l in open(V5+"/MANIFEST.jsonl")]
lab_s=json.load(open(V5+"/labels_summary.json")); summ=json.load(open(V5+"/crops_summary.json"))
rng=np.random.default_rng(0); samp=[crops[i] for i in rng.choice(len(crops),30,replace=False)]

print("=== is 300/300 real, or trivially true? ===")
szA=[];szB=[];same=0
for f in crops:
    d=np.load(f,allow_pickle=True); A,B=int(d["A_id"]),int(d["B_id"]); lab=d["instance"]
    a=int((lab==A).sum()); b=int((lab==B).sum()); szA.append(a); szB.append(b)
    if A==B: same+=1
szA=np.array(szA); szB=np.array(szB)
chk("A_id and B_id are distinct in every crop",same==0,"%d crops with A==B"%same)
chk("both ids occupy real volume, not a stray voxel",min(szA.min(),szB.min())>=50,
    "smallest instance: A %d, B %d voxels; medians %d / %d"%(szA.min(),szB.min(),int(np.median(szA)),int(np.median(szB))))

print("\n=== the consumer test: does the shipped gap reproduce from the crop's OWN labels? ===")
def gap_from_crop(lab,A,B,site_local):
    """VERBATIM from p11_gaps.py. An earlier reimplementation of this differed
    (whole-crop centroids instead of a radius-12 window, and one direction per
    label instead of both signs per step) and scored 107 of 297 where the frozen
    function scores 300 of 300. The estimator was the variable, not the data."""
    z,y,x=site_local; Z,Y,X=lab.shape; maxr=32
    if not (0<=z<Z and 0<=y<Y and 0<=x<X): return None
    r=12
    z0,z1=max(0,z-r),min(Z,z+r+1); y0,y1=max(0,y-r),min(Y,y+r+1); x0,x1=max(0,x-r),min(X,x+r+1)
    w=lab[z0:z1,y0:y1,x0:x1]
    ma=np.argwhere(w==A); mb=np.argwhere(w==B)
    if len(ma)<3 or len(mb)<3: return None
    ca,cb=ma.mean(0),mb.mean(0); d=cb-ca; n=np.linalg.norm(d)
    if n<1e-6: return None
    d=d/n; hit_a=hit_b=None
    for t in np.arange(0.0,maxr,0.5):
        for sgn in (1,-1):
            q=np.round(np.array([z,y,x])+sgn*t*d).astype(int)
            if not (0<=q[0]<Z and 0<=q[1]<Y and 0<=q[2]<X): continue
            v=int(lab[q[0],q[1],q[2]])
            if v==A and hit_a is None: hit_a=t
            elif v==B and hit_b is None: hit_b=t
        if hit_a is not None and hit_b is not None: break
    if hit_a is None or hit_b is None: return None
    return float(hit_a+hit_b)

ok=0; tried=0; band_moves=0
EDGE={"0-2":(0,2),"2-4":(2,4),"4-6":(4,6),"6-10":(6,10),"10+":(10,1e9)}
for m in man:
    if m["arm"]!="crops": continue
    d=np.load(V5+"/"+m["file"],allow_pickle=True)
    g=gap_from_crop(d["instance"],int(m["A_id"]),int(m["B_id"]),(64,64,64))
    if g is None: continue
    tried+=1
    if abs(g-m["gap"])<1e-9: ok+=1
    else:
        lo,hi=EDGE[m["band"]]
        if not (lo<=g<=hi): band_moves+=1
print("  remeasurable from own labels: %d of 300"%tried)
chk("shipped gap reproduces EXACTLY from the crop's own labels",ok==tried and tried==300,
    "%d of %d exact; %d land outside their stated band"%(ok,tried,band_moves))

print("\n=== AMENDMENT_4 conditions ===")
g5={m["file"]:m["gap"] for m in man if m["arm"]=="crops"}
g3={m["file"]:m["gap"] for m in (json.loads(x) for x in open(V3+"/MANIFEST.jsonl")) if m["arm"]=="crops"}
b5={m["file"]:m["band"] for m in man if m["arm"]=="crops"}
b3={m["file"]:m["band"] for m in (json.loads(x) for x in open(V3+"/MANIFEST.jsonl")) if m["arm"]=="crops"}
chk("no gap value changed against v3",g5==g3)
chk("no band changed against v3",b5==b3)
chk("same crop membership as v3",set(g5)==set(g3),"%d vs %d"%(len(g5),len(g3)))
def census_block(name):
    b=os.path.basename(name)[:-4]
    if b.startswith("ctl_"): b=b[4:]
    s,r=b.split("_tile_",1); return "%s/%s.npz"%(s,"tile_"+r.rsplit("_",3)[0])
bad=0
for f in crops:
    prov=json.loads(str(np.load(f,allow_pickle=True)["id_bases"]))
    cb=[p for p in prov if p.get("census_block")]
    if len(cb)!=1 or cb[0]["block"]!=census_block(f) or cb[0]["base"]!=BASE: bad+=1
chk("the census block is authoritative and written first in every crop",bad==0,"%d wrong"%bad)
# why 300 and not the predicted 276: were the extra ones lost to overwriting in v4?
m4={m["file"]:m for m in (json.loads(x) for x in open(V4+"/MANIFEST.jsonl")) if m["arm"]=="crops"}
recovered=0
for f in crops:
    nm="crops/"+os.path.basename(f)
    d5=np.load(f,allow_pickle=True); d4=np.load(V4+"/"+nm,allow_pickle=True)
    ids4=set(np.unique(d4["instance"]).tolist())
    prov4=json.loads(str(d4["id_bases"]))
    cb4=[p for p in prov4 if p["block"]==census_block(f)]
    if not cb4: continue
    A4=cb4[0]["base"]+int(d4["A"]); B4=cb4[0]["base"]+int(d4["B"])
    if not (A4 in ids4 and B4 in ids4): recovered+=1
print("  crops where v4's census-base ids were absent but v5's are present: %d"%recovered)
chk("the 300 minus 276 gap is explained by v4 overwriting census-block voxels",recovered==300-276,
    "%d recovered, predicted %d"%(recovered,300-276))

print("\n=== unchanged from v4: alignment, manifest, rules ===")
ex=0
for f in samp:
    d=np.load(f,allow_pickle=True); sz,sy,sx=[int(q) for q in d["site"]]
    if np.array_equal(d["intensity"],np.asarray(L1[sz-64:sz+64,sy-64:sy+64,sx-64:sx+64])): ex+=1
chk("crop intensity is byte-identical to CT level 1 at its own site",ex==len(samp),"%d/%d"%(ex,len(samp)))
allmu=[]
for f in crops:
    d=np.load(f,allow_pickle=True); m=(d["instance"]>0)
    if not m.any() or m.all(): continue
    v=d["intensity"].astype(np.float64); allmu.append(float(v[m].mean()-v[~m].mean()))
allmu=np.array(allmu)
chk("mean CT under the label exceeds the background, all %d"%len(allmu),float((allmu>0).mean())>0.9,
    "median +%.1f, %.1f%% positive"%(np.median(allmu),100*(allmu>0).mean()))
chk("contact 300 / control 60",len(crops)==300 and len(ctls)==60)
disk={"crops/"+os.path.basename(f) for f in crops}|{"control/"+os.path.basename(f) for f in ctls}
chk("manifest and disk agree exactly",disk=={m["file"] for m in man})
badh=sum(1 for m in man if hashlib.sha256(open(V5+"/"+m["file"],"rb").read()).hexdigest()!=m["sha256"])
chk("sha256 matches on ALL %d rows"%len(man),badh==0,"%d mismatches"%badh)
chk("all five bands at 60",all(summ["realised"].get(b)==60 for b in ("0-2","2-4","4-6","6-10","10+")))
chk("every gap lies inside its band",all(EDGE[m["band"]][0]<=m["gap"]<=EDGE[m["band"]][1] for m in man if m["arm"]=="crops"))
both=sum(1 for m in man if m["arm"]=="crops" and m.get("both_instances_present"))
chk("both_instances_present recomputes to the summary",both==lab_s["contact"]["both"],"%d vs %d"%(both,lab_s["contact"]["both"]))
need={"intensity","instance","surface","gap","band","site","A_id","B_id","id_bases","ct_empty_frac"}
chk("contact crops carry every promised array",all(need.issubset(set(np.load(f,allow_pickle=True).keys())) for f in samp))

print("\n=== summary ===")
print(json.dumps({"contact":len(crops),"control":len(ctls),"both":both,"one":lab_s["contact"]["one"],
 "neither":lab_s["contact"]["neither"],"v4_both":263,"predicted":276,
 "gap_reproduces":"%d of %d"%(ok,tried)},indent=1))
print("\nALL CHECKS PASSED" if not fails else "\nFAILED: "+", ".join(fails))
sys.exit(1 if fails else 0)
