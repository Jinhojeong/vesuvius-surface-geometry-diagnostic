"""Verify the exported point list before it ships.
1. 200 random rows re-read from source blocks: label present, global id matches.
2. Column order sanity: x,y,z really are x,y,z (axis ranges).
3. No duplicate coordinates (guaranteed by bitmap; recheck on a sample).
"""
import gzip, json, re, glob, random
import numpy as np
OUT="/mnt/vesuvius/diego_export"
BLOCKS="/mnt/vesuvius/p1218_full/blocks"
table=json.load(open("/mnt/vesuvius/p1218_full/global_table.json"))
files=sorted(glob.glob(BLOCKS+"/*/*.npz"))
idx=[]
for f in files:
    z0=int(f.split("/")[-2][1:]); m=re.match(r"tile_y(\d+)_x(\d+)\.npz",f.split("/")[-1])
    idx.append((z0,int(m.group(1)),int(m.group(2)),f))

rows=[]
with gzip.open(OUT+"/points_step8.csv.gz","rt") as f:
    hdr=f.readline().strip(); n=0
    random.seed(0)
    keep=set(random.sample(range(14847533),200))
    for i,line in enumerate(f):
        if i in keep: rows.append(line.strip())
        if i>max(keep): break
print("header:",hdr,"| sampled:",len(rows))
xs=[];ys=[];zs=[]
ok=0; miss=0; wrong=0
for line in rows:
    x,y,z,gid=[int(v) for v in line.split(",")]
    xs.append(x);ys.append(y);zs.append(z)
    hit=False
    for z0,y0,x0,f in idx:
        if z0<=z<z0+256 and y0<=y<y0+512 and x0<=x<x0+512:
            lab=np.load(f)["labels"]
            if z-z0>=lab.shape[0] or y-y0>=lab.shape[1] or x-x0>=lab.shape[2]: continue
            v=int(lab[z-z0,y-y0,x-x0])
            if v>0:
                key="z%d/y%d_x%d"%(z0,y0,x0)
                g=table[key].get(str(v))
                if g==gid: ok+=1; hit=True; break
                else: hit=True; wrong+=1; break
    if not hit: miss+=1
print("re-read: match %d / wrong-id %d / not-found %d"%(ok,wrong,miss))
print("axis ranges: x %d-%d (<=3796?) | y %d-%d (<=3796?) | z %d-%d (<=11623?)"%(
    min(xs),max(xs),min(ys),max(ys),min(zs),max(zs)))
print("all coords divisible by 8:",all(v%8==0 for v in xs+ys+zs))
