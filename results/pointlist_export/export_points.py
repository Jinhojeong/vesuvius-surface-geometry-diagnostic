"""Subsampled (x, y, z, instance_id) export of the stitched PHerc1218 labels.

For Diego's route B (PR #1). Rules, stated so the manifest can repeat them:
  - lattice: every 8th voxel of the global L1 grid (coords divisible by 8)
  - ids: block-local labels mapped through global_table.json; unmapped ids
    dropped and counted
  - overlaps: blocks visited in plain lexicographic order of their full paths,
    which is what sorted(glob(...)) gives and is NOT numeric z0,y0,x0 order
    (z0, z10080, z10304, ..., z1120, ...). First writer wins on the lattice.
    Since overlapping blocks map to the same global id when the merge is right,
    an id disagreement on overlap is a merge fault, measured separately below
    rather than silently resolved, and the visit order only decides the winner
    on those disagreeing voxels.
"""
import glob, gzip, json, re, os
import numpy as np

BLOCKS="/mnt/vesuvius/p1218_full/blocks"
GT="/mnt/vesuvius/p1218_full/global_table.json"
OUT="/mnt/vesuvius/diego_export"
STEP=8; SHAPE=(11624,3797,3797)
os.makedirs(OUT,exist_ok=True)

table=json.load(open(GT))
files=sorted(glob.glob(BLOCKS+"/*/*.npz"))
print("blocks:",len(files),"| table keys:",len(table),flush=True)

def block_key(f):
    z0=f.split("/")[-2][1:]
    m=re.match(r"tile_y(\d+)_x(\d+)\.npz",os.path.basename(f))
    return "z%s/y%s_x%s"%(z0,m.group(1),m.group(2)),int(z0),int(m.group(1)),int(m.group(2))

grid=tuple((s+STEP-1)//STEP for s in SHAPE)
bitmap=np.zeros(grid,dtype=bool)
print("lattice grid:",grid,"bitmap %.0f MB"%(bitmap.nbytes/1e6),flush=True)

rows=0; dropped_unmapped=0; skipped_taken=0; missing_key=0
ids_seen=set()
f_out=gzip.open(OUT+"/points_step8.csv.gz","wt",newline="")
f_out.write("x,y,z,instance_id\n")
for i,f in enumerate(files):
    key,z0,y0,x0=block_key(f)
    src=table.get(key)
    if src is None: missing_key+=1; continue
    lab=np.load(f)["labels"]
    lmax=int(lab.max())
    lut=np.zeros(lmax+1,dtype=np.int64)
    for k,v in src.items():
        if int(k)<=lmax: lut[int(k)]=v
    a,b,c=(-z0)%STEP,(-y0)%STEP,(-x0)%STEP
    sub=lab[a::STEP,b::STEP,c::STEP]
    gids=lut[sub]
    zz,yy,xx=np.nonzero(gids)
    if len(zz):
        gz=(z0+a)//STEP+zz; gy=(y0+b)//STEP+yy; gx=(x0+c)//STEP+xx
        inb=(gz<grid[0])&(gy<grid[1])&(gx<grid[2])
        gz,gy,gx,vals=gz[inb],gy[inb],gx[inb],gids[zz,yy,xx][inb]
        free=~bitmap[gz,gy,gx]
        skipped_taken+=int((~free).sum())
        gz,gy,gx,vals=gz[free],gy[free],gx[free],vals[free]
        bitmap[gz,gy,gx]=True
        dropped_unmapped+=int((sub>0).sum()-len(zz))
        buf="\n".join("%d,%d,%d,%d"%(x*STEP,y*STEP,z*STEP,v)
                      for z,y,x,v in zip(gz.tolist(),gy.tolist(),gx.tolist(),vals.tolist()))
        if buf: f_out.write(buf+"\n")
        rows+=len(vals)
        ids_seen.update(np.unique(vals).tolist())
    if (i+1)%100==0: print("  %d/%d blocks, %.1fM rows"%(i+1,len(files),rows/1e6),flush=True)
f_out.close()

# overlap agreement: same physical lattice point from two blocks should map to
# the same global id. sample 15 y-adjacent pairs.
pairs=[]; import collections
by=collections.defaultdict(dict)
for f in files:
    key,z0,y0,x0=block_key(f); by[(z0,x0)][y0]=f
for (z0,x0),m in by.items():
    ys=sorted(m)
    for j in range(len(ys)-1):
        pairs.append((m[ys[j]],m[ys[j+1]]))
        if len(pairs)>=15: break
    if len(pairs)>=15: break
agree=tot=0
for fa,fb in pairs:
    ka,za,ya,xa=block_key(fa); kb,zb,yb,xb=block_key(fb)
    la=np.load(fa)["labels"]; lb=np.load(fb)["labels"]
    ta,tb=table[ka],table[kb]
    lua=np.zeros(int(la.max())+1,np.int64); [lua.__setitem__(int(k),v) for k,v in ta.items() if int(k)<=la.max()]
    lub=np.zeros(int(lb.max())+1,np.int64); [lub.__setitem__(int(k),v) for k,v in tb.items() if int(k)<=lb.max()]
    ov0,ov1=yb,ya+la.shape[1]          # overlap in global y
    if ov1<=ov0: continue
    ga=lua[la[::16,ov0-ya:ov1-ya:16,::16]]; gb=lub[lb[::16,ov0-yb:ov1-yb:16,::16]]
    m2=(ga>0)&(gb>0)
    tot+=int(m2.sum()); agree+=int((ga[m2]==gb[m2]).sum())

summ=dict(rows=rows,distinct_ids=len(ids_seen),blocks=len(files),
          missing_table_key=missing_key,dropped_unmapped_vox=dropped_unmapped,
          overlap_lattice_repeats=skipped_taken,
          overlap_id_agreement="%d/%d = %.4f"%(agree,tot,agree/max(tot,1)),
          step=STEP,frame="L1 voxels, same grid and axes as the crossing table",
          source="p1218_full/blocks + global_table.json, first-writer-wins on the step-8 lattice")
json.dump(summ,open(OUT+"/export_summary.json","w"),indent=1)
print(json.dumps(summ,indent=1),flush=True)
