"""Replay the frozen selection rule in one clean pass to fix authoritative membership.

The octant test depends only on CT, so a site that passed it in either earlier run
passed it, and a crop file exists for exactly those sites. Replaying rule order plus
the overlap exclusion against a fresh accepted list therefore reproduces the single
clean pass the frozen document specifies.
"""
import glob, json, os
import numpy as np

SRC="/mnt/vesuvius/tightgap1218"
N=128; TARGET=60
rows=json.load(open(SRC+"/sites_gaps.json"))
present={os.path.basename(f)[:-4] for f in glob.glob(SRC+"/crops/*.npz")}
print("files on disk:", len(present))

def origin(slab,tile):
    p=tile.split("_"); return int(slab[1:]), int(p[1][1:]), int(p[2][1:])

accepted=[]; chosen=[]; per={}
for r in rows:
    b=r["band"]
    if per.get(b,0)>=TARGET: continue
    tag="%s_%s_%d_%d_%d"%(r["slab"],r["tile"],r["z"],r["y"],r["x"])
    if tag not in present: continue          # failed octant/oob in both runs
    oz,oy,ox=origin(r["slab"],r["tile"])
    z0,y0,x0=oz+r["z"]-N//2, oy+r["y"]-N//2, ox+r["x"]-N//2
    ov=False
    for a in accepted:
        inter=1.0
        for k,dd in enumerate((abs(z0-a[0]),abs(y0-a[1]),abs(x0-a[2]))): inter*=max(0,N-dd)
        if inter>0.25*N**3: ov=True; break
    if ov: continue
    accepted.append((z0,y0,x0)); chosen.append(tag); per[b]=per.get(b,0)+1

import collections
cb=collections.Counter(t for t in chosen for _ in [0])
bands={}
for r in rows:
    tag="%s_%s_%d_%d_%d"%(r["slab"],r["tile"],r["z"],r["y"],r["x"])
    if tag in set(chosen): bands[tag]=r["band"]
print("replayed selection:", len(chosen), dict(collections.Counter(bands.values())))
orphans=sorted(present-set(chosen))
print("orphans (on disk, not selected by the rule):", len(orphans))
ob=collections.Counter()
for t in orphans:
    d=np.load(SRC+"/crops/%s.npz"%t, allow_pickle=True); ob[str(d["band"])]+=1
print("orphans by band:", dict(ob))
for t in orphans: print("   ", t)
json.dump(dict(authoritative=sorted(chosen), orphans=orphans,
               by_band=dict(collections.Counter(bands.values()))),
          open(SRC+"/membership_replay.json","w"), indent=1)
