"""Exact per-band rejection breakdown, replayed from the frozen loop in p11_crops.py.

No CT reads are needed. p11_crops.py tests in a fixed order (target gate, oob,
overlap, read, octant) and skipped["read"] is 0, so for any census row that was
not accepted the reason is determined by: oob if the cube leaves the volume,
else overlap if it overlaps an already-accepted cube by more than a quarter,
else octant. The accepted list is exactly the 254 published crops in census
order, so the overlap test can be replayed as arithmetic.
"""
import glob, json, re
import numpy as np
from collections import Counter, defaultdict

OUT="/mnt/vesuvius/tightgap1218"; BASE="/mnt/vesuvius/kaggle_tightgap1218"
N=128; TARGET=60; Z,Y,X=23247,7593,7593

rows=json.load(open(OUT+"/sites_gaps.json"))
pub={f.split("/")[-1][:-4] for f in glob.glob(BASE+"/crops/*.npz")}
print("census rows %d | published crops %d"%(len(rows),len(pub)),flush=True)

def origin(slab,tile):
    p=tile.split("_"); return int(slab[1:]),int(p[1][1:]),int(p[2][1:])

per_band=Counter(); accepted=[]
reason=defaultdict(Counter)
for r in rows:
    b=r["band"]
    if per_band[b]>=TARGET:
        reason[b]["gated (band already full)"]+=1; continue
    oz,oy,ox=origin(r["slab"],r["tile"])
    cz,cy,cx=oz+r["z"],oy+r["y"],ox+r["x"]
    z0,y0,x0=cz-N//2,cy-N//2,cx-N//2
    tag="%s_%s_%d_%d_%d"%(r["slab"],r["tile"],r["z"],r["y"],r["x"])
    if z0<0 or y0<0 or x0<0 or z0+N>Z or y0+N>Y or x0+N>X:
        reason[b]["oob"]+=1; continue
    ov=False
    for a in accepted:
        inter=1.0
        for k,v in enumerate((abs(z0-a[0]),abs(y0-a[1]),abs(x0-a[2]))): inter*=max(0,N-v)
        if inter>0.25*N**3: ov=True; break
    if ov: reason[b]["overlap"]+=1; continue
    if tag in pub:
        accepted.append((z0,y0,x0)); per_band[b]+=1; reason[b]["accepted"]+=1
    else:
        reason[b]["octant"]+=1

print("\nreplayed realised:",dict(per_band),"total",sum(per_band.values()))
tot=Counter()
print("\n%-8s %8s %8s %8s %8s %8s"%("band","census","accept","octant","overlap","oob"))
for b in ("0-2","2-4","4-6","6-10","10+"):
    c=reason[b]; n=sum(c.values()); tot.update(c)
    print("%-8s %8d %8d %8d %8d %8d"%(b,n,c["accepted"],c["octant"],c["overlap"],c["oob"]))
print("%-8s %8d %8d %8d %8d %8d"%("all",sum(sum(v.values()) for v in reason.values()),
      tot["accepted"],tot["octant"],tot["overlap"],tot["oob"]))
print("\ncrops_summary.json says: accepted 254, octant 15840, overlap 679, oob 24")
c=reason["0-2"]
tried=c["accepted"]+c["octant"]+c["overlap"]+c["oob"]
print("\n0-2 band: %d census sites, all tried (the band never reached the 60 target)."%tried)
print("          octant rejects %d of %d tried = %.1f%%, ships %d."%(c["octant"],tried,100*c["octant"]/tried,c["accepted"]))
