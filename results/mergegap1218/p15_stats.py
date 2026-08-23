"""Uncertainty on the band merge rates, and the preregistered flat-curve question."""
import json
import numpy as np
d=json.load(open("/mnt/vesuvius/mergegap1218/merge_by_gap_v2.json"))
rows=d["rows"]
BANDS=[("0-2","under 34.6"),("2-4","34.6-69.1"),("4-6","69.1-103.7"),
       ("6-10","103.7-172.8"),("10+","above 172.8")]
rng=np.random.default_rng(0)
print("%-8s %-14s %5s %8s %-18s %8s %-18s"%("band","microns","n","median","median 95% CI","pooled","pooled 95% CI"))
med={}; poo={}
for b,um in BANDS:
    v=np.array([r["merge_frac"] for r in rows if r["band"]==b])
    m=np.array([r["merged"] for r in rows if r["band"]==b])
    t=np.array([r["n_points"] for r in rows if r["band"]==b])
    bm=[float(np.median(rng.choice(v,len(v)))) for _ in range(4000)]
    idx=np.arange(len(v))
    bp=[]
    for _ in range(4000):
        s=rng.choice(idx,len(idx))
        bp.append(m[s].sum()/max(t[s].sum(),1))
    med[b]=(float(np.median(v)),np.percentile(bm,2.5),np.percentile(bm,97.5))
    poo[b]=(m.sum()/t.sum(),np.percentile(bp,2.5),np.percentile(bp,97.5))
    print("%-8s %-14s %5d %8.3f  [%.3f, %.3f]      %8.3f  [%.3f, %.3f]"%(
        b,um,len(v),med[b][0],med[b][1],med[b][2],poo[b][0],poo[b][1],poo[b][2]))

print("\nPREREG Q: is the curve flat across all five bands?")
lo=min(poo[b][0] for b,_ in BANDS); hi=max(poo[b][0] for b,_ in BANDS)
print("  pooled range %.3f to %.3f, spread %.3f"%(lo,hi,hi-lo))
print("  do the tightest four overlap each other?")
four=[b for b,_ in BANDS[:4]]
print("   ", " ".join("%s[%.3f,%.3f]"%(b,poo[b][1],poo[b][2]) for b in four))
allpair=all(poo[a][1]<=poo[c][2] and poo[c][1]<=poo[a][2] for a in four for c in four)
print("    every pair of the tightest four overlaps:",allpair)
print("  loosest band vs the tightest four:")
for b in four:
    ov = poo["10+"][1]<=poo[b][2] and poo[b][1]<=poo["10+"][2]
    print("    10+ [%.3f,%.3f] vs %s [%.3f,%.3f] -> %s"%(poo["10+"][1],poo["10+"][2],b,poo[b][1],poo[b][2],"overlap" if ov else "SEPARATE"))

print("\nmonotone-decreasing prediction (my prereg expectation)?")
seq=[poo[b][0] for b,_ in BANDS]
print("  pooled sequence:", " ".join("%.3f"%x for x in seq))
print("  monotone decreasing:", all(seq[i]>=seq[i+1] for i in range(4)))
# gap in microns vs merge, point-biserial style over crops
g=np.array([r["gap_um"] for r in rows]); mf=np.array([r["merge_frac"] for r in rows])
from math import sqrt
r_=np.corrcoef(g,mf)[0,1]
bs=[np.corrcoef(g[s],mf[s])[0,1] for s in (rng.choice(len(g),len(g)) for _ in range(4000))]
print("\ncorrelation of per-crop merge fraction with gap in microns: r=%.3f  95%% CI [%.3f, %.3f]"%(
    r_,np.percentile(bs,2.5),np.percentile(bs,97.5)))
print("crops per band:", {b:int((np.array([r["band"] for r in rows])==b).sum()) for b,_ in BANDS})
