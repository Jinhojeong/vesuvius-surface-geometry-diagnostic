import json
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
d=json.load(open("/mnt/vesuvius/mergegap1218/merge_by_gap_v2.json")); rows=d["rows"]
B=[("0-2","<34.6"),("2-4","34.6-69"),("4-6","69-104"),("6-10","104-173"),("10+",">173")]
rng=np.random.default_rng(0)
pool=[];lo=[];hi=[];per=[0.390,0.395,0.393,0.385,0.302]
for b,_ in B:
    m=np.array([r["merged"] for r in rows if r["band"]==b]); t=np.array([r["n_points"] for r in rows if r["band"]==b])
    idx=np.arange(len(m)); bs=[]
    for _ in range(4000):
        s=rng.choice(idx,len(idx)); bs.append(m[s].sum()/max(t[s].sum(),1))
    pool.append(m.sum()/t.sum()); lo.append(np.percentile(bs,2.5)); hi.append(np.percentile(bs,97.5))
x=np.arange(5)
fig,ax=plt.subplots(figsize=(9.5,4.6))
ax.errorbar(x[:4],pool[:4],yerr=[np.array(pool[:4])-lo[:4],np.array(hi[:4])-np.array(pool[:4])],
            fmt="o-",color="#c0392b",lw=2,ms=8,capsize=5,label="merge rate, contact regime")
ax.errorbar(x[4:],pool[4:],yerr=[[pool[4]-lo[4]],[hi[4]-pool[4]]],fmt="s",color="#999",ms=8,capsize=5,
            label="widest band (longer intervals, not comparable)")
ax.plot(x,per,"^--",color="#2980b9",lw=1.5,ms=7,label="per-step positive rate (length-free)")
ax.axhspan(min(lo[:4]),max(hi[:4]),color="#c0392b",alpha=0.07)
ax.set_xticks(x); ax.set_xticklabels([u for _,u in B])
ax.set_xlabel("measured inter-sheet gap, microns")
ax.set_ylabel("rate")
ax.set_title("Does the published m7 surface prediction merge touching sheets?\n"
             "300 PHerc1218 contact crops; flat across the contact regime",fontsize=11,weight="bold")
ax.set_ylim(0,0.5); ax.legend(fontsize=8,loc="upper right"); ax.grid(alpha=0.25)
fig.tight_layout(); fig.savefig("/mnt/vesuvius/mergegap1218/merge_by_gap.png",dpi=130,facecolor="white")
print("saved")
