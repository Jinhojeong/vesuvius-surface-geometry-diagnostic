"""Bottom-end collapse audit. Does the histogram LUT pile mass at 0 the way the
affine arm's clip does, and is that pile faithful to real non-located volumes?"""
import json, statistics as st
from pathlib import Path
import numpy as np, tifffile
D=Path("/mnt/vesuvius/experiments/histmatch_confound"); IM=Path("/mnt/vesuvius/kaggle892/images")
SIZE=256; LO,HI=0.0,212.0
rows={}
for f in ("loc60_4arm.jsonl","loc60_affine.jsonl"):
    for ln in (D/f).read_text().splitlines():
        if ln.strip():
            r=json.loads(ln); rows.setdefault(r["sample"],{}).update(r)
luts=json.loads(Path("/mnt/vesuvius/kaggle892/histmatch_check.json").read_text())["luts"]
def blk(nm):
    ct=np.asarray(tifffile.imread(str(IM/f"{nm}.tif"))); o=(ct.shape[0]-SIZE)//2
    return ct[o:o+SIZE,o:o+SIZE,o:o+SIZE], ct
out=[]
for s in sorted(json.loads((D/"located60.json").read_text())):
    b,ct=blk(s); n=b.size; r=rows[s]
    L=np.asarray(luts[s],dtype=np.uint8); m=L[ct]; o=(ct.shape[0]-SIZE)//2
    mb=m[o:o+SIZE,o:o+SIZE,o:o+SIZE]
    aff=b.astype(np.float32)*r["affine_a"]+r["affine_b"]
    affc=np.clip(aff,LO,HI)
    out.append({"sample":s,
      "orig_at0":float((b==0).sum())/n,
      "orig_le2":float((b<=2).sum())/n,
      "match_at0":float((mb==0).sum())/n,
      "match_le2":float((mb<=2).sum())/n,
      "affclip_at0":float((affc<=0).sum())/n,
      "affclip_le2":float((affc<=2).sum())/n,
      "orig_at_top":float((b>=HI).sum())/n,
      "match_at_top":float((mb>=HI).sum())/n,
      "affclip_at_top":float((affc>=HI).sum())/n,
      # entropy-ish: how many distinct values survive in the bottom decile of mass
      "match_distinct":int(np.unique(mb).size), "orig_distinct":int(np.unique(b).size),
      "affq_distinct":int(np.unique(np.rint(affc)).size)})
nl=[]
for s in sorted(json.loads((D/"nonloc60.json").read_text())):
    b,_=blk(s); n=b.size
    nl.append({"sample":s,"at0":float((b==0).sum())/n,"le2":float((b<=2).sum())/n,
               "at_top":float((b>=HI).sum())/n,"mean":float(b.mean()),"std":float(b.std())})
res={"located60":out,"nonlocated60_orig":nl,
 "medians_located":{k:round(st.median([r[k] for r in out]),6) for k in out[0] if k!="sample"},
 "medians_nonlocated_orig":{k:round(st.median([r[k] for r in nl]),6) for k in nl[0] if k!="sample"}}
Path(D/"zerocheck.json").write_text(json.dumps(res,indent=1))
print(json.dumps(res["medians_located"],indent=1)); print(json.dumps(res["medians_nonlocated_orig"],indent=1))
print("match_at0 range", round(min(r["match_at0"] for r in out),5), round(max(r["match_at0"] for r in out),5))
print("affclip_at0 range", round(min(r["affclip_at0"] for r in out),5), round(max(r["affclip_at0"] for r in out),5))
print("nonloc at0 range", round(min(r["at0"] for r in nl),5), round(max(r["at0"] for r in nl),5))
