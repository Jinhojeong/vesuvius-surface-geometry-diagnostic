"""Why do 44 contact crops carry neither split id? Measure, do not assume."""
import glob, json
import numpy as np

OUT = "/mnt/vesuvius/tightgap1218"
res = {"neither": [], "both": []}
for f in sorted(glob.glob(OUT + "/crops/*.npz")):
    d = np.load(f, allow_pickle=True)
    if "A_id" not in d:
        continue
    A, B = int(d["A_id"]), int(d["B_id"])
    lab = d["instance"]
    ids = set(np.unique(lab).tolist())
    bases = json.loads(str(d["id_bases"]))
    site = [int(v) for v in d["site"]]
    centre = lab[60:68, 60:68, 60:68]
    entry = dict(band=str(d["band"]), gap=float(d["gap"]),
                 nblocks=len(bases), n_ids=len(ids - {0}),
                 centre_labelled=float((centre > 0).mean()),
                 A_local=int(d["A"]), B_local=int(d["B"]))
    key = "both" if (A in ids and B in ids) else ("neither" if (A not in ids and B not in ids) else None)
    if key:
        res[key].append(entry)

for k in ("both", "neither"):
    v = res[k]
    if not v:
        continue
    print("== %s: n=%d" % (k, len(v)))
    print("   median blocks per crop %.1f | median distinct ids %.0f | median centre-labelled %.3f"
          % (np.median([e["nblocks"] for e in v]), np.median([e["n_ids"] for e in v]),
             np.median([e["centre_labelled"] for e in v])))
    import collections
    print("   by band:", dict(collections.Counter(e["band"] for e in v)))
