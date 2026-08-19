"""P11 label fix: instance ids in the repair blocks are block-local, so a crop
that spans blocks needs a disambiguated id space.

Each contributing block gets a unique base offset inside the crop, and the two
ids the split separated are remapped through the base of the block the site
belongs to. Without this, an id from a neighbouring block can collide with A
or B and make a contact look present when it is not. The first pass had that
defect; its both_present count is superseded by the one written here.
"""
import glob, json
import numpy as np

OUT = "/mnt/vesuvius/tightgap1218"
BLOCKS = "/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired"
N = 128
BASE = 1000000

idx = []
for f in sorted(glob.glob(BLOCKS + "/*/*.npz")):
    z = np.load(f)
    idx.append((int(z["z0"]), int(z["y0"]), int(z["x0"]), z["labels"].shape, f))
print("blocks indexed:", len(idx), flush=True)


def stitch(z0, y0, x0, site):
    out = np.zeros((N, N, N), np.int32)
    prov = []
    site_base = None
    rank = 0
    for bz, by, bx, sh, f in idx:
        if bz >= z0 + N or bz + sh[0] <= z0:
            continue
        if by >= y0 + N or by + sh[1] <= y0:
            continue
        if bx >= x0 + N or bx + sh[2] <= x0:
            continue
        z = np.load(f)
        lab = z["labels"]
        az0, az1 = max(z0, bz), min(z0 + N, bz + sh[0])
        ay0, ay1 = max(y0, by), min(y0 + N, by + sh[1])
        ax0, ax1 = max(x0, bx), min(x0 + N, bx + sh[2])
        sub = lab[az0 - bz:az1 - bz, ay0 - by:ay1 - by, ax0 - bx:ax1 - bx].astype(np.int32)
        rank += 1
        base = BASE * rank
        m = sub > 0
        out[az0 - z0:az1 - z0, ay0 - y0:ay1 - y0, ax0 - x0:ax1 - x0][m] = sub[m] + base
        prov.append({"block": "/".join(f.split("/")[-2:]), "base": base})
        sz, sy, sx = site
        if bz <= sz < bz + sh[0] and by <= sy < by + sh[1] and bx <= sx < bx + sh[2]:
            site_base = base
    return out, prov, site_base


stats = {"contact": {"n": 0, "both": 0, "one": 0, "neither": 0,
                     "no_site_block": 0, "multi_block": 0},
         "control": {"n": 0, "multi_block": 0}}
lf = {"contact": [], "control": []}
for arm in ("crops", "control"):
    key = "contact" if arm == "crops" else "control"
    for f in sorted(glob.glob("%s/%s/*.npz" % (OUT, arm))):
        d = dict(np.load(f, allow_pickle=True))
        site = [int(v) for v in d["site"]]
        z0, y0, x0 = [v - N // 2 for v in site]
        lab, prov, sb = stitch(z0, y0, x0, site)
        d["instance"] = lab
        d["surface"] = (lab > 0).astype(np.uint8)
        d["id_bases"] = json.dumps(prov)
        stats[key]["n"] += 1
        lf[key].append(float((lab > 0).mean()))
        if len(prov) > 1:
            stats[key]["multi_block"] += 1
        if key == "contact":
            if sb is None:
                stats["contact"]["no_site_block"] += 1
                d["A_id"] = -1
                d["B_id"] = -1
            else:
                A = sb + int(d["A"])
                B = sb + int(d["B"])
                d["A_id"] = A
                d["B_id"] = B
                ids = set(np.unique(lab).tolist())
                if A in ids and B in ids:
                    stats["contact"]["both"] += 1
                elif A in ids or B in ids:
                    stats["contact"]["one"] += 1
                else:
                    stats["contact"]["neither"] += 1
        np.savez_compressed(f, **d)
    print("done", arm, flush=True)

for k in lf:
    stats[k]["labelled_frac_median"] = round(float(np.median(lf[k])), 4)
stats["note"] = ("ids are block-local in the source, so each contributing block is "
                 "offset by a unique base inside the crop; A_id and B_id are the "
                 "split pair remapped through the site block base")
json.dump(stats, open(OUT + "/labels_summary.json", "w"), indent=1)
print(json.dumps(stats, indent=1))
