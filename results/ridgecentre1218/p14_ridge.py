"""Measure ridge-minus-centre on the published tight-contact crops.

Follows PREREGISTRATION.md (frozen ad186ffbf127216b). The normal is oriented
A to B and then flipped so positive is outward from the scroll axis, so a
negative offset means the CT ridge sits inward of the label-run centre.
"""
import glob, json, os
import numpy as np

SRC = "/mnt/vesuvius/kaggle_tightgap1218/crops"
OUT = "/mnt/vesuvius/ridgecentre1218"
os.makedirs(OUT, exist_ok=True)
# scroll axis in level-0 voxels: volume is 23247 x 7593 x 7593, axis at y/x centre
AXIS = np.array([7593 / 2.0, 7593 / 2.0])
STEP = 0.25


def profile(vol, p, n, half):
    ts = np.arange(-half, half + 1e-9, STEP)
    out = []
    for t in ts:
        q = p + t * n
        i = np.floor(q).astype(int)
        f = q - i
        if np.any(i < 0) or np.any(i + 1 >= vol.shape):
            return None, None
        c = vol[i[0]:i[0] + 2, i[1]:i[1] + 2, i[2]:i[2] + 2].astype(np.float32)
        w = np.array([[[(1 - f[0]) * (1 - f[1]) * (1 - f[2]), (1 - f[0]) * (1 - f[1]) * f[2]],
                       [(1 - f[0]) * f[1] * (1 - f[2]), (1 - f[0]) * f[1] * f[2]]],
                      [[f[0] * (1 - f[1]) * (1 - f[2]), f[0] * (1 - f[1]) * f[2]],
                       [f[0] * f[1] * (1 - f[2]), f[0] * f[1] * f[2]]]])
        out.append(float((c * w).sum()))
    return ts, np.array(out)


def run_centre(lab, p, n, ids, half=20):
    """midpoint of the labelled run through p along n"""
    lo = hi = 0.0
    for sgn, store in ((1, "hi"), (-1, "lo")):
        t = 0.0
        while t < half:
            t += 0.5
            q = np.round(p + sgn * t * n).astype(int)
            if np.any(q < 0) or np.any(q >= lab.shape):
                break
            if int(lab[q[0], q[1], q[2]]) not in ids:
                break
        if store == "hi":
            hi = t - 0.5
        else:
            lo = t - 0.5
    if hi + lo < 1.0:
        return None
    return (hi - lo) / 2.0   # offset of run centre from p, along n


rows = []
disc = {"corridor": 0, "no_run": 0, "flat": 0, "no_ids": 0}
for f in sorted(glob.glob(SRC + "/*.npz")):
    d = np.load(f, allow_pickle=True)
    A, B = int(d["A_id"]), int(d["B_id"])
    lab, vol = d["instance"], d["intensity"]
    ids = set(np.unique(lab).tolist())
    if A not in ids or B not in ids:
        disc["no_ids"] += 1
        continue
    ca = np.argwhere(lab == A).mean(0)
    cb = np.argwhere(lab == B).mean(0)
    n = cb - ca
    nn = np.linalg.norm(n)
    if nn < 1e-6:
        disc["no_run"] += 1
        continue
    n = n / nn
    # orient outward from the scroll axis, using the crop's global position
    site = np.array([float(v) for v in d["site"]])
    radial = site[1:] - AXIS
    if np.dot(n[1:], radial) < 0:
        n = -n
    p0 = np.array([64.0, 64.0, 64.0])
    c = run_centre(lab, p0, n, {A, B})
    if c is None:
        disc["no_run"] += 1
        continue
    centre = p0 + c * n
    rec = {"file": os.path.basename(f), "band": str(d["band"]), "gap": float(d["gap"])}
    okany = False
    for half in (3, 4, 8):
        ts, prof = profile(vol, centre, n, half)
        if ts is None:
            continue
        if float(prof.max() - prof.min()) < 1.0:
            disc["flat"] += 1
            continue
        rec["off_%d" % half] = float(ts[int(np.argmax(prof))])
        okany = True
    if not okany:
        disc["corridor"] += 1
        continue
    rows.append(rec)

def summarise(key):
    v = np.array([r[key] for r in rows if key in r])
    if not len(v):
        return None
    rng = np.random.default_rng(0)
    boot = [float(np.median(rng.choice(v, len(v)))) for _ in range(2000)]
    return dict(n=len(v), median=round(float(np.median(v)), 4),
                median_abs=round(float(np.median(np.abs(v))), 4),
                q10=round(float(np.percentile(v, 10)), 4),
                q90=round(float(np.percentile(v, 90)), 4),
                frac_negative=round(float((v < 0).mean()), 4),
                ci95=[round(float(np.percentile(boot, 2.5)), 4),
                      round(float(np.percentile(boot, 97.5)), 4)])

out = dict(prereg="ad186ffbf127216b", crops_used=len(rows), discarded=disc,
           corridors={k: summarise("off_%d" % k) for k in (3, 4, 8)})
json.dump({"summary": out, "rows": rows}, open(OUT + "/ridge_offsets.json", "w"), indent=1)
print(json.dumps(out, indent=1))
