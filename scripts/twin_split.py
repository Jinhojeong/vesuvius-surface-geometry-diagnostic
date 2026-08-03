"""Splitter identity accuracy on the twin, against true turn labels.

Simulates fusion by merging label pairs (the weld pair plus two natural
contact pairs), then runs the production repair path (p1218_repair_c.
repair_cluster: neighbour anchors -> ADL-RW solve -> assignment) at sites on
the true interface, and scores every assigned voxel against the true turn id.
Correct assignment: the voxel of the inner true sheet goes to the inner
neighbour id, outer to outer. First ground-truthed identity accuracy for the
continuity splitter.
"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, "/mnt/vesuvius")
import p1218_repair_c as R

# The twin is contact-rich and z-uniform, so a 64^3 crop's union component
# exceeds the production near-solid guard (110k) everywhere. Same solve at
# 40^3 keeps the union around 34k; only the crop extent changes.
R.HALF = 20

RUN = "/mnt/vesuvius/twin/run1"
PAIRS = [(22, 23), (11, 12), (28, 29)]   # (inner_id, outer_id) merged
N_PER_PAIR = 150

true_lab = np.load(f"{RUN}/labels.npy")
fused = true_lab.copy()
for a, b in PAIRS:
    fused[true_lab == a] = b

nz = true_lab.shape[0]
zlo, zhi = R.HALF + 2, nz - R.HALF - 2


def sites_for(a, b):
    """Interface voxels of the merged pair, central z band, subsampled."""
    am = true_lab == a
    bm = true_lab == b
    touch = am & (
        np.roll(bm, 1, 1) | np.roll(bm, -1, 1)
        | np.roll(bm, 1, 2) | np.roll(bm, -1, 2))
    touch[:zlo] = False
    touch[zhi:] = False
    pts = np.argwhere(touch)
    if len(pts) == 0:
        return []
    stride = max(len(pts) // N_PER_PAIR, 1)
    return [tuple(int(v) for v in p) for p in pts[::stride][:N_PER_PAIR]]


def job(args):
    a, b, p0 = args
    rep = (p0, b, 0.0, 2.0)
    dec, rec, assigns = R.repair_cluster(fused, rep, [rep])
    if dec != "SPLIT":
        return (a, b, dec, 0, 0, None)
    inner_nb, outer_nb = a - 1, b + 1
    ok = bad = other = 0
    for (gz, gy, gx), nid in assigns:
        t = int(true_lab[gz, gy, gx])
        if t not in (a, b):
            continue
        want = inner_nb if t == a else outer_nb
        if nid == want:
            ok += 1
        elif nid in (inner_nb, outer_nb):
            bad += 1
        else:
            other += 1
    return (a, b, "SPLIT", ok, bad, rec["conf"] if rec else None)


def main():
    jobs = []
    for a, b in PAIRS:
        ss = sites_for(a, b)
        print(f"pair ({a},{b}): {len(ss)} interface sites", flush=True)
        jobs += [(a, b, p) for p in ss]
    res = []
    with ProcessPoolExecutor(max_workers=4) as ex:
        for k, r in enumerate(ex.map(job, jobs)):
            res.append(r)
            if (k + 1) % 50 == 0:
                print(f"[{k+1}/{len(jobs)}]", flush=True)
    from collections import Counter
    dec = Counter(r[2] for r in res)
    print("decisions:", dict(dec))
    tot_ok = sum(r[3] for r in res)
    tot_bad = sum(r[4] for r in res)
    n_split = sum(1 for r in res if r[2] == "SPLIT")
    acc = 100.0 * tot_ok / max(tot_ok + tot_bad, 1)
    print(f"SPLIT at {n_split}/{len(res)} sites "
          f"({100*n_split/len(res):.1f}%)")
    print(f"IDENTITY ACCURACY over assigned voxels: {acc:.1f}% "
          f"({tot_ok:,} ok / {tot_bad:,} wrong-side)")
    per_pair = {}
    for a, b in PAIRS:
        rr = [r for r in res if r[0] == a and r[2] == "SPLIT"]
        o = sum(r[3] for r in rr)
        w = sum(r[4] for r in rr)
        per_pair[f"{a}-{b}"] = {
            "n_split": len(rr),
            "acc_pct": round(100.0 * o / max(o + w, 1), 2)}
        print(f"  pair {a}-{b}: split {len(rr)}, "
              f"acc {100.0*o/max(o+w,1):.1f}%")
    json.dump({"decisions": dict(dec), "n_sites": len(res),
               "n_split": n_split, "identity_acc_pct": round(acc, 2),
               "ok": tot_ok, "wrong": tot_bad, "per_pair": per_pair},
              open(f"{RUN}/split_verdict.json", "w"), indent=1)


main()
