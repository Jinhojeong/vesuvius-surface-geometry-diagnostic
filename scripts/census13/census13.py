"""Diego's 1.3x question: stratified low-threshold census pass.

Same detector as census8k (p1218_repair_c.detect_fused) with THICK_F dropped to
1.3, on a deterministic 100-tile stratified sample (every 13th tile by sorted
order, which walks the whole z range). Retains the full per-site ratio value so
the 1.3-1.6 band Diego asked about is separable from the >=1.6 population the
main census already enumerated. CPU only, resumable, one JSON per tile.
"""
import glob, json, os, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
sys.path.insert(0, "/mnt/vesuvius")
import p1218_repair_c as R

R.MAX_SOLVES = 10 ** 9
R.THICK_F = 1.3
SRC = "/mnt/vesuvius/p1218_repair_v3/blocks_repaired"
OUT = "/mnt/vesuvius/census13"
NPTS = 8000


def tile_job(path):
    slab, name = path.split("/")[-2], path.split("/")[-1][:-4]
    od = os.path.join(OUT, slab); os.makedirs(od, exist_ok=True)
    op = os.path.join(od, name + ".json")
    if os.path.exists(op):
        return "skip"
    t0 = time.time()
    with np.load(path) as d:
        lab = d["labels"].astype(np.int32)
    mask = lab > 0
    if mask.sum() < 200000:
        json.dump({"tile": f"{slab}/{name}", "thin": True, "clusters": []}, open(op, "w"))
        return "thin"
    nrm = R.normals_of(mask)
    R.NPTS = NPTS
    sites, med = R.detect_fused(lab, mask, nrm, np.random.default_rng(7))
    clusters = R.cluster_sites(sites)
    rec = {"tile": f"{slab}/{name}", "npts": NPTS, "thick_f": R.THICK_F,
           "median_run": float(med), "n_raw_sites": int(len(sites)),
           "n_clusters": int(len(clusters)),
           "clusters": [{"z": int(r[0][0]), "y": int(r[0][1]), "x": int(r[0][2]),
                         "inst": int(r[1]), "th": round(float(r[2]), 2),
                         "ratio": round(float(r[3]), 3), "n_sites": len(g)}
                        for r, g in clusters],
           "seconds": round(time.time() - t0, 1)}
    json.dump(rec, open(op, "w"))
    return "ok"


def main():
    tiles = sorted(glob.glob(f"{SRC}/*/*.npz"))[::13][:100]
    print(f"{len(tiles)} sampled tiles, THICK_F={R.THICK_F}, NPTS={NPTS}", flush=True)
    t0 = time.time(); done = 0
    with ProcessPoolExecutor(max_workers=8) as ex:
        for _ in ex.map(tile_job, tiles):
            done += 1
            if done % 10 == 0:
                print(f"[{done}] {(time.time()-t0)/60:.0f}min", flush=True)
    print("DONE", flush=True)


main()
