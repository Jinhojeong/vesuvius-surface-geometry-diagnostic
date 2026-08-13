#!/usr/bin/env python3
"""Full PHerc0332 topological repair run (v2.0 machinery, frozen constants).

Scope and design, as approved:
  - all 109,280 census clusters (sites_0332.csv.gz); frozen 1218 repair
    constants (HALF=32, TAU=0.10, APPLY_R=14, stock anchor gates);
    group=[rep] adaptation as piloted (census stores representatives only).
  - clusters in tiles not owned by any of the 137 grid crops: recorded
    UNOWNED, not solved.
  - pass 1: solve every owned cluster in its owning crop (crop-local ids).
  - pass 2 (OOB recovery): every cluster whose site lacks the full HALF=32
    margin in its owning crop is re-solved in the crop that contains it with
    full margin, chosen deterministically (nearest crop centre by Euclidean
    distance; ties by lexicographic centre). No candidate -> OOB_UNRECOVERABLE.
  - repairs applied to a COPY of each crop array, pass-1 assigns first then
    pass-2, each ratio-descending. Originals untouched. Invariants asserted
    per crop: identical foreground mask, no new instance id.
  - retile to 256^3 blocks (same layout/min-vox rule as kaggle_p0332_sep)
    under /mnt/vesuvius/p0332_repair_v1/blocks_repaired.
  - validation per crop while in memory: ray recast on every SPLIT site
    (both passes, reported separately) + double-thickness column share
    before/after (census RESCALE ray constants, NPTS=12000, deterministic
    per-crop seed = crop centre). column/ids_in_run/fused_share follow
    validate_c.py with the float32 mask hoisted out of column().
  - resource discipline: site-solve workers = 8 while labelcov1218's
    coverage_run.py is alive, else 12; checked before each crop.
  - resumable per crop (crops/<tag>.json is written last); a structurally
    failed crop is recorded and skipped; more than 5 failed crops aborts.

Usage:
  python3 full_run.py            # full run (resumable)
  python3 full_run.py --only TAG # single crop, e.g. --only 1152-1600-1152
"""
import csv
import gzip
import hashlib
import json
import os
import subprocess
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy import ndimage as ndi

sys.path.insert(0, "/mnt/vesuvius")
import p1218_repair_c as R
from unmerge_v0 import normals_of

ROOT = "/mnt/vesuvius/p5_pilot0332/full"
CROPS_DIR = os.path.join(ROOT, "crops")
BLOCKS = "/mnt/vesuvius/p0332_repair_v1/blocks_repaired"
GRID = "/mnt/vesuvius/pilot0332/crop_grid.json"
CENSUS = "/mnt/vesuvius/kaggle_p0332_sep/census/sites_0332.csv.gz"
NPY_DIR = "/mnt/vesuvius/pilot0332/vesuvius-sheet-tools/output"

HALF, SIZE = 32, 512
SPAN, STEP = 7.5, 0.3125
K = int(2 * SPAN / STEP) + 1
THICK_F = 1.6
NPTS = 12000
MIN_VOX = 20000
MAX_FAILED = 5

# ---------------------------------------------------------------- assignment

def tag_of(cen):
    return f"{cen[0]}-{cen[1]}-{cen[2]}"


def npy_of(cen):
    return os.path.join(NPY_DIR,
                        f"sheets_scroll3_L1_c{cen[0]}-{cen[1]}-{cen[2]}.npy")


def build_assignments():
    grid = json.load(open(GRID))
    centres = [tuple(c["centre_l1"]) for c in grid["crops"]]
    centre_set = set(centres)

    def owner_of_corner(corner):
        cands = [()]
        for c in corner:
            cands = [t + (cc,) for t in cands for cc in (c + 256, c)]
        hits = [t for t in cands if t in centre_set]
        return hits[0] if len(hits) == 1 else None

    assign = {tag_of(c): {"centre": list(c),
                          "origin": [v - 256 for v in c],
                          "pass1": [], "pass2": []} for c in centres}
    unowned = oob_unrec = oob_total = 0
    for r in csv.DictReader(gzip.open(CENSUS, "rt")):
        slab, tile = r["tile"].split("/")
        corner = (int(slab[1:]), int(tile[1:6]), int(tile[8:13]))
        own = owner_of_corner(corner)
        if own is None:
            unowned += 1
            continue
        g = (int(r["gz"]), int(r["gy"]), int(r["gx"]))
        cl = {"inst": int(r["inst"]), "th": float(r["thickness_vox"]),
              "ratio": float(r["ratio"]), "n_sites": int(r["cluster_sites"]),
              "g": list(g)}
        o = [v - 256 for v in own]
        p = [g[i] - o[i] for i in range(3)]
        assign[tag_of(own)]["pass1"].append(
            dict(cl, z=p[0], y=p[1], x=p[2]))
        if any(p[i] < HALF or p[i] >= SIZE - HALF for i in range(3)):
            oob_total += 1
            cands = []
            for cen in centres:
                q = [g[i] - (cen[i] - 256) for i in range(3)]
                if all(HALF <= q[i] < SIZE - HALF for i in range(3)):
                    cands.append(cen)
            if not cands:
                oob_unrec += 1
                continue
            best = min(cands, key=lambda c: (
                sum((g[i] - c[i]) ** 2 for i in range(3)), c))
            ob = [v - 256 for v in best]
            assign[tag_of(best)]["pass2"].append(
                dict(cl, z=g[0] - ob[0], y=g[1] - ob[1], x=g[2] - ob[2]))
    for e in assign.values():
        for key in ("pass1", "pass2"):
            e[key].sort(key=lambda c: (-c["ratio"], c["z"], c["y"], c["x"]))
    meta = {"unowned": unowned, "oob_geometric": oob_total,
            "oob_unrecoverable": oob_unrec,
            "n_pass1": sum(len(e["pass1"]) for e in assign.values()),
            "n_pass2": sum(len(e["pass2"]) for e in assign.values())}
    json.dump({"meta": meta, "crops": assign},
              open(os.path.join(ROOT, "assignments.json"), "w"))
    print("assignments:", json.dumps(meta), flush=True)
    return meta, assign


# ---------------------------------------------------------------- validator

def column(maskf, nrm, p0):
    offs = np.linspace(-SPAN, SPAN, K)
    nv = nrm[:, p0[0], p0[1], p0[2]]
    cc = (np.array(p0, float)[None, :] + offs[:, None] * nv[None, :]).T
    on = ndi.map_coordinates(maskf, cc, order=0, mode="constant") > 0.5
    c = K // 2
    if not on[c]:
        return None
    a = c
    while a > 0 and on[a - 1]:
        a -= 1
    b = c
    while b < K - 1 and on[b + 1]:
        b += 1
    return a, b, nv


def ids_in_run(lab, maskf, nrm, p0):
    col = column(maskf, nrm, p0)
    if col is None:
        return None
    a, b, nv = col
    ids = set()
    for kk in range(a, b + 1):
        q = np.round(np.array(p0, float) + (kk * STEP - SPAN) * nv).astype(int)
        if ((q >= 0) & (q < np.array(lab.shape))).all():
            v = int(lab[q[0], q[1], q[2]])
            if v > 0:
                ids.add(v)
    return ids


def fused_share(lab, mask, maskf, nrm, rng):
    pts = np.argwhere(mask)
    pts = pts[rng.permutation(len(pts))[:NPTS * 3]]
    m_ = int(SPAN) + 2
    pts = pts[((pts >= m_) & (pts < np.array(mask.shape) - m_)).all(1)][:NPTS]
    rows = []
    for p0 in pts:
        col = column(maskf, nrm, p0)
        if col is None:
            continue
        a, b, nv = col
        rows.append((p0, nv, a, b, (b - a + 1) * STEP))
    if len(rows) < 200:
        return None
    med = float(np.median([r[4] for r in rows]))
    thick = fused = 0
    for p0, nv, a, b, th in rows:
        if th <= THICK_F * med:
            continue
        ids = set()
        for kk in range(a, b + 1):
            q = np.round(np.array(p0, float) + (kk * STEP - SPAN) * nv).astype(int)
            if ((q >= 0) & (q < np.array(mask.shape))).all():
                v = int(lab[q[0], q[1], q[2]])
                if v > 0:
                    ids.add(v)
        if not ids:
            continue
        thick += 1
        fused += (len(ids) == 1)
    return thick, fused


# ---------------------------------------------------------------- solve

LAB = None


def solve_one(cl):
    t0 = time.time()
    p0 = (cl["z"], cl["y"], cl["x"])
    if not all(0 <= p0[i] < LAB.shape[i] for i in range(3)):
        return "OUT_OF_BOUNDS", None, [], time.time() - t0
    lab_id = int(LAB[p0])
    if lab_id <= 0:
        return "SITE_NOT_IN_MASK", None, [], time.time() - t0
    rep = (p0, lab_id, cl["th"], cl["ratio"])
    verdict, rec, assigns = R.repair_cluster(LAB, rep, [rep])
    return verdict, rec, assigns, time.time() - t0


def labelcov_running():
    return subprocess.run(["pgrep", "-f", "coverage_run.py"],
                          stdout=subprocess.DEVNULL).returncode == 0


def run_crop(tag, entry):
    global LAB
    t0 = time.time()
    workers = 8 if labelcov_running() else 12
    LAB = np.load(npy_of(entry["centre"])).astype(np.int32)
    if LAB.shape != (SIZE, SIZE, SIZE):
        raise RuntimeError(f"bad shape {LAB.shape}")
    out = LAB.copy()
    dec = {1: {}, 2: {}}
    recs = {1: [], 2: []}
    solve_s = 0.0
    for pno, key in ((1, "pass1"), (2, "pass2")):
        todo = entry[key]
        if not todo:
            continue
        with ProcessPoolExecutor(max_workers=workers) as ex:
            results = list(ex.map(solve_one, todo, chunksize=4))
        for cl, (v, rec, assigns, dt) in zip(todo, results):
            dec[pno][v] = dec[pno].get(v, 0) + 1
            solve_s += dt
            if rec is not None:
                rec["pass"] = pno
                rec["g"] = cl["g"]
                recs[pno].append(rec)
                for (gz_, gy_, gx_), nid in assigns:
                    out[gz_, gy_, gx_] = nid
    # invariants
    mask_ok = bool(((out > 0) == (LAB > 0)).all())
    new_ids = np.setdiff1d(np.unique(out), np.unique(LAB))
    applied = int((out != LAB).sum())
    if not mask_ok or new_ids.size:
        raise RuntimeError(f"invariant broken: mask_ok={mask_ok} "
                           f"new_ids={new_ids.tolist()[:5]}")
    # validation
    mask = LAB > 0
    maskf = mask.astype(np.float32)
    nrm = normals_of(mask)
    val = {}
    for pno in (1, 2):
        checked = fixed = 0
        for rec in recs[pno]:
            after = ids_in_run(out, maskf, nrm, tuple(rec["site"]))
            if after is None:
                continue
            checked += 1
            fixed += (len(after) >= 2)
        val[pno] = {"checked": checked, "fixed": fixed}
    seed = list(entry["centre"])
    b = fused_share(LAB, mask, maskf, nrm, np.random.default_rng(seed))
    a = (b if applied == 0 else
         fused_share(out, mask, maskf, nrm, np.random.default_rng(seed)))
    # retile
    o = entry["origin"]
    tiles_written = 0
    for dz in (0, 256):
        for dy in (0, 256):
            for dx in (0, 256):
                sub = out[dz:dz + 256, dy:dy + 256, dx:dx + 256]
                if int((sub > 0).sum()) < MIN_VOX:
                    continue
                slab = f"z{o[0] + dz:05d}"
                name = f"y{o[1] + dy:05d}_x{o[2] + dx:05d}"
                td = os.path.join(BLOCKS, slab)
                os.makedirs(td, exist_ok=True)
                np.savez_compressed(os.path.join(td, name + ".npz"),
                                    labels=sub, z0=o[0] + dz, y0=o[1] + dy,
                                    x0=o[2] + dx)
                tiles_written += 1
    res = {"tag": tag, "centre": entry["centre"], "origin": o,
           "n_pass1": len(entry["pass1"]), "n_pass2": len(entry["pass2"]),
           "decisions_pass1": dec[1], "decisions_pass2": dec[2],
           "validation": {str(k): v for k, v in val.items()},
           "thick_before": b[0] if b else None,
           "fused_before": b[1] if b else None,
           "thick_after": a[0] if a else None,
           "fused_after": a[1] if a else None,
           "applied_voxels": applied, "mask_invariant": mask_ok,
           "no_new_ids": bool(new_ids.size == 0),
           "tiles_written": tiles_written, "workers": workers,
           "solve_core_s": round(solve_s, 1),
           "wall_s": round(time.time() - t0, 1),
           "repairs": recs[1] + recs[2]}
    json.dump(res, open(os.path.join(CROPS_DIR, tag + ".json"), "w"))
    LAB = None
    return res


# ---------------------------------------------------------------- provenance

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def build_manifest(assign):
    mp = os.path.join(ROOT, "input_manifest.json")
    if os.path.exists(mp):
        return json.load(open(mp))
    t0 = time.time()
    man = {"census_csv": {"path": CENSUS, "sha256": sha256_file(CENSUS)},
           "crop_npys": {}}
    for tag, e in sorted(assign.items()):
        p = npy_of(e["centre"])
        man["crop_npys"][tag] = {"path": p, "sha256": sha256_file(p)}
    json.dump(man, open(mp, "w"), indent=1)
    print(f"manifest: hashed {len(man['crop_npys'])} npys + census "
          f"({(time.time()-t0)/60:.1f}m)", flush=True)
    return man


# ---------------------------------------------------------------- aggregate

def aggregate(meta, assign, manifest, failed, t_start):
    hist = {"UNOWNED": meta["unowned"],
            "OOB_UNRECOVERABLE": meta["oob_unrecoverable"]}
    val = {"1": {"checked": 0, "fixed": 0}, "2": {"checked": 0, "fixed": 0}}
    tb = fb = ta = fa = 0
    per_crop = {}
    applied_total = 0
    shares = []
    for tag in sorted(assign):
        p = os.path.join(CROPS_DIR, tag + ".json")
        if not os.path.exists(p):
            continue
        r = json.load(open(p))
        for k, v in r["decisions_pass1"].items():
            hist[k] = hist.get(k, 0) + v
        for k, v in r["decisions_pass2"].items():
            kk = k + "_OOB_RECOVERED"
            hist[kk] = hist.get(kk, 0) + v
        for pno in ("1", "2"):
            val[pno]["checked"] += r["validation"][pno]["checked"]
            val[pno]["fixed"] += r["validation"][pno]["fixed"]
        if r["thick_before"]:
            tb += r["thick_before"]
            fb += r["fused_before"]
            ta += r["thick_after"]
            fa += r["fused_after"]
            shares.append({"tag": tag,
                           "before": round(r["fused_before"] / r["thick_before"], 4),
                           "after": round(r["fused_after"] / r["thick_after"], 4)})
        applied_total += r["applied_voxels"]
        per_crop[tag] = {"decisions_pass1": r["decisions_pass1"],
                         "decisions_pass2": r["decisions_pass2"],
                         "applied_voxels": r["applied_voxels"],
                         "mask_invariant": r["mask_invariant"],
                         "no_new_ids": r["no_new_ids"],
                         "wall_s": r["wall_s"]}
    deltas = sorted(s["after"] - s["before"] for s in shares)
    q = (lambda f: round(deltas[int(f * (len(deltas) - 1))], 4)) if deltas \
        else (lambda f: None)
    checked_all = val["1"]["checked"] + val["2"]["checked"]
    fixed_all = val["1"]["fixed"] + val["2"]["fixed"]
    results = {
        "verdict_histogram": dict(sorted(hist.items())),
        "validation": {
            "pass1": dict(val["1"], rate=round(
                val["1"]["fixed"] / max(val["1"]["checked"], 1), 4)),
            "pass2_oob_recovery": dict(val["2"], rate=round(
                val["2"]["fixed"] / max(val["2"]["checked"], 1), 4)),
            "overall": {"checked": checked_all, "fixed": fixed_all,
                        "rate": round(fixed_all / max(checked_all, 1), 4)}},
        "thick_column_share": {
            "pooled_before": round(fb / max(tb, 1), 4),
            "pooled_after": round(fa / max(ta, 1), 4),
            "thick_cols_sampled": tb,
            "crops_measured": len(shares),
            "delta_min": q(0.0), "delta_p25": q(0.25), "delta_median": q(0.5),
            "delta_p75": q(0.75), "delta_max": q(1.0),
            "per_crop": shares},
        "applied_voxels_total": applied_total,
        "invariants_all_ok": all(c["mask_invariant"] and c["no_new_ids"]
                                 for c in per_crop.values()),
        "crops_done": len(per_crop), "crops_failed": failed,
        "runtime_h": round((time.time() - t_start) / 3600.0, 2),
    }
    json.dump(results, open(os.path.join(ROOT, "results_full.json"), "w"),
              indent=1)
    prov = {
        "run": "p0332 full repair v1 (frozen v2.0 constants)",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "inputs": manifest,
        "constants": {
            "repair": {"HALF": R.HALF, "TAU": R.TAU, "APPLY_R": R.APPLY_R,
                       "CLUSTER_R": R.CLUSTER_R,
                       "anchor_gates": "count>=30, |mean_off|>0.5, "
                                       "inplane<=18, |s_off|<=10 (stock)"},
            "validation_rays": {"SPAN": SPAN, "STEP": STEP, "K": K,
                                "THICK_F": THICK_F, "NPTS": NPTS,
                                "seed": "per-crop centre [cz, cy, cx]"}},
        "pass_definitions": {
            "pass1": "each owned census cluster solved in its owning crop "
                     "(tile corner -> unique crop on the 448-step lattice); "
                     "group=[rep]",
            "pass2_oob_recovery": "sites lacking the HALF=32 margin in the "
                                  "owning crop re-solved in the nearest-"
                                  "centre crop with full margin (ties: "
                                  "lexicographic centre)",
            "UNOWNED": "census rows in tiles owned by no grid crop; not "
                       "solved",
            "OOB_UNRECOVERABLE": "no crop contains the site with full "
                                 "margin"},
        "assignment_meta": meta,
        "per_crop": per_crop,
        "outputs": {"blocks": BLOCKS, "crop_reports": CROPS_DIR,
                    "results": os.path.join(ROOT, "results_full.json")},
    }
    pp = os.path.join(ROOT, "provenance.json")
    json.dump(prov, open(pp, "w"), indent=1)
    print("PROVENANCE", pp, sha256_file(pp), flush=True)
    print(json.dumps({k: v for k, v in results.items()
                      if k != "thick_column_share"}, indent=1), flush=True)
    print("pooled share:", results["thick_column_share"]["pooled_before"],
          "->", results["thick_column_share"]["pooled_after"], flush=True)
    print("FULL RUN COMPLETE", flush=True)


# ---------------------------------------------------------------- main

def main():
    only = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv \
        else None
    os.makedirs(CROPS_DIR, exist_ok=True)
    os.makedirs(BLOCKS, exist_ok=True)
    t_start = time.time()
    ap = os.path.join(ROOT, "assignments.json")
    if os.path.exists(ap):
        d = json.load(open(ap))
        meta, assign = d["meta"], d["crops"]
    else:
        meta, assign = build_assignments()
    manifest = build_manifest(assign)
    tags = [only] if only else sorted(assign)
    failed = []
    for i, tag in enumerate(tags):
        if os.path.exists(os.path.join(CROPS_DIR, tag + ".json")):
            continue
        try:
            r = run_crop(tag, assign[tag])
            print(f"[{i+1}/{len(tags)}] {tag}: p1 {r['decisions_pass1']} "
                  f"p2 {r['decisions_pass2']} applied {r['applied_voxels']} "
                  f"({r['wall_s']:.0f}s, w{r['workers']})", flush=True)
        except Exception:
            failed.append(tag)
            print(f"[{i+1}/{len(tags)}] {tag}: FAILED\n"
                  f"{traceback.format_exc()[-600:]}", flush=True)
            if len(failed) > MAX_FAILED:
                print(f"ABORT: {len(failed)} crops failed: {failed}",
                      flush=True)
                sys.exit(2)
    if only:
        print("single-crop run done (no aggregate)", flush=True)
        return
    aggregate(meta, assign, manifest, failed, t_start)


if __name__ == "__main__":
    main()
