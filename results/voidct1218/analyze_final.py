#!/usr/bin/env python3
"""Final analysis for voidct1218.

1. Rescue pass: for sampled void windows whose 64-cube had n_pos < 500,
   pick the CT chunk overlapping the window with the most m7 positives
   (m7 from local cache), fetch that single chunk, redo stats on the
   overlap region.
2. Classification of all 60 void windows by the declared rule.
3. Pool chunk-existence pass: list all CT L1 chunks via S3 ListObjectsV2,
   then for every pool window (1174) count m7 positives falling in absent
   CT chunks (m7 read from the coverage run's local cache; no CT download).
4. Write windows_class.csv and final summary.json.

Rule (declared after inspecting the 90-window distributions, which are
bimodal in zfrac and tight in control nonzero mean):
  masked_empty (non-papyrus-like): zfrac_pos >= 0.9
  papyrus_like: zfrac_pos <= 0.5 and |nzmean_pos - MU_L| <= 2*SD_L
  ambiguous: otherwise
where zfrac_pos = fraction of m7-positive voxels with CT == 0 in the cube,
nzmean_pos = mean CT over the nonzero positives, MU_L/SD_L = mean/sd of the
30 control-window labeled-voxel nonzero means.
"""
import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from xml.etree import ElementTree

import numpy as np

sys.path.insert(0, "/mnt/vesuvius/overlap_step2")
from zarr_http import RemoteZarrLevel  # noqa: E402

OUT_DIR = "/mnt/vesuvius/voidct1218"
CT_URL = ("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/"
          "PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr")
M7_URL = ("https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc1218/"
          "representations/predictions/surfaces/"
          "20250521120456-surface-20260413222639-surface-m7-L0-th0.2.zarr")
CT_CACHE = os.path.join(OUT_DIR, "ct_cache")
M7_CACHE = "/mnt/vesuvius/hazard_zarr_smoke/m7L1_cache"
COV_CSV = "/mnt/vesuvius/labelcov1218/coverage_tiles.csv"
LIST_BASE = ("https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/"
             "?list-type=2&prefix=PHerc1218/volumes/"
             "20250521120456-8.640um-1.2m-116keV-masked.zarr/1/")
MIN_PRED = 100000
PRED_TH = 128
ZF_EMPTY = 0.9
ZF_LIVE = 0.5
NSD = 2.0
MIN_POS_CUBE = 500

t0 = time.time()


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def nz_stats(h):
    n = int(h.sum())
    if n == 0:
        return None, None
    zf = float(h[0]) / n
    nn = n - int(h[0])
    if nn == 0:
        return zf, None
    v = np.arange(256)
    return zf, float((h[1:] * v[1:]).sum() / nn)


H = np.load(os.path.join(OUT_DIR, "window_hists.npz"))
rows = list(csv.DictReader(open(os.path.join(OUT_DIR, "windows.csv"))))
cov = {r["path"]: r for r in csv.DictReader(open(COV_CSV))}


def hkey(path, name):
    return path.replace("/", "__").replace(".npz", "") + "__" + name


# ---- control reference ---------------------------------------------------
ctrl = [r for r in rows if r["group"] == "control"]
lab_nzm, lab_zf, unlab_nzm = [], [], []
for r in ctrl:
    zf, m = nz_stats(H[hkey(r["path"], "lab")])
    if m is not None:
        lab_zf.append(zf)
        lab_nzm.append(m)
    zfu, mu = nz_stats(H[hkey(r["path"], "unlab")])
    if mu is not None and int(r["n_pos_unlab"] or 0) >= 100:
        unlab_nzm.append(mu)
MU_L = float(np.mean(lab_nzm))
SD_L = float(np.std(lab_nzm, ddof=1))
log(f"control labeled ref: MU_L={MU_L:.1f} SD_L={SD_L:.1f} "
    f"max lab zfrac={max(lab_zf):.4f}; unlab-in-ctrl mean="
    f"{np.mean(unlab_nzm):.1f} sd={np.std(unlab_nzm, ddof=1):.1f} "
    f"(n={len(unlab_nzm)})")

# ---- rescue pass ---------------------------------------------------------
ct_lvl = RemoteZarrLevel(CT_URL, 1, cache_dir=CT_CACHE)
m7_lvl = RemoteZarrLevel(M7_URL, 1, cache_dir=M7_CACHE)
ct_bytes0 = sum(os.path.getsize(os.path.join(r_, f))
                for r_, _, fs in os.walk(CT_CACHE) for f in fs)

voidr = [r for r in rows if r["group"] == "void"]
need_rescue = [r for r in voidr if int(r["n_pos"]) < MIN_POS_CUBE]
log(f"rescue pass: {len(need_rescue)} windows")


def rescue(r):
    c = cov[r["path"]]
    lo = tuple(int(c[f"win_lo_{a}"]) for a in "zyx")
    hi = tuple(int(c[f"win_hi_{a}"]) for a in "zyx")
    m7w = m7_lvl.read_crop(lo, hi)
    pos = m7w >= PRED_TH
    best, best_n = None, -1
    for ci in range(lo[0] // 128, (hi[0] - 1) // 128 + 1):
        for cj in range(lo[1] // 128, (hi[1] - 1) // 128 + 1):
            for ck in range(lo[2] // 128, (hi[2] - 1) // 128 + 1):
                cl = (ci * 128, cj * 128, ck * 128)
                ov_lo = tuple(max(l, c_) for l, c_ in zip(lo, cl))
                ov_hi = tuple(min(h, c_ + 128) for h, c_ in zip(hi, cl))
                sl = tuple(slice(a - o, b - o)
                           for a, b, o in zip(ov_lo, ov_hi, lo))
                n = int(pos[sl].sum())
                if n > best_n:
                    best, best_n = (ov_lo, ov_hi), n
    ov_lo, ov_hi = best
    ctc = ct_lvl.read_crop(ov_lo, ov_hi)
    sl = tuple(slice(a - o, b - o) for a, b, o in zip(ov_lo, ov_hi, lo))
    p = pos[sl]
    n_pos = int(p.sum())
    if n_pos == 0:
        return r["path"], None, None, 0, ov_lo
    vals = ctc[p]
    h = np.bincount(vals, minlength=256)
    zf, nzm = nz_stats(h)
    return r["path"], zf, nzm, n_pos, ov_lo


rescued = {}
with ThreadPoolExecutor(max_workers=3) as ex:
    for path, zf, nzm, n_pos, ov_lo in ex.map(rescue, need_rescue):
        rescued[path] = (zf, nzm, n_pos, ov_lo)
        log(f"  rescue {path}: n_pos={n_pos} zf={zf} nzm={nzm} lo={ov_lo}")

# ---- classification ------------------------------------------------------
def classify(zf, nzm, n_pos):
    if n_pos == 0 or zf is None:
        return "no_positives_in_probe"
    if zf >= ZF_EMPTY:
        return "masked_empty"
    if zf <= ZF_LIVE and nzm is not None and abs(nzm - MU_L) <= NSD * SD_L:
        return "papyrus_like"
    return "ambiguous"


out_rows = []
for r in rows:
    zf, nzm = nz_stats(H[hkey(r["path"], "pos")])
    n_pos = int(r["n_pos"])
    rescue_used = 0
    if r["group"] == "void" and n_pos < MIN_POS_CUBE and r["path"] in rescued:
        zf, nzm, n_pos, _ = rescued[r["path"]]
        rescue_used = 1
    cls = classify(zf, nzm, n_pos) if r["group"] == "void" else ""
    out_rows.append(dict(r, zfrac_pos=("" if zf is None else round(zf, 4)),
                         nzmean_pos=("" if nzm is None else round(nzm, 1)),
                         n_pos_probe=n_pos, rescue=rescue_used, cls=cls))

with open(os.path.join(OUT_DIR, "windows_class.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
    w.writeheader()
    for o in out_rows:
        w.writerow(o)

vc = [o for o in out_rows if o["group"] == "void"]
counts = {}
for o in vc:
    counts[o["cls"]] = counts.get(o["cls"], 0) + 1
log(f"void classification: {counts}")

by_z, by_r = {}, {}
for o in vc:
    by_z.setdefault(o["zband"], {}).setdefault(o["cls"], 0)
    by_z[o["zband"]][o["cls"]] += 1
    by_r.setdefault(o["rthird"], {}).setdefault(o["cls"], 0)
    by_r[o["rthird"]][o["cls"]] += 1

# ---- pool chunk-existence pass ------------------------------------------
log("listing CT L1 chunks via S3 ...")
prefix = "PHerc1218/volumes/20250521120456-8.640um-1.2m-116keV-masked.zarr/1/"
existing = set()
token = None
n_req = 0
while True:
    url = LIST_BASE + "&max-keys=1000"
    if token:
        url += "&continuation-token=" + urllib.parse.quote(token)
    req = urllib.request.Request(url, headers={"User-Agent": "voidct-list"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        tree = ElementTree.fromstring(resp.read())
    n_req += 1
    ns = {"s3": tree.tag.split("}")[0].strip("{")}
    for k in tree.findall(".//s3:Contents/s3:Key", ns):
        rel = k.text[len(prefix):]
        if rel and not rel.startswith("."):
            existing.add(tuple(int(x) for x in rel.split("/")))
    trunc = tree.find("s3:IsTruncated", ns)
    if trunc is None or trunc.text != "true":
        break
    token = tree.find("s3:NextContinuationToken", ns).text
    time.sleep(0.1)
log(f"S3 listing: {len(existing)} existing chunks in {n_req} requests")

pool = [c for c in cov.values()
        if c["status"] == "ok" and c["clamped"] == "0"
        and int(c["n_pred"]) >= MIN_PRED]
log(f"pool pass over {len(pool)} windows (m7 from local cache)")


def pool_window(c):
    lo = tuple(int(c[f"win_lo_{a}"]) for a in "zyx")
    hi = tuple(int(c[f"win_hi_{a}"]) for a in "zyx")
    m7w = m7_lvl.read_crop(lo, hi)
    pos = m7w >= PRED_TH
    n_pos = int(pos.sum())
    n_absent = 0
    for ci in range(lo[0] // 128, (hi[0] - 1) // 128 + 1):
        for cj in range(lo[1] // 128, (hi[1] - 1) // 128 + 1):
            for ck in range(lo[2] // 128, (hi[2] - 1) // 128 + 1):
                if (ci, cj, ck) in existing:
                    continue
                cl = (ci * 128, cj * 128, ck * 128)
                ov_lo = tuple(max(l, x) for l, x in zip(lo, cl))
                ov_hi = tuple(min(h, x + 128) for h, x in zip(hi, cl))
                sl = tuple(slice(a - o, b - o)
                           for a, b, o in zip(ov_lo, ov_hi, lo))
                n_absent += int(pos[sl].sum())
    return c["path"], n_pos, n_absent, float(c["coverage_near4"])


pool_rows = []
with ThreadPoolExecutor(max_workers=3) as ex:
    for i, res in enumerate(ex.map(pool_window, pool)):
        pool_rows.append(res)
        if (i + 1) % 200 == 0:
            log(f"  pool {i + 1}/{len(pool)}")

tot_pos = sum(n for _, n, _, _ in pool_rows)
tot_absent = sum(a for _, _, a, _ in pool_rows)
void_pool = [(p, n, a) for p, n, a, cv in pool_rows if cv < 0.1]
ctrl_pool = [(p, n, a) for p, n, a, cv in pool_rows if cv > 0.98]
mid_pool = [(p, n, a) for p, n, a, cv in pool_rows if 0.1 <= cv <= 0.98]
pool_stats = {
    "n_windows": len(pool_rows),
    "s3_list_requests": n_req,
    "n_ct_chunks_existing_L1": len(existing),
    "total_m7_pos": tot_pos,
    "m7_pos_on_absent_ct_chunks": tot_absent,
    "frac_on_absent": round(tot_absent / tot_pos, 4),
    "void_windows(cov<0.1)": {
        "n": len(void_pool),
        "pos": sum(n for _, n, _ in void_pool),
        "on_absent": sum(a for _, _, a in void_pool),
        "frac": round(sum(a for _, _, a in void_pool)
                      / max(1, sum(n for _, n, _ in void_pool)), 4)},
    "control_windows(cov>0.98)": {
        "n": len(ctrl_pool),
        "pos": sum(n for _, n, _ in ctrl_pool),
        "on_absent": sum(a for _, _, a in ctrl_pool),
        "frac": round(sum(a for _, _, a in ctrl_pool)
                      / max(1, sum(n for _, n, _ in ctrl_pool)), 4)},
    "mid_windows(0.1..0.98)": {
        "n": len(mid_pool),
        "pos": sum(n for _, n, _ in mid_pool),
        "on_absent": sum(a for _, _, a in mid_pool),
        "frac": round(sum(a for _, _, a in mid_pool)
                      / max(1, sum(n for _, n, _ in mid_pool)), 4)},
}
log(f"pool: {tot_absent}/{tot_pos} = {tot_absent / tot_pos:.3f} of m7 "
    f"positives sit on absent CT chunks")

# ---- summary update ------------------------------------------------------
ct_bytes1 = sum(os.path.getsize(os.path.join(r_, f))
                for r_, _, fs in os.walk(CT_CACHE) for f in fs)
with open(os.path.join(OUT_DIR, "summary.json")) as f:
    summary = json.load(f)

examples = {}
for cls in ("masked_empty", "papyrus_like", "ambiguous",
            "no_positives_in_probe"):
    exs = []
    for o in vc:
        if o["cls"] != cls or len(exs) >= 5:
            continue
        c = cov[o["path"]]
        exs.append({
            "path": o["path"],
            "window_lo_zyx": [int(c["win_lo_z"]), int(c["win_lo_y"]),
                              int(c["win_lo_x"])],
            "cube_lo_zyx": [int(o["cube_lo_z"]), int(o["cube_lo_y"]),
                            int(o["cube_lo_x"])],
            "zband": int(o["zband"]), "rthird": int(o["rthird"]),
            "n_pos_probe": o["n_pos_probe"],
            "zfrac_pos": o["zfrac_pos"], "nzmean_pos": o["nzmean_pos"],
            "rescue": o["rescue"],
        })
    examples[cls] = exs

summary["analysis"] = {
    "rule": {
        "masked_empty": f"zfrac_pos >= {ZF_EMPTY}",
        "papyrus_like": f"zfrac_pos <= {ZF_LIVE} and |nzmean_pos - "
                        f"{MU_L:.1f}| <= {NSD}*{SD_L:.1f}",
        "ambiguous": "otherwise",
        "probe": f"64-cube if n_pos >= {MIN_POS_CUBE}, else rescue: the CT "
                 "chunk-overlap region of the window with most m7 positives",
        "note": "In this masked volume air and case are zeroed and fully "
                "empty chunks are absent (404 -> fill 0), so CT==0 means "
                "'no material in the canonical masked reconstruction'. "
                "zfrac_pos is the fraction of m7-positive voxels on CT==0.",
    },
    "control_reference": {
        "labeled_window_nzmeans_mean": round(MU_L, 1),
        "labeled_window_nzmeans_sd": round(SD_L, 1),
        "labeled_zfrac_max": round(max(lab_zf), 4),
        "unlab_in_control_nzmean_mean": round(float(np.mean(unlab_nzm)), 1),
        "unlab_in_control_nzmean_sd":
            round(float(np.std(unlab_nzm, ddof=1)), 1),
        "n_unlab_windows_used": len(unlab_nzm),
    },
    "void_classification_counts": counts,
    "void_by_zband": {str(k): v for k, v in sorted(by_z.items())},
    "void_by_radial_third": {str(k): v for k, v in sorted(by_r.items())},
    "examples": examples,
    "pool_chunk_existence": pool_stats,
    "rescue_windows": len(need_rescue),
    "extra_ct_download_bytes": ct_bytes1 - ct_bytes0,
    "analysis_runtime_s": round(time.time() - t0, 1),
    "limitations": [
        "Intensity is a partial discriminant, not ground truth: dense case "
        "material or minerals can mimic papyrus intensity, and the "
        "papyrus_like class only says the CT there is consistent with sheet "
        "material, not that a sheet is truly present.",
        "The volume is the canonical masked reconstruction; if the masking "
        "pipeline itself removed real faint papyrus, this method inherits "
        "that miss and would call it masked_empty.",
        "One 64-cube (or one chunk-overlap) probes each 128-window, not the "
        "whole window; the pool chunk-existence pass is chunk-granular "
        "(absent chunks only), so its fraction is a lower bound on "
        "positives without CT support.",
        "Windows are tile-center samples from the coverage run, not a "
        "whole-volume integral.",
    ],
}
with open(os.path.join(OUT_DIR, "summary.json"), "w") as f:
    json.dump(summary, f, indent=1)
log(f"analysis done in {time.time() - t0:.0f}s, extra CT dl "
    f"{(ct_bytes1 - ct_bytes0) / 1e6:.1f} MB")
