#!/usr/bin/env python3
"""Backfill evidence gallery: before/after figures from existing artifacts.

Writes 5 PNGs to /mnt/vesuvius/evidence_gallery/ and prints a JSON summary
(site choices, CT status, verification) on stdout between SUMMARY markers.

All inputs are read-only. Runs single-process; BLAS capped to 2 threads.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

import colorsys
import csv
import gzip
import json
import socket
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

socket.setdefaulttimeout(30)

OUT = "/mnt/vesuvius/evidence_gallery"
os.makedirs(OUT, exist_ok=True)

V3 = "/mnt/vesuvius/p1218_repair_v3/blocks_repaired"
V4 = "/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired"
SITES_CSV = "/mnt/vesuvius/kaggle_p1218_repair_v2/validation/v4_recast_sites.csv.gz"
PILOT = "/mnt/vesuvius/pilot0332"
SHEETS = f"{PILOT}/vesuvius-sheet-tools/output"
COVCSV = "/mnt/vesuvius/labelcov1218/coverage_tiles.csv"
ZONES = "/mnt/vesuvius/audit0139/thickness_map_1203/zones.npz"
DEMO = "/mnt/vesuvius/hazard_zarr_smoke/demo_sites.json"
OFFS = "/mnt/vesuvius/hazard_zarr_smoke/offseeds_placement_v2.json"

CT0332 = ("https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc0332/"
          "volumes/20251211183505-2.399um-0.2m-78keV-masked.zarr")

PHI = 0.6180339887498949
summary = {}

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "figure.facecolor": "white",
    "savefig.facecolor": "white",
})


def gr_color(i):
    h = (int(i) * PHI) % 1.0
    return colorsys.hsv_to_rgb(h, 0.65, 0.95)


def label_rgb(sl):
    out = np.zeros(sl.shape + (3,), np.float32)
    for uid in np.unique(sl):
        if uid == 0:
            continue
        out[sl == uid] = gr_color(uid)
    return out


def cross(ax, x, y, color="yellow", r=8):
    ax.plot([x - r, x + r], [y, y], color=color, lw=1.0)
    ax.plot([x, x], [y - r, y + r], color=color, lw=1.0)


def clean(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_edgecolor("0.7")
        s.set_linewidth(0.6)


def load_tile(tree, slab, tile):
    f = np.load(os.path.join(tree, slab, tile + ".npz"))
    return f["labels"], int(f["z0"]), int(f["y0"]), int(f["x0"])


# ---------------------------------------------------------------- FIG 1
def fig1():
    rows = []
    with gzip.open(SITES_CSV, "rt") as fh:
        for r in csv.DictReader(fh):
            if r["recast_ok"] != "1":
                continue
            rows.append({
                "slab": r["slab"], "tile": r["tile"],
                "z": int(r["z"]), "y": int(r["y"]), "x": int(r["x"]),
                "conf": float(r["conf"]),
                "applied": int(r["assigned_A"]) + int(r["assigned_B"]),
            })
    rows.sort(key=lambda r: (-r["applied"], r["slab"], r["tile"],
                             r["z"], r["y"], r["x"]))
    chosen, used_tiles, checked = [], set(), 0
    for r in rows:
        if len(chosen) == 3:
            break
        key = (r["slab"], r["tile"])
        if key in used_tiles:
            continue
        p3 = os.path.join(V3, r["slab"], r["tile"] + ".npz")
        p4 = os.path.join(V4, r["slab"], r["tile"] + ".npz")
        if not (os.path.exists(p3) and os.path.exists(p4)):
            continue
        checked += 1
        lab3, z0, y0, x0 = load_tile(V3, r["slab"], r["tile"])
        lab4, _, _, _ = load_tile(V4, r["slab"], r["tile"])
        s3, s4 = lab3[r["z"]].copy(), lab4[r["z"]].copy()
        del lab3, lab4
        fused = int(s3[r["y"], r["x"]])
        if fused == 0:  # nearest nonzero within 3 px
            yy0, yy1 = max(0, r["y"] - 3), min(512, r["y"] + 4)
            xx0, xx1 = max(0, r["x"] - 3), min(512, r["x"] + 4)
            nb = s3[yy0:yy1, xx0:xx1]
            vals = nb[nb > 0]
            if vals.size == 0:
                continue
            fused = int(np.bincount(vals).argmax())
        h = 80
        cy0, cy1 = max(0, r["y"] - h), min(512, r["y"] + h)
        cx0, cx1 = max(0, r["x"] - h), min(512, r["x"] + h)
        c3, c4 = s3[cy0:cy1, cx0:cx1], s4[cy0:cy1, cx0:cx1]
        region = c3 == fused
        ids4 = c4[region]
        ids4 = ids4[ids4 > 0]
        if ids4.size == 0:
            continue
        u, cnt = np.unique(ids4, return_counts=True)
        order = np.argsort(-cnt)
        if len(u) < 2 or cnt[order[1]] < 20:
            continue  # not a visually meaningful split on this slice
        used_tiles.add(key)
        chosen.append(dict(r, fused_id=fused, idA=int(u[order[0]]),
                           idB=int(u[order[1]]), crop=(cy0, cy1, cx0, cx1),
                           c3=c3, c4=c4, z0=z0, y0=y0, x0=x0))
    if len(chosen) < 3:
        raise RuntimeError(f"only {len(chosen)} usable split sites found")

    fig, axes = plt.subplots(3, 2, figsize=(13.0, 19.0), dpi=100)
    for i, site in enumerate(chosen):
        cy0, cy1, cx0, cx1 = site["crop"]
        c3, c4 = site["c3"], site["c4"]
        gz = site["z0"] + site["z"]
        gy = site["y0"] + site["y"]
        gx = site["x0"] + site["x"]
        # before: other labels gray, fused instance red
        img = np.zeros(c3.shape + (3,), np.float32)
        img[c3 > 0] = (0.42, 0.42, 0.42)
        img[c3 == site["fused_id"]] = (0.86, 0.30, 0.25)
        ax = axes[i, 0]
        ax.imshow(img, interpolation="nearest")
        cross(ax, site["x"] - cx0, site["y"] - cy0)
        ax.set_title(f"before (v3 labels)  site {i+1}: "
                     f"{site['slab']}/{site['tile']}", fontsize=9)
        ax.set_xlabel(f"fused instance id {site['fused_id']} (red), axial slice "
                      f"z={gz} (L1), crop {c3.shape[1]}x{c3.shape[0]} at "
                      f"y={gy} x={gx}", fontsize=7.5)
        clean(ax)
        # after: others gray, split ids A/B colored
        img = np.zeros(c4.shape + (3,), np.float32)
        img[c4 > 0] = (0.42, 0.42, 0.42)
        img[c4 == site["idA"]] = (0.25, 0.51, 0.77)
        img[c4 == site["idB"]] = (0.95, 0.61, 0.19)
        ax = axes[i, 1]
        ax.imshow(img, interpolation="nearest")
        cross(ax, site["x"] - cx0, site["y"] - cy0)
        ax.set_title("after (v2.0 repair, v4 labels)", fontsize=9)
        ax.set_xlabel(f"split ids {site['idA']} (blue) / {site['idB']} (orange), "
                      f"applied voxels {site['applied']}, conf {site['conf']:.3f}, "
                      f"recast_ok=1", fontsize=7.5)
        clean(ax)
    fig.suptitle("PHerc1218 v2.0 topological repair: fused-sheet sites before/after "
                 "(3 validated splits, largest applied-voxel counts)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    p = f"{OUT}/repair_1218_before_after.png"
    fig.savefig(p)
    plt.close(fig)
    summary["fig1"] = {
        "path": p, "n_candidates_checked": checked,
        "sites": [{k: s[k] for k in ("slab", "tile", "z", "y", "x", "conf",
                                     "applied", "fused_id", "idA", "idB")}
                  for s in chosen],
    }


# ---------------------------------------------------------------- FIG 2
def fig2():
    grid = json.load(open(f"{PILOT}/crop_grid.json"))
    scored = []
    for c in grid["crops"]:
        cz, cy, cx = c["centre_l1"]
        f = f"{SHEETS}/sheets_scroll3_L1_c{cz}-{cy}-{cx}.npy"
        if not os.path.exists(f):
            continue
        sl = np.load(f, mmap_mode="r")[256]
        frac = float((np.asarray(sl) > 0).mean())
        scored.append((-frac, os.path.basename(f), c, f))
    scored.sort(key=lambda t: (t[0], t[1]))
    picked = [(c, f) for _, _, c, f in scored[:2]]
    summary["fig2_selection"] = [
        {"file": b, "midz_nonzero_frac": round(-s, 4)}
        for s, b, _, _ in scored[:2]]
    ct_status = {}
    ct_slices = {}
    try:
        import zarr
        ct = zarr.open_array(f"{CT0332}/3", mode="r")
        for c, f in picked:
            cz, cy, cx = c["centre_l1"]
            ys, xs = cy - 256, cx - 256
            z = min(max(cz, 0), ct.shape[0] - 1)
            sl = np.asarray(ct[z, max(0, ys):ys + 512, max(0, xs):xs + 512]
                            ).astype(np.float32)
            if sl.size and sl.max() > 0:
                ct_slices[f] = sl
                ct_status[f] = f"ok CT L3 z={z} from {CT0332}"
            else:
                ct_status[f] = "empty slice, fell back to mask silhouette"
    except Exception as e:  # noqa: BLE001
        for _, f in picked:
            ct_status[f] = f"CT fetch failed ({type(e).__name__}: {e}); silhouette"

    fig, axes = plt.subplots(2, 2, figsize=(14.0, 14.6), dpi=100)
    for i, (c, f) in enumerate(picked):
        cz, cy, cx = c["centre_l1"]
        lab = np.load(f, mmap_mode="r")[256].copy()  # mid-z slice
        ax = axes[i, 0]
        if f in ct_slices:
            sl = ct_slices[f]
            lo, hi = np.percentile(sl, (1, 99))
            ax.imshow(np.clip((sl - lo) / max(hi - lo, 1e-6), 0, 1),
                      cmap="gray", interpolation="nearest")
            left_cap = ("streamed CT slice, masked volume L3 z={} "
                        "(grid matches pred L1 to ~px)".format(cz))
        else:
            ax.imshow((lab > 0).astype(np.float32), cmap="gray",
                      interpolation="nearest")
            left_cap = "label>0 silhouette (CT not fetched)"
        ax.set_title(f"crop {i+1}: centre L1 z={cz} y={cy} x={cx} "
                     f"(fill {c['fill_l2']:.3f})", fontsize=9)
        ax.set_xlabel(left_cap, fontsize=7.5)
        clean(ax)
        ax = axes[i, 1]
        ax.imshow(label_rgb(lab), interpolation="nearest")
        nid = len(np.unique(lab)) - 1
        ax.set_title("separated sheet instances, same slice", fontsize=9)
        ax.set_xlabel(f"{nid} instance ids in slice, golden-ratio hue palette, "
                      "mid-z axial slice of 512^3 crop", fontsize=7.5)
        clean(ax)
    fig.suptitle("PHerc0332 sheet separation: 512^3 L1 crops with densest "
                 "mid-z slices (left: input context, right: instances)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.98), h_pad=2.5)
    p = f"{OUT}/separation_0332_slices.png"
    fig.savefig(p)
    plt.close(fig)
    summary["fig2"] = {
        "path": p,
        "crops": [{"centre_l1": c["centre_l1"], "fill_l2": c["fill_l2"],
                   "file": f, "ct": ct_status[f]} for c, f in picked],
    }


# ---------------------------------------------------------------- FIG 3
def fig3():
    rows = []
    with open(COVCSV) as fh:
        for r in csv.DictReader(fh):
            if r["status"] != "ok" or int(r["n_pred"]) < 1000:
                continue
            rows.append((int(r["z0"]), int(r["y0"]), int(r["x0"]),
                         float(r["coverage_near4"])))
    z0s = np.array([r[0] for r in rows])
    zmax = 11624
    bands = [(0, zmax // 3), (zmax // 3, 2 * zmax // 3), (2 * zmax // 3, zmax + 1)]
    names = ["low z", "mid z", "high z"]
    ys = sorted({r[1] for r in rows})
    xs = sorted({r[2] for r in rows})
    yi = {v: i for i, v in enumerate(ys)}
    xi = {v: i for i, v in enumerate(xs)}

    fig, axes = plt.subplots(2, 2, figsize=(14.0, 12.4), dpi=100,
                             layout="constrained")
    im = None
    for b, (blo, bhi) in enumerate(bands):
        acc = np.zeros((len(ys), len(xs)), np.float64)
        n = np.zeros((len(ys), len(xs)), np.int32)
        for z0, y0, x0, cov in rows:
            if blo <= z0 < bhi:
                acc[yi[y0], xi[x0]] += cov
                n[yi[y0], xi[x0]] += 1
        mean = np.where(n > 0, acc / np.maximum(n, 1), np.nan)
        ax = axes.flat[b]
        im = ax.imshow(np.ma.masked_invalid(mean), cmap="viridis",
                       vmin=0, vmax=1, interpolation="nearest",
                       extent=(xs[0], xs[-1] + 448, ys[-1] + 448, ys[0]))
        ax.set_title(f"{names[b]}: z0 in [{blo}, {bhi})  "
                     f"({int((z0s >= blo).sum() - (z0s >= bhi).sum())} tiles)",
                     fontsize=9)
        ax.set_xlabel(f"tile x0+256 (L1 vox); mean coverage_near4 over slabs "
                      f"in band; white = no tile", fontsize=7.5)
        ax.set_ylabel("tile y0+256 (L1 vox)", fontsize=7.5)
        ax.tick_params(labelsize=7)
    cb = fig.colorbar(im, ax=[axes[0, 0], axes[0, 1], axes[1, 0]],
                      shrink=0.55, location="right", pad=0.01)
    cb.set_label("coverage_near4", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax = axes.flat[3]
    covs = [r[3] for r in rows]
    ax.hist(covs, bins=np.linspace(0, 1, 51), color="0.35")
    ax.set_title(f"per-tile coverage_near4 histogram ({len(covs)} tiles, "
                 "n_pred >= 1000)", fontsize=9)
    ax.set_xlabel("coverage_near4 (fraction of predicted voxels within "
                  "Chebyshev 4 of a label)", fontsize=7.5)
    ax.set_ylabel("tiles", fontsize=7.5)
    ax.tick_params(labelsize=7)
    fig.suptitle("PHerc1218 label coverage vs m7 surface prediction: per-tile "
                 "128^3 center windows (outer-wrap void ring visible)", fontsize=10)
    p = f"{OUT}/labelcov_1218_map.png"
    fig.savefig(p)
    plt.close(fig)
    summary["fig3"] = {"path": p, "n_tiles": len(covs), "bands": bands}


# ---------------------------------------------------------------- FIG 4
def fig4():
    z = np.load(ZONES)
    fs = z["fused_share"]
    nr = z["n_runs"]
    vs = z["valid_share"]
    zsum = np.nan_to_num(vs, nan=0.0).sum(axis=(1, 2))
    zsel = sorted(np.argsort(-zsum)[:6].tolist())
    edge_runs = float(np.nansum(nr[:, 35, :]) + np.nansum(nr[:, :, 35]))

    cmap = matplotlib.colormaps["viridis"].copy()
    cmap.set_bad("0.82")
    fig, axes = plt.subplots(2, 3, figsize=(14.0, 10.3), dpi=100,
                             layout="constrained")
    im = None
    for i, zi in enumerate(zsel):
        ax = axes.flat[i]
        sl = np.ma.masked_where(nr[zi] == 0, fs[zi])
        im = ax.imshow(sl, cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
        ax.add_patch(Rectangle((34.52, -0.48), 0.96, 35.96, fill=False,
                               edgecolor="red", lw=0.9, clip_on=False))
        ax.add_patch(Rectangle((-0.48, 34.52), 35.96, 0.96, fill=False,
                               edgecolor="red", lw=0.9, clip_on=False))
        ax.set_title(f"z-cell {zi} (L1 z {3936 + zi * 96}-{3936 + (zi + 1) * 96})",
                     fontsize=9)
        ax.set_xlabel(f"{int((nr[zi] > 0).sum())} sampled cells",
                      fontsize=7.5)
        clean(ax)
    cb = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.7, pad=0.01)
    cb.set_label("fused_share (fraction of thick single-crossing runs)",
                 fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.suptitle("PHerc1203 thickness map v2 (audit0139): fused_share on the "
                 "21x36x36 zone grid, 6 z-cells with most valid data.\n"
                 "Gray = 0 sampled runs; red outline = dead edge row/col 35 "
                 "(grid-edge padding, 0 runs everywhere)", fontsize=10)
    p = f"{OUT}/thickness1203_v2_map.png"
    fig.savefig(p)
    plt.close(fig)
    summary["fig4"] = {"path": p, "z_cells": zsel,
                       "edge35_total_runs": edge_runs}


# ---------------------------------------------------------------- FIG 5
def fig5():
    demo = json.load(open(DEMO))
    offs = json.load(open(OFFS))
    assert len(demo) == len(offs) == 8, (len(demo), len(offs))
    panels = []  # (kind, idx, slab, tile, lz, ly, lx, gz, gy, gx, extra)
    tiles = {}
    for i, (d, o) in enumerate(zip(demo, offs)):
        assert d["slab"] == o["slab"] and d["tile"] == o["tile"], i
        key = (d["slab"], d["tile"])
        if key not in tiles:
            lab, z0, y0, x0 = load_tile(V3, *key)
            tiles[key] = (lab, z0, y0, x0)
        lab, z0, y0, x0 = tiles[key]
        panels.append(("ON", i, key, d["gz"] - z0, d["gy"] - y0, d["gx"] - x0,
                       d["gz"], d["gy"], d["gx"], f"n_sites={d['n_sites']}"))
        panels.append(("OFF", i, key, o["gz"] - z0, o["gy"] - y0, o["gx"] - x0,
                       o["gz"], o["gy"], o["gx"],
                       f"cheb>={o['min_chebyshev_to_cluster']}"))

    fig, axes = plt.subplots(4, 4, figsize=(15.0, 16.4), dpi=100)
    fracs = []
    half = 112
    for p_i, (kind, i, key, lz, ly, lx, gz, gy, gx, extra) in enumerate(panels):
        lab = tiles[key][0]
        lz = min(max(lz, 0), 255)
        sl = lab[lz]
        cy0, cy1 = max(0, ly - half), min(512, ly + half)
        cx0, cx1 = max(0, lx - half), min(512, lx + half)
        crop = sl[cy0:cy1, cx0:cx1]
        frac = float((crop > 0).mean())
        fracs.append({"site": i + 1, "kind": kind, "label_px_frac": round(frac, 4)})
        ax = axes.flat[p_i]
        ax.imshow(label_rgb(crop), interpolation="nearest")
        cross(ax, lx - cx0, ly - cy0, r=10)
        ax.set_title(f"S{i+1} {kind} seed", fontsize=9)
        ax.set_xlabel(f"z={gz} y={gy} x={gx} (L1)\nlabels {100*frac:.1f}% of px, "
                      f"{extra}", fontsize=7)
        clean(ax)
    for lab, *_ in tiles.values():
        del lab
    fig.suptitle("PHerc1218 placement A/B: ON seeds sit in labeled sheet "
                 "structure, OFF seeds (placement v2) sit in label-void space; "
                 "axial v3-label slices, 224^2 crops, seed at crosshair",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    p = f"{OUT}/placement_offseed_context.png"
    fig.savefig(p)
    plt.close(fig)
    summary["fig5"] = {"path": p, "panels": fracs}


# ---------------------------------------------------------------- run + verify
def main():
    todo = sys.argv[1:] or ["fig1", "fig2", "fig3", "fig4", "fig5"]
    for fn in (fig1, fig2, fig3, fig4, fig5):
        if fn.__name__ not in todo:
            continue
        try:
            fn()
            print(f"[done] {fn.__name__}", file=sys.stderr)
        except Exception as e:  # noqa: BLE001
            import traceback
            traceback.print_exc()
            summary[fn.__name__] = {"error": f"{type(e).__name__}: {e}"}
    from PIL import Image
    ver = []
    for name in sorted(os.listdir(OUT)):
        if not name.endswith(".png"):
            continue
        fp = os.path.join(OUT, name)
        im = Image.open(fp)
        arr = np.asarray(im.convert("L"), np.float32)
        ver.append({"file": name, "bytes": os.path.getsize(fp),
                    "size": list(im.size), "std": round(float(arr.std()), 2)})
    summary["verify"] = ver
    print("===SUMMARY===")
    print(json.dumps(summary, indent=1))
    print("===END===")


if __name__ == "__main__":
    main()
