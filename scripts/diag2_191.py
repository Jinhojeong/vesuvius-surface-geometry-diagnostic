"""#191 diagnostic v2 — validated harness (surface_recto_059_redo, dice~0.8 baseline).

Per-patch: inference -> per-voxel recall stratified by
  (a) curvature = meso-scale normal dispersion (1 - |mean normal| in r-ball),
  (b) compression = distance to the NEXT sheet along +/- normal (voxels).
Both computed from GT geometry only. Plus per-scroll breakdown and worst gallery.
"""
import os, glob, sys, json, time, csv
import numpy as np, tifffile, torch
from scipy import ndimage as ndi

OUT = "/mnt/vesuvius/surf191_rand"
DIAG = os.path.join(OUT, "diag2"); os.makedirs(DIAG, exist_ok=True)
CKPT = "/mnt/vesuvius/models/surface_recto_059_redo/Model_epoch499.pth"
TH = 0.4

sys.path.insert(0, "/mnt/vesuvius")
from loader059 import load_059, predict  # validated loader + ct-norm sliding window

def geometry(gt):
    """normals from smoothed SDF; sign-invariant orientation moments for curvature."""
    d_out = ndi.distance_transform_edt(~gt); d_in = ndi.distance_transform_edt(gt)
    sdf = ndi.gaussian_filter((d_out - d_in).astype(np.float32), 3.0)
    gz, gy, gx = np.gradient(sdf)
    mag = np.sqrt(gz**2 + gy**2 + gx**2) + 1e-6
    n = np.stack([gz/mag, gy/mag, gx/mag])          # [3,Z,Y,X] unit normals
    # SDF gradient flips sign across the sheet midplane, so a mean-vector
    # dispersion saturates at 1. Use the structure tensor of orientations
    # instead: smoothed second moments of n (sign-invariant).
    idx = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
    mom = np.stack([ndi.uniform_filter(n[i]*n[j], size=13) for i, j in idx])
    return n, mom

def curvature_at(points, mom):
    """1 - lambda_max of the local orientation tensor at [K,3] points.
    0 = coherent single orientation (flat sheet), higher = bent/complex."""
    pz, py, px = points[:, 0], points[:, 1], points[:, 2]
    m = [mom[c][pz, py, px] for c in range(6)]
    T = np.empty((len(points), 3, 3), np.float32)
    T[:, 0, 0], T[:, 1, 1], T[:, 2, 2] = m[0], m[1], m[2]
    T[:, 0, 1] = T[:, 1, 0] = m[3]
    T[:, 0, 2] = T[:, 2, 0] = m[4]
    T[:, 1, 2] = T[:, 2, 1] = m[5]
    lam = np.linalg.eigvalsh(T)[:, -1]              # largest eigenvalue
    return (1.0 - lam).clip(0, 1)

def spacing_at(points, gt, n):
    """next-sheet distance along +/-normal for [K,3] integer points."""
    Z, Y, X = gt.shape
    K = len(points)
    sp = np.full(K, 40.0, np.float32)   # cap = 40 vox (no neighbor found)
    pz, py, px = points[:, 0], points[:, 1], points[:, 2]
    nz, ny, nx = n[0][pz, py, px], n[1][pz, py, px], n[2][pz, py, px]
    for sign in (+1, -1):
        left_own = np.zeros(K, bool); found = np.zeros(K, bool)
        for t in range(1, 40):
            qz = np.clip((pz + sign*t*nz).round().astype(int), 0, Z-1)
            qy = np.clip((py + sign*t*ny).round().astype(int), 0, Y-1)
            qx = np.clip((px + sign*t*nx).round().astype(int), 0, X-1)
            hit = gt[qz, qy, qx]
            left_own |= ~hit
            new = left_own & hit & ~found
            sp[new] = np.minimum(sp[new], t); found |= new
            if found.all(): break
    return sp

def main():
    imgs = sorted(glob.glob(os.path.join(OUT, "imagesTr", "*_0000.tif")))
    print(f"{len(imgs)} patches", flush=True)
    net, norm, props = load_059()
    all_cv = []; all_sp = []; all_hit = []
    rows = []; worst = []
    rng = np.random.default_rng(0)
    for i, f in enumerate(imgs):
        lf = os.path.join(OUT, "labelsTr", os.path.basename(f).replace("_0000.tif", ".tif"))
        if not os.path.exists(lf): continue
        try:
            img = tifffile.imread(f); lbl = tifffile.imread(lf)
        except Exception as e:
            print("skip unreadable", os.path.basename(f), e); continue
        if img.ndim != 3 or lbl.shape != img.shape: continue
        gt = lbl > 0
        if gt.sum() < 5000: continue
        t0 = time.time()
        prob = predict(net, img, 192, norm, props)
        pred = prob > TH
        inter = (pred & gt).sum()
        dice = 2*inter/(pred.sum()+gt.sum()+1e-6); rec = inter/gt.sum()
        n, mom = geometry(gt)
        # sample GT voxels for stratification (speed)
        idx = np.argwhere(gt)
        take = rng.choice(len(idx), size=min(30000, len(idx)), replace=False)
        pts = idx[take]
        cv = curvature_at(pts, mom)
        sp = spacing_at(pts, gt, n)
        hit = pred[pts[:, 0], pts[:, 1], pts[:, 2]]
        all_cv.append(cv); all_sp.append(sp); all_hit.append(hit)
        scroll = os.path.basename(f).split("_")[0]
        rows.append((os.path.basename(f), scroll, float(dice), float(rec),
                     float(np.median(cv)), float(np.median(sp))))
        worst.append((dice, f, lf))
        print(f"[{i+1}/{len(imgs)}] {os.path.basename(f)[:30]:30s} {scroll} "
              f"dice={dice:.3f} rec={rec:.3f} medcurv={np.median(cv):.3f} "
              f"medsp={np.median(sp):.1f} ({time.time()-t0:.0f}s)", flush=True)
    with open(os.path.join(DIAG, "per_patch.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["patch","scroll","dice","recall","med_curv","med_spacing"])
        w.writerows(rows)
    cv = np.concatenate(all_cv); sp = np.concatenate(all_sp)
    hit = np.concatenate(all_hit)
    def table(vals, bins, name):
        print(f"\n=== recall by {name} ==="); out = []
        for b in range(len(bins)-1):
            m = (vals >= bins[b]) & (vals < bins[b+1])
            nn = int(m.sum()); r = float(hit[m].mean()) if nn else float("nan")
            out.append([float(bins[b]), float(bins[b+1]), nn, r])
            print(f"  [{bins[b]:.3f},{bins[b+1]:.3f})  n={nn:>9,}  recall={r:.3f}")
        return out
    # global decile bins for curvature (distribution unknown a priori)
    cbins = np.unique(np.quantile(cv, np.linspace(0, 1, 11))); cbins[-1] += 1e-6
    sbins = np.array([0, 4, 6, 8, 11, 15, 22, 39, 41])
    ct = table(cv, cbins, "curvature (orientation-tensor dispersion, decile bins)")
    st = table(sp, sbins, "spacing to next sheet (vox; 39-41=none found)")
    # per-scroll summary
    print("\n=== per-scroll dice ===")
    import collections
    g = collections.defaultdict(list)
    for _, s, d, *_ in rows: g[s].append(d)
    scr = {s: [float(np.mean(v)), float(np.median(v)), len(v)] for s, v in g.items()}
    for s, (m, md, k) in sorted(scr.items()): print(f"  {s}: mean={m:.3f} median={md:.3f} n={k}")
    json.dump({"curvature": ct, "spacing": st, "per_scroll": scr, "th": TH},
              open(os.path.join(DIAG, "strata2.json"), "w"), indent=1)
    # plots + worst gallery
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, t, nm in [(axes[0], ct, "curvature (normal dispersion)"),
                      (axes[1], st, "spacing to next sheet (vox)")]:
        xs = range(len(t))
        ax.bar(xs, [r[3] for r in t], color="#33638d")
        ax.set_xticks(list(xs))
        ax.set_xticklabels([f"{r[0]:g}-{r[1]:g}" for r in t], rotation=45, ha="right", fontsize=7)
        ax.set_ylim(0, 1); ax.set_ylabel(f"voxel recall @th={TH}"); ax.set_xlabel(nm)
    fig.tight_layout(); fig.savefig(os.path.join(DIAG, "strata2.png"), dpi=130)
    worst.sort()
    fig, axes = plt.subplots(3, 4, figsize=(14, 10))
    for j, (d, f, lf) in enumerate(worst[:4]):
        img = tifffile.imread(f); lbl = tifffile.imread(lf)
        prob = predict(net, img, 192, norm, props)
        z = img.shape[0]//2
        axes[0][j].imshow(img[z], cmap="gray"); axes[0][j].set_title(f"{os.path.basename(f)[:18]}\ndice={d:.2f}", fontsize=8)
        axes[1][j].imshow(lbl[z], cmap="viridis")
        axes[2][j].imshow(prob[z], cmap="magma", vmin=0, vmax=1)
        for r in range(3): axes[r][j].axis("off")
    fig.tight_layout(); fig.savefig(os.path.join(DIAG, "worst_gallery.png"), dpi=120)
    print("\nDIAG2 COMPLETE ->", DIAG, flush=True)

if __name__ == "__main__":
    main()
