"""Label-completeness QC for segmentation-derived surface GT (villa #192/#193).

Finds CANDIDATE UNLABELED SHEETS: high-confidence model predictions that are
(a) >3 vox outside the GT, (b) sheet-shaped (planar PCA), (c) CT-bright.
Outputs per-patch incompleteness score, correlation with dice, worst overlays.
"""
import os, glob, sys, json, time, csv
import numpy as np, tifffile, torch
from scipy import ndimage as ndi

OUT = "/mnt/vesuvius/surf191_rand"
QC = os.path.join(OUT, "qc192"); os.makedirs(QC, exist_ok=True)
TH_CONF = 0.65          # high-precision operating point
MIN_COMP = 300          # voxels
PLANARITY_MIN = 0.55    # (l1-l2)/l1 high & thin: use thin-ness below
THIN_MAX = 0.25         # l3/l1 small => thin sheet-like
CT_PCT_MIN = 52         # component mean CT percentile must exceed background-ish

sys.path.insert(0, "/mnt/vesuvius")
from loader059 import load_059, predict

def sheetlike(comp_pts):
    """PCA eigenvalue shape test for a component's coordinates [K,3]."""
    if len(comp_pts) < 30: return False, 0.0
    c = comp_pts - comp_pts.mean(0)
    cov = (c.T @ c) / len(c)
    lam = np.linalg.eigvalsh(cov)[::-1]  # l1 >= l2 >= l3
    thin = lam[2] / (lam[0] + 1e-9)
    return thin < THIN_MAX, float(thin)

def main():
    imgs = sorted(glob.glob(os.path.join(OUT, "imagesTr", "*_0000.tif")))
    net, norm, props = load_059()
    # dice per patch from diag2 (for the score<->dice correlation)
    dice_map = {}
    with open(os.path.join(OUT, "diag2", "per_patch.csv")) as fh:
        for r in csv.DictReader(fh):
            dice_map[r["patch"]] = float(r["dice"])
    rows = []; worst = []
    for i, f in enumerate(imgs):
        base = os.path.basename(f)
        lf = os.path.join(OUT, "labelsTr", base.replace("_0000.tif", ".tif"))
        if not os.path.exists(lf): continue
        try:
            img = tifffile.imread(f); lbl = tifffile.imread(lf)
        except Exception:
            continue
        if img.ndim != 3 or lbl.shape != img.shape: continue
        gt = lbl > 0
        if gt.sum() < 5000: continue
        t0 = time.time()
        prob = predict(net, img, 192, norm, props)
        conf = prob > TH_CONF
        edt = ndi.distance_transform_edt(~gt)
        cand = conf & (edt > 3)
        # CT percentile map
        qs = np.quantile(img, np.linspace(0, 1, 101))
        pct = np.searchsorted(qs, img).astype(np.float32)
        labels_arr, ncomp = ndi.label(cand)
        accepted = np.zeros_like(gt)
        n_acc = 0
        if ncomp:
            sizes = np.bincount(labels_arr.ravel())
            for cid in np.nonzero(sizes >= MIN_COMP)[0]:
                if cid == 0: continue
                m = labels_arr == cid
                pts = np.argwhere(m)
                ok_shape, thin = sheetlike(pts[:: max(1, len(pts)//4000)])
                bright = pct[m].mean()
                if ok_shape and bright >= CT_PCT_MIN:
                    accepted |= m; n_acc += 1
        miss_vox = int(accepted.sum())
        score = miss_vox / (gt.sum() + miss_vox + 1e-9)  # incompleteness fraction
        d = dice_map.get(base, float("nan"))
        rows.append((base, float(score), miss_vox, n_acc, d))
        worst.append((score, f, lf, accepted.copy() if score > 0.05 else None))
        print(f"[{i+1}/{len(imgs)}] {base[:30]:30s} incompleteness={score:.3f} "
              f"({n_acc} sheets, {miss_vox:,} vox) dice={d:.3f} ({time.time()-t0:.0f}s)",
              flush=True)
    with open(os.path.join(QC, "incompleteness.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["patch", "incompleteness", "missed_vox", "n_sheets", "dice_059redo"])
        w.writerows(rows)
    sc = np.array([r[1] for r in rows]); dc = np.array([r[4] for r in rows])
    ok = ~np.isnan(dc)
    r_pear = float(np.corrcoef(sc[ok], dc[ok])[0, 1])
    from scipy import stats
    r_spear = float(stats.spearmanr(sc[ok], dc[ok]).statistic)
    print(f"\n=== SCORE<->DICE correlation (n={ok.sum()}) ===")
    print(f"pearson r = {r_pear:.3f}   spearman rho = {r_spear:.3f}")
    print(f"mean incompleteness = {sc.mean():.3f}  median = {np.median(sc):.3f}  "
          f"p90 = {np.quantile(sc, .9):.3f}")
    json.dump({"pearson": r_pear, "spearman": r_spear,
               "mean": float(sc.mean()), "median": float(np.median(sc)),
               "p90": float(np.quantile(sc, .9)), "th_conf": TH_CONF},
              open(os.path.join(QC, "summary.json"), "w"), indent=1)
    # overlays of the 4 worst
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    worst.sort(key=lambda x: -x[0])
    fig, axes = plt.subplots(3, 4, figsize=(14, 10))
    for j, (s, f, lf, acc) in enumerate(worst[:4]):
        img = tifffile.imread(f); lbl = tifffile.imread(lf)
        z = img.shape[0] // 2
        axes[0][j].imshow(img[z], cmap="gray")
        axes[0][j].set_title(f"{os.path.basename(f)[:18]}\nincompleteness={s:.2f}", fontsize=8)
        axes[1][j].imshow(lbl[z], cmap="viridis"); axes[1][j].set_ylabel("GT")
        ov = np.zeros((*img[z].shape, 3), np.uint8)
        ov[..., 1] = (lbl[z] > 0) * 200                    # green = GT
        if acc is not None: ov[..., 0] = acc[z] * 255      # red = candidate unlabeled sheet
        axes[2][j].imshow(ov)
        for r_ in range(3): axes[r_][j].axis("off")
    fig.tight_layout(); fig.savefig(os.path.join(QC, "worst_incomplete.png"), dpi=120)
    print("\nQC192 COMPLETE ->", QC, flush=True)

if __name__ == "__main__":
    main()
