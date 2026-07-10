"""#191 diag4 — threshold-free two-model comparison: per-geometry-bin AUC.

For each patch: sample 30k GT surface points + 30k near-sheet background points
(2-10 vox from the nearest sheet). Record model prob at each point for BOTH
models (059_redo and nnUNet surface_recto). Per curvature/spacing bin, AUC =
P(prob_gt > prob_bg) between GT and background points of that bin (background
points inherit geometry from their nearest GT voxel). Firing-rate independent.
"""
import os, glob, sys, json, time
import numpy as np, tifffile, torch
from scipy import ndimage as ndi

OUT = "/mnt/vesuvius/surf191_rand"
DIAG = os.path.join(OUT, "diag4"); os.makedirs(DIAG, exist_ok=True)
K = 30000

sys.path.insert(0, "/mnt/vesuvius")
from diag2_191 import geometry, curvature_at, spacing_at
from loader059 import load_059, predict as predict059
from diag3_d058 import make_predictor

def auc(pos, neg):
    """rank-based AUC, robust to ties."""
    if len(pos) == 0 or len(neg) == 0: return float("nan")
    allv = np.concatenate([pos, neg])
    order = allv.argsort(kind="mergesort")
    ranks = np.empty(len(allv)); ranks[order] = np.arange(1, len(allv)+1)
    # average ranks for ties
    sv = allv[order]
    i = 0
    while i < len(sv):
        j = i
        while j+1 < len(sv) and sv[j+1] == sv[i]: j += 1
        if j > i: ranks[order[i:j+1]] = ranks[order[i:j+1]].mean()
        i = j+1
    rp = ranks[:len(pos)].sum()
    return float((rp - len(pos)*(len(pos)+1)/2) / (len(pos)*len(neg)))

def main():
    imgs = sorted(glob.glob(os.path.join(OUT, "imagesTr", "*_0000.tif")))
    print(f"{len(imgs)} patches", flush=True)
    netA, norm, props = load_059()
    predB = make_predictor()
    nnprops = {"spacing": [1.0, 1.0, 1.0]}
    rng = np.random.default_rng(0)
    recs = {"curv": [], "sp": [], "is_gt": [], "pA": [], "pB": []}
    for i, f in enumerate(imgs):
        lf = os.path.join(OUT, "labelsTr", os.path.basename(f).replace("_0000.tif", ".tif"))
        if not os.path.exists(lf): continue
        try:
            img = tifffile.imread(f); lbl = tifffile.imread(lf)
        except Exception:
            continue
        if img.ndim != 3 or lbl.shape != img.shape: continue
        gt = lbl > 0
        if gt.sum() < 5000: continue
        t0 = time.time()
        probA = predict059(netA, img, 192, norm, props)
        probB = predB.predict_single_npy_array(img[None].astype(np.float32),
                                               nnprops, None, None, True)[1][1]
        n, mom = geometry(gt)
        # GT sample
        idx = np.argwhere(gt)
        take = rng.choice(len(idx), size=min(K, len(idx)), replace=False)
        pts_gt = idx[take]
        # near-sheet background sample: 2-10 vox from nearest sheet
        edt, inds = ndi.distance_transform_edt(~gt, return_indices=True)
        shell = (edt >= 2) & (edt <= 10)
        idxb = np.argwhere(shell)
        takeb = rng.choice(len(idxb), size=min(K, len(idxb)), replace=False)
        pts_bg = idxb[takeb]
        # background points inherit geometry from nearest GT voxel
        nz = inds[0][pts_bg[:, 0], pts_bg[:, 1], pts_bg[:, 2]]
        ny = inds[1][pts_bg[:, 0], pts_bg[:, 1], pts_bg[:, 2]]
        nx = inds[2][pts_bg[:, 0], pts_bg[:, 1], pts_bg[:, 2]]
        anchor = np.stack([nz, ny, nx], 1)
        cv_gt = curvature_at(pts_gt, mom); cv_bg = curvature_at(anchor, mom)
        sp_gt = spacing_at(pts_gt, gt, n); sp_bg = spacing_at(anchor, gt, n)
        for pts, cv, sp, isg in ((pts_gt, cv_gt, sp_gt, 1), (pts_bg, cv_bg, sp_bg, 0)):
            recs["curv"].append(cv); recs["sp"].append(sp)
            recs["is_gt"].append(np.full(len(pts), isg, np.uint8))
            recs["pA"].append(probA[pts[:, 0], pts[:, 1], pts[:, 2]])
            recs["pB"].append(probB[pts[:, 0], pts[:, 1], pts[:, 2]])
        print(f"[{i+1}/{len(imgs)}] {os.path.basename(f)[:28]:28s} ({time.time()-t0:.0f}s)",
              flush=True)
    d = {k: np.concatenate(v) for k, v in recs.items()}
    np.savez_compressed(os.path.join(DIAG, "points.npz"), **d)
    isg = d["is_gt"] == 1
    cbins = np.unique(np.quantile(d["curv"][isg], np.linspace(0, 1, 11))); cbins[-1] += 1e-6
    sbins = np.array([0, 4, 6, 8, 11, 15, 22, 39, 41])
    out = {}
    for name, vals, bins in (("curvature", d["curv"], cbins),
                             ("spacing", d["sp"], sbins)):
        print(f"\n=== per-bin AUC by {name} (A=059_redo, B=surface_recto) ===")
        tab = []
        for b in range(len(bins)-1):
            m = (vals >= bins[b]) & (vals < bins[b+1])
            pos = m & isg; neg = m & ~isg
            aA = auc(d["pA"][pos], d["pA"][neg])
            aB = auc(d["pB"][pos], d["pB"][neg])
            tab.append([float(bins[b]), float(bins[b+1]), int(pos.sum()), int(neg.sum()),
                        aA, aB])
            print(f"  [{bins[b]:.3f},{bins[b+1]:.3f})  nP={int(pos.sum()):>8,} "
                  f"nN={int(neg.sum()):>8,}  AUC_A={aA:.3f}  AUC_B={aB:.3f}")
        out[name] = tab
    json.dump(out, open(os.path.join(DIAG, "auc.json"), "w"), indent=1)
    print("\nDIAG4 COMPLETE ->", DIAG, flush=True)

if __name__ == "__main__":
    main()
