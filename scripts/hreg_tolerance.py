"""H-REG: is the compressed-bin plateau (~0.65) a label-REGISTRATION artifact?

Score each GT point by the MAX prob within +/-2 vox along the sheet normal
(tolerant-positive), bg points stay exact. If labels are locally mis-registered
in tight regions, tolerant AUC should jump there; if the model genuinely can't
separate sheets, tolerance adds little (bg stays low either way).
Models: ft_full and baseline. Same sampling protocol as diag4 (self-contained).
"""
import os, glob, sys, json, time
import numpy as np, tifffile, torch
from scipy import ndimage as ndi

OUT = "/mnt/vesuvius/surf191_rand"
FULL = "/mnt/vesuvius/experiments/FT191/ckpt_ft_full.pth"
K = 20000
OFFS = (-2, -1, 0, 1, 2)

sys.path.insert(0, "/mnt/vesuvius")
from diag2_191 import geometry, spacing_at
from loader059 import load_059, predict
from diag4_auc import auc

def probs_at(prob, pts):
    return prob[pts[:, 0], pts[:, 1], pts[:, 2]]

def tolerant_probs(prob, pts, n):
    """max prob over +/-2 vox along the local normal (trilinear)."""
    pz, py, px = pts[:, 0], pts[:, 1], pts[:, 2]
    nz = n[0][pz, py, px]; ny = n[1][pz, py, px]; nx = n[2][pz, py, px]
    best = np.full(len(pts), -1.0, np.float32)
    for t in OFFS:
        cz = pz + t * nz; cy = py + t * ny; cx = px + t * nx
        v = ndi.map_coordinates(prob, np.stack([cz, cy, cx]), order=1, mode="nearest")
        best = np.maximum(best, v)
    return best

def main():
    netB, norm, props = load_059()
    netF, _, _ = load_059()
    ck = torch.load(FULL, map_location="cpu", weights_only=False)
    netF.load_state_dict(ck["model"], strict=True); netF.cuda().eval()
    imgs = sorted(glob.glob(os.path.join(OUT, "imagesTr", "*_0000.tif")))
    rng = np.random.default_rng(3)
    rec = {k: [] for k in ("sp", "isg", "bx", "btol", "fx", "ftol")}
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
        pB = predict(netB, img, 192, norm, props)
        pF = predict(netF, img, 192, norm, props)
        n, mom = geometry(gt)
        idx = np.argwhere(gt)
        pts_gt = idx[rng.choice(len(idx), size=min(K, len(idx)), replace=False)]
        edt, inds = ndi.distance_transform_edt(~gt, return_indices=True)
        shell = (edt >= 2) & (edt <= 10)
        idxb = np.argwhere(shell)
        pts_bg = idxb[rng.choice(len(idxb), size=min(K, len(idxb)), replace=False)]
        anchor = np.stack([inds[c][pts_bg[:, 0], pts_bg[:, 1], pts_bg[:, 2]]
                           for c in range(3)], 1)
        # GT points: exact + tolerant. BG points: exact only (tolerance would
        # let shell points reach the sheet and poison the comparison).
        rec["sp"].append(spacing_at(pts_gt, gt, n)); rec["isg"].append(np.ones(len(pts_gt), np.uint8))
        rec["bx"].append(probs_at(pB, pts_gt)); rec["btol"].append(tolerant_probs(pB, pts_gt, n))
        rec["fx"].append(probs_at(pF, pts_gt)); rec["ftol"].append(tolerant_probs(pF, pts_gt, n))
        rec["sp"].append(spacing_at(anchor, gt, n)); rec["isg"].append(np.zeros(len(pts_bg), np.uint8))
        for k, p in (("bx", pB), ("btol", pB), ("fx", pF), ("ftol", pF)):
            rec[k].append(probs_at(p, pts_bg))
        print(f"[{i+1}/{len(imgs)}] ({time.time()-t0:.0f}s)", flush=True)
    d = {k: np.concatenate(v) for k, v in rec.items()}
    isg = d["isg"] == 1; sp = d["sp"]
    sbins = np.array([0, 4, 6, 8, 11, 15, 41])
    print("\n=== H-REG: exact vs tolerant(+/-2vox) AUC by spacing ===")
    print("  bin        base_ex  base_tol   full_ex  full_tol")
    tab = []
    for b in range(len(sbins)-1):
        m = (sp >= sbins[b]) & (sp < sbins[b+1])
        row = [auc(d[k][m & isg], d[k][m & ~isg]) for k in ("bx", "btol", "fx", "ftol")]
        tab.append([float(sbins[b]), float(sbins[b+1])] + [float(x) for x in row])
        print(f"  [{sbins[b]:2.0f},{sbins[b+1]:2.0f})    {row[0]:.3f}    {row[1]:.3f}     "
              f"{row[2]:.3f}    {row[3]:.3f}")
    json.dump(tab, open(os.path.join(OUT, "eval_ft", "hreg_tolerance.json"), "w"), indent=1)
    print("\nHREG COMPLETE", flush=True)

if __name__ == "__main__":
    main()
