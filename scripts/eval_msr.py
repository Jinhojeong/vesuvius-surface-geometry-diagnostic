"""Evaluate ft_msr: (1) per-bin AUC at points aligned with eval_ft's points3.npz,
(2) official leaderboard score at th 0.4/0.5/0.6 on the standard 40 patches."""
import os, glob, sys, json, time
import numpy as np, tifffile, torch
from scipy import ndimage as ndi

OUT = "/mnt/vesuvius/surf191_rand"
CKPT = "/mnt/vesuvius/experiments/FT191/ckpt_ft_msr.pth"
K = 30000

sys.path.insert(0, "/mnt/vesuvius")
from diag2_191 import geometry, curvature_at, spacing_at
from loader059 import load_059, predict
from diag4_auc import auc
from official_score import official

def main():
    net, norm, props = load_059()
    ck = torch.load(CKPT, map_location="cpu", weights_only=False)
    net.load_state_dict(ck["model"], strict=True); net.cuda().eval()
    imgs = sorted(glob.glob(os.path.join(OUT, "imagesTr", "*_0000.tif")))

    # ---- part 1: aligned AUC strata (same rng protocol as eval_ft) ----
    rng = np.random.default_rng(0)
    allp = []; allgt = []; allsp = []
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
        prob = predict(net, img, 192, norm, props)
        n, mom = geometry(gt)
        idx = np.argwhere(gt)
        take = rng.choice(len(idx), size=min(K, len(idx)), replace=False)
        pts_gt = idx[take]
        edt, inds = ndi.distance_transform_edt(~gt, return_indices=True)
        shell = (edt >= 2) & (edt <= 10)
        idxb = np.argwhere(shell)
        takeb = rng.choice(len(idxb), size=min(K, len(idxb)), replace=False)
        pts_bg = idxb[takeb]
        anchor = np.stack([inds[c][pts_bg[:, 0], pts_bg[:, 1], pts_bg[:, 2]]
                           for c in range(3)], 1)
        for pts, apts, isg in ((pts_gt, pts_gt, 1), (pts_bg, anchor, 0)):
            allsp.append(spacing_at(apts, gt, n))
            allgt.append(np.full(len(pts), isg, np.uint8))
            allp.append(prob[pts[:, 0], pts[:, 1], pts[:, 2]])
        print(f"[auc {i+1}/{len(imgs)}] ({time.time()-t0:.0f}s)", flush=True)
    p = np.concatenate(allp); isg = np.concatenate(allgt) == 1
    sp = np.concatenate(allsp)
    d3 = np.load(os.path.join(OUT, "eval_ft", "points3.npz"))
    assert len(d3["is_gt"]) == len(isg) and (d3["is_gt"] == isg.astype(np.uint8)).all(), \
        "alignment with eval_ft points failed"
    sbins = np.array([0, 4, 6, 8, 11, 15, 22, 39, 41])
    print("\n=== per-spacing-bin AUC: baseline / ft_ctrl / ft_msr ===")
    tab = []
    for b in range(len(sbins)-1):
        m = (sp >= sbins[b]) & (sp < sbins[b+1])
        row = [auc(d3["p_baseline"][m & isg], d3["p_baseline"][m & ~isg]),
               auc(d3["p_ft_ctrl"][m & isg], d3["p_ft_ctrl"][m & ~isg]),
               auc(p[m & isg], p[m & ~isg])]
        tab.append([float(sbins[b]), float(sbins[b+1])] + [float(x) for x in row])
        print(f"  [{sbins[b]:4.0f},{sbins[b+1]:4.0f})  base={row[0]:.3f}  "
              f"ctrl={row[1]:.3f}  msr={row[2]:.3f}")
    json.dump(tab, open(os.path.join(OUT, "eval_ft", "msr_auc.json"), "w"), indent=1)

    # ---- part 2: official score at 3 thresholds, standard 40 patches ----
    rng7 = np.random.default_rng(7)
    picks = sorted(rng7.choice(len(imgs), size=min(40, len(imgs)), replace=False))
    acc = {th: [] for th in (0.4, 0.5, 0.6)}
    for j, i in enumerate(picks):
        f = imgs[i]
        lf = os.path.join(OUT, "labelsTr", os.path.basename(f).replace("_0000.tif", ".tif"))
        if not os.path.exists(lf): continue
        img = tifffile.imread(f); gt = tifffile.imread(lf) > 0
        if gt.sum() < 5000: continue
        t0 = time.time()
        prob = predict(net, img, 192, norm, props)
        for th in acc:
            acc[th].append(official(gt, prob > th))
        print(f"[off {j+1}/{len(picks)}] ({time.time()-t0:.0f}s)", flush=True)
    print("\n=== ft_msr OFFICIAL (means) — compare base@0.4: 0.578 / ft_ctrl@0.6: 0.567 ===")
    out = {}
    for th, v in acc.items():
        a = np.array(v).mean(0)
        out[f"msr@{th}"] = [float(x) for x in a]
        print(f"msr@{th}: score={a[0]:.3f} topo={a[1]:.3f} sdice={a[2]:.3f} voi={a[3]:.3f}")
    json.dump(out, open(os.path.join(OUT, "eval_ft", "msr_official.json"), "w"), indent=1)
    print("EVAL_MSR COMPLETE", flush=True)

if __name__ == "__main__":
    main()
