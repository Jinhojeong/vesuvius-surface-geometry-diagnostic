"""3-arm evaluation: baseline vs ft-ctrl vs ft-ign, per-geometry-bin AUC
on the held-out 200 eval patches (same protocol as diag4)."""
import os, glob, sys, json, time
import numpy as np, tifffile, torch
from scipy import ndimage as ndi

OUT = "/mnt/vesuvius/surf191_rand"
DIAG = os.path.join(OUT, "eval_ft"); os.makedirs(DIAG, exist_ok=True)
FTDIR = "/mnt/vesuvius/experiments/FT191"
K = 30000

sys.path.insert(0, "/mnt/vesuvius")
from diag2_191 import geometry, curvature_at, spacing_at
from loader059 import load_059, predict
from diag4_auc import auc

def load_arm(ckpt_path):
    net, norm, props = load_059()   # baseline arch + weights
    if ckpt_path:
        ck = torch.load(ckpt_path, map_location="cpu", weights_only=False)
        net.load_state_dict(ck["model"], strict=True)
        net.cuda().eval()
    return net, norm, props

def main():
    arms = [("baseline", None),
            ("ft_ctrl", os.path.join(FTDIR, "ckpt_ft_ctrl.pth")),
            ("ft_ign", os.path.join(FTDIR, "ckpt_ft_ign.pth"))]
    nets = {}
    for name, p in arms:
        nets[name], norm, props = load_arm(p)
        print("loaded", name, flush=True)
    imgs = sorted(glob.glob(os.path.join(OUT, "imagesTr", "*_0000.tif")))
    rng = np.random.default_rng(0)
    recs = {"curv": [], "sp": [], "is_gt": []}
    for name, _ in arms: recs["p_" + name] = []
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
        probs = {name: predict(nets[name], img, 192, norm, props) for name, _ in arms}
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
            recs["curv"].append(curvature_at(apts, mom))
            recs["sp"].append(spacing_at(apts, gt, n))
            recs["is_gt"].append(np.full(len(pts), isg, np.uint8))
            for name, _ in arms:
                recs["p_" + name].append(probs[name][pts[:, 0], pts[:, 1], pts[:, 2]])
        print(f"[{i+1}/{len(imgs)}] {os.path.basename(f)[:26]:26s} ({time.time()-t0:.0f}s)",
              flush=True)
    d = {k: np.concatenate(v) for k, v in recs.items()}
    np.savez_compressed(os.path.join(DIAG, "points3.npz"), **d)
    isg = d["is_gt"] == 1
    cbins = np.unique(np.quantile(d["curv"][isg], np.linspace(0, 1, 11))); cbins[-1] += 1e-6
    sbins = np.array([0, 4, 6, 8, 11, 15, 22, 39, 41])
    out = {}
    for gname, vals, bins in (("curvature", d["curv"], cbins),
                              ("spacing", d["sp"], sbins)):
        print(f"\n=== per-bin AUC by {gname} ===")
        hdr = "  bin              " + "".join(f"{nm:>10s}" for nm, _ in arms)
        print(hdr); tab = []
        for b in range(len(bins)-1):
            m = (vals >= bins[b]) & (vals < bins[b+1])
            pos = m & isg; neg = m & ~isg
            row = [float(bins[b]), float(bins[b+1]), int(pos.sum()), int(neg.sum())]
            cells = ""
            for nm, _ in arms:
                a = auc(d["p_"+nm][pos], d["p_"+nm][neg])
                row.append(a); cells += f"{a:>10.3f}"
            tab.append(row)
            print(f"  [{bins[b]:6.3f},{bins[b+1]:6.3f}) {cells}")
        out[gname] = {"bins": tab, "arms": [nm for nm, _ in arms]}
    json.dump(out, open(os.path.join(DIAG, "auc3.json"), "w"), indent=1)
    print("\nEVAL_FT COMPLETE ->", DIAG, flush=True)

if __name__ == "__main__":
    main()
