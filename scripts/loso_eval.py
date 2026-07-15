"""LOSO external-validity eval: does the compressed-region repair transfer to a
scroll the fine-tune NEVER saw?

Arms: base (released 059) vs loso_s1 (fine-tuned on SCROLL 1 patches only).
Eval: SCROLL 4 patches only -- unseen by loso_s1 by construction.
(ft_full is deliberately excluded: it trained on s4 patches, so it would be
contaminated on this eval set.)
Reports per-spacing-bin AUC and the official topometrics blend (paired).
"""
import os, glob, sys, json, time, subprocess, tempfile
import numpy as np, tifffile, torch
from scipy import ndimage as ndi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from loader059 import load_059, predict
from diag2_191 import geometry, curvature_at, spacing_at
from diag4_auc import auc

RAND = "/mnt/vesuvius/surf191_rand"          # held out from every fine-tune
FTP = "/mnt/vesuvius/surf191_ft"             # original ft pool (s4 unseen by loso_s1)
CKPT = "/mnt/vesuvius/experiments/FT191/ckpt_ft_loso_s1.pth"
OUTD = "/mnt/vesuvius/experiments/FT191"
K = 30000
N_EXTRA_S4 = 150                              # top up from the ft pool
N_OFFICIAL = 40                               # official metric is expensive
SBINS = np.array([0, 4, 6, 8, 11, 15, 22, 39, 41])

def official(gt, pred):
    with tempfile.TemporaryDirectory() as td:
        gp = f"{td}/gt.npy"; pp = f"{td}/pred.npy"
        np.save(gp, gt); np.save(pp, pred)
        try:
            r = subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), "official_one.py"), gp, pp],
                               capture_output=True, text=True, timeout=180)
            return [float(x) for x in r.stdout.split()]
        except Exception:
            return [float("nan")]*4

def s4_list():
    a = [(f, os.path.join(RAND, "labelsTr",
                          os.path.basename(f).replace("_0000.tif", ".tif")))
         for f in sorted(glob.glob(os.path.join(RAND, "imagesTr", "s4_*_0000.tif")))]
    b = [(f, os.path.join(FTP, "labelsTr",
                          os.path.basename(f).replace("_0000.tif", ".tif")))
         for f in sorted(glob.glob(os.path.join(FTP, "imagesTr", "s4_*_0000.tif")))]
    rng = np.random.default_rng(7)
    b = [b[i] for i in rng.permutation(len(b))[:N_EXTRA_S4]]
    out = [(f, l) for f, l in a + b if os.path.exists(l)]
    print(f"s4 eval patches: {len(a)} held-out + {len(b)} from ft-pool "
          f"(all UNSEEN by loso_s1) = {len(out)}", flush=True)
    return out

def main():
    base, norm, props = load_059()
    loso, _, _ = load_059()
    ck = torch.load(CKPT, map_location="cpu", weights_only=False)
    assert ck.get("train_scrolls") == ["s1"], f"ckpt not s1-only: {ck.get('train_scrolls')}"
    loso.load_state_dict(ck["model"], strict=True); loso.cuda().eval()
    arms = {"base": base, "loso_s1": loso}
    files = s4_list()
    rng = np.random.default_rng(0)
    recs = {"sp": [], "is_gt": [], "p_base": [], "p_loso_s1": []}
    offi = {"base@0.4": [], "loso@0.5": [], "loso@0.6": []}
    for i, (f, lf) in enumerate(files):
        try:
            img = tifffile.imread(f); lbl = tifffile.imread(lf)
        except Exception:
            continue
        if img.ndim != 3 or lbl.shape != img.shape: continue
        gt = lbl > 0
        if gt.sum() < 5000: continue
        t0 = time.time()
        probs = {k: predict(v, img, 192, norm, props) for k, v in arms.items()}
        n, mom = geometry(gt)
        idx = np.argwhere(gt)
        pts_gt = idx[rng.choice(len(idx), size=min(K, len(idx)), replace=False)]
        edt, inds = ndi.distance_transform_edt(~gt, return_indices=True)
        shell = (edt >= 2) & (edt <= 10)
        idxb = np.argwhere(shell)
        pts_bg = idxb[rng.choice(len(idxb), size=min(K, len(idxb)), replace=False)]
        anchor = np.stack([inds[c][pts_bg[:, 0], pts_bg[:, 1], pts_bg[:, 2]]
                           for c in range(3)], 1)
        for pts, apts, isg in ((pts_gt, pts_gt, 1), (pts_bg, anchor, 0)):
            recs["sp"].append(spacing_at(apts, gt, n))
            recs["is_gt"].append(np.full(len(pts), isg, np.uint8))
            for k in arms:
                recs["p_"+k].append(probs[k][pts[:, 0], pts[:, 1], pts[:, 2]])
        if i < N_OFFICIAL:
            offi["base@0.4"].append(official(gt, probs["base"] > 0.4))
            offi["loso@0.5"].append(official(gt, probs["loso_s1"] > 0.5))
            offi["loso@0.6"].append(official(gt, probs["loso_s1"] > 0.6))
        print(f"[{i+1}/{len(files)}] {os.path.basename(f)[:24]:24s} ({time.time()-t0:.0f}s)",
              flush=True)
    d = {k: np.concatenate(v) for k, v in recs.items()}
    np.savez_compressed(os.path.join(OUTD, "loso_points.npz"), **d)

    isg = d["is_gt"] == 1
    print("\n=== LOSO: per-SPACING-bin AUC on SCROLL 4 (unseen by the fine-tune) ===")
    print(f"{'spacing':>12s}  {'n_pos':>7s}  {'base':>6s}  {'loso_s1':>7s}  {'delta':>6s}")
    table = []
    for b in range(len(SBINS)-1):
        m = (d["sp"] >= SBINS[b]) & (d["sp"] < SBINS[b+1])
        pos = m & isg; neg = m & ~isg
        if pos.sum() < 50 or neg.sum() < 50:
            continue
        ab = auc(d["p_base"][pos], d["p_base"][neg])
        al = auc(d["p_loso_s1"][pos], d["p_loso_s1"][neg])
        print(f"  [{SBINS[b]:2d},{SBINS[b+1]:2d})  {int(pos.sum()):7d}  "
              f"{ab:6.3f}  {al:7.3f}  {al-ab:+6.3f}")
        table.append({"bin": [int(SBINS[b]), int(SBINS[b+1])], "n_pos": int(pos.sum()),
                      "base": float(ab), "loso_s1": float(al), "delta": float(al-ab)})

    print("\n=== LOSO: official topometrics blend on SCROLL 4 (paired) ===")
    off = {}
    for k, v in offi.items():
        if not v: continue
        a = np.nanmean(np.array(v), 0)
        off[k] = [float(x) for x in a]
        print(f"  {k:10s} score={a[0]:.3f}  topo={a[1]:.3f}  sdice={a[2]:.3f}  voi={a[3]:.3f}"
              f"  (n={len(v)})")
    json.dump({"auc_by_spacing": table, "official": off,
               "n_patches": len(files), "train_scrolls": ["s1"], "eval_scroll": "s4"},
              open(os.path.join(OUTD, "loso_result.json"), "w"), indent=1)
    tight = [t for t in table if t["bin"][0] == 0]
    if tight:
        t = tight[0]
        print(f"\nVERDICT (<4vox, unseen scroll): base {t['base']:.3f} -> "
              f"loso_s1 {t['loso_s1']:.3f} ({t['delta']:+.3f})")
        print("TRANSFERS" if t["delta"] > 0.05 else "DOES NOT TRANSFER (in-distribution only)")
    print("LOSO EVAL COMPLETE", flush=True)

if __name__ == "__main__":
    main()
