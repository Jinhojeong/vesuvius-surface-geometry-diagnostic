"""Trade-off check: did CE+dice fine-tuning buy AUC at the cost of topology?

Scores baseline vs ft_ctrl on official-metric COMPONENTS (approximation of the
Kaggle blend; exact wrapper lives in the Kaggle dataset):
  - SurfaceDice @ tolerance 2.0 (surface_distance lib)  [official weight .35]
  - VOI score = 1/(1+VOI), connectivity-26 CCs           [official weight .35]
  - Betti-matching topo error (matched fraction)         [official weight .30]
40 eval patches, threshold 0.4 both arms.
"""
import os, glob, sys, json, time
import numpy as np, tifffile, torch
from scipy import ndimage as ndi

sys.path.insert(0, "/mnt/vesuvius")
sys.path.insert(0, "/mnt/vesuvius/Betti-Matching-3D/build")
from loader059 import load_059, predict
import betti_matching
from surface_distance import compute_surface_distances, compute_surface_dice_at_tolerance
from skimage.metrics import variation_of_information

OUT = "/mnt/vesuvius/surf191_rand"
FT = "/mnt/vesuvius/experiments/FT191/ckpt_ft_ctrl.pth"
N = 40; TH = 0.4
S26 = np.ones((3, 3, 3), bool)

def cc(vol):
    lab, _ = ndi.label(vol, structure=S26)
    return lab

def topo_matched_frac(gt, pred):
    """Betti-matching: fraction of GT topological features matched in pred.
    Uses distance-transform filtrations of the binary masks (downsampled 2x
    for tractability)."""
    g = gt[::2, ::2, ::2]; p = pred[::2, ::2, ::2]
    fg = ndi.distance_transform_edt(~g) - ndi.distance_transform_edt(g)
    fp = ndi.distance_transform_edt(~p) - ndi.distance_transform_edt(p)
    r = betti_matching.compute_matching(fg.astype(np.float64), fp.astype(np.float64))
    matched = float(np.sum(r.num_matched))
    un1 = float(np.sum(r.num_unmatched_input1))
    un2 = float(np.sum(r.num_unmatched_input2))
    return matched / max(matched + 0.5 * (un1 + un2), 1)

def main():
    netA, norm, props = load_059()
    netF, _, _ = load_059()
    ck = torch.load(FT, map_location="cpu", weights_only=False)
    netF.load_state_dict(ck["model"], strict=True); netF.cuda().eval()
    imgs = sorted(glob.glob(os.path.join(OUT, "imagesTr", "*_0000.tif")))
    rng = np.random.default_rng(7)
    picks = rng.choice(len(imgs), size=min(N, len(imgs)), replace=False)
    rows = []
    for j, i in enumerate(sorted(picks)):
        f = imgs[i]
        lf = os.path.join(OUT, "labelsTr", os.path.basename(f).replace("_0000.tif", ".tif"))
        if not os.path.exists(lf): continue
        img = tifffile.imread(f); gt = tifffile.imread(lf) > 0
        if gt.sum() < 5000: continue
        t0 = time.time()
        res = {}
        for name, net in (("base", netA), ("ft", netF)):
            pred = predict(net, img, 192, norm, props) > TH
            sd = compute_surface_distances(gt, pred, spacing_mm=(1., 1., 1.))
            sdice = compute_surface_dice_at_tolerance(sd, 2.0)
            v_under, v_over = variation_of_information(cc(gt), cc(pred))
            voi_score = 1.0 / (1.0 + v_under + v_over)
            topo = topo_matched_frac(gt, pred)
            res[name] = (sdice, voi_score, topo)
        rows.append((os.path.basename(f), *res["base"], *res["ft"]))
        b, ftv = res["base"], res["ft"]
        print(f"[{j+1}/{len(picks)}] {os.path.basename(f)[:24]:24s} "
              f"sdice {b[0]:.3f}->{ftv[0]:.3f}  voi {b[1]:.3f}->{ftv[1]:.3f}  "
              f"topo {b[2]:.3f}->{ftv[2]:.3f}  ({time.time()-t0:.0f}s)", flush=True)
    a = np.array([r[1:] for r in rows], np.float64)
    names = ["sdice", "voi", "topo"]
    print("\n=== MEANS (baseline -> ft_ctrl) ===")
    blend_b = blend_f = 0.0
    W = {"sdice": .35, "voi": .35, "topo": .30}
    for k, nme in enumerate(names):
        mb, mf = a[:, k].mean(), a[:, k + 3].mean()
        blend_b += W[nme] * mb; blend_f += W[nme] * mf
        print(f"  {nme:6s}: {mb:.3f} -> {mf:.3f}  (delta {mf-mb:+.3f})")
    print(f"  BLEND (.3 topo/.35 sdice/.35 voi): {blend_b:.3f} -> {blend_f:.3f} "
          f"(delta {blend_f-blend_b:+.3f})")
    import csv
    with open(os.path.join(OUT, "eval_ft", "topo_check.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["patch", "b_sdice", "b_voi", "b_topo", "f_sdice", "f_voi", "f_topo"])
        w.writerows(rows)
    print("\nTOPO CHECK COMPLETE", flush=True)

if __name__ == "__main__":
    main()
