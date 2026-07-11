"""Threshold-fairness sweep: baseline vs ft_ctrl at th in {0.3,0.4,0.5,0.6},
same 40 patches, same 3 components. Compare each arm at its own best blend."""
import os, glob, sys, json, time, csv
import numpy as np, tifffile, torch
from scipy import ndimage as ndi

sys.path.insert(0, "/mnt/vesuvius")
sys.path.insert(0, "/mnt/vesuvius/Betti-Matching-3D/build")
from loader059 import load_059, predict
from topo_check import cc, topo_matched_frac
from surface_distance import compute_surface_distances, compute_surface_dice_at_tolerance
from skimage.metrics import variation_of_information

OUT = "/mnt/vesuvius/surf191_rand"
FT = "/mnt/vesuvius/experiments/FT191/ckpt_ft_ctrl.pth"
THS = (0.3, 0.4, 0.5, 0.6)
N = 40
W = (.30, .35, .35)  # topo, sdice, voi

def comps(gt, pred):
    sd = compute_surface_distances(gt, pred, spacing_mm=(1., 1., 1.))
    sdice = compute_surface_dice_at_tolerance(sd, 2.0)
    vu, vo = variation_of_information(cc(gt), cc(pred))
    return topo_matched_frac(gt, pred), sdice, 1.0 / (1.0 + vu + vo)

def main():
    netA, norm, props = load_059()
    netF, _, _ = load_059()
    ck = torch.load(FT, map_location="cpu", weights_only=False)
    netF.load_state_dict(ck["model"], strict=True); netF.cuda().eval()
    imgs = sorted(glob.glob(os.path.join(OUT, "imagesTr", "*_0000.tif")))
    rng = np.random.default_rng(7)
    picks = sorted(rng.choice(len(imgs), size=min(N, len(imgs)), replace=False))
    acc = {}  # (arm, th) -> list of (topo, sdice, voi)
    for j, i in enumerate(picks):
        f = imgs[i]
        lf = os.path.join(OUT, "labelsTr", os.path.basename(f).replace("_0000.tif", ".tif"))
        if not os.path.exists(lf): continue
        img = tifffile.imread(f); gt = tifffile.imread(lf) > 0
        if gt.sum() < 5000: continue
        t0 = time.time()
        for arm, net in (("base", netA), ("ft", netF)):
            prob = predict(net, img, 192, norm, props)
            for th in THS:
                acc.setdefault((arm, th), []).append(comps(gt, prob > th))
        print(f"[{j+1}/{len(picks)}] {os.path.basename(f)[:22]:22s} ({time.time()-t0:.0f}s)",
              flush=True)
    print("\n=== SWEEP MEANS: blend = .30*topo + .35*sdice + .35*voi ===")
    print("arm   th    topo   sdice   voi    BLEND")
    best = {}
    out = {}
    for (arm, th), v in sorted(acc.items()):
        a = np.array(v).mean(0)
        blend = W[0]*a[0] + W[1]*a[1] + W[2]*a[2]
        out[f"{arm}@{th}"] = [float(x) for x in a] + [float(blend)]
        print(f"{arm:5s} {th:.1f}  {a[0]:.3f}  {a[1]:.3f}  {a[2]:.3f}  {blend:.3f}")
        if arm not in best or blend > best[arm][1]:
            best[arm] = (th, blend)
    print(f"\nbest base: th={best['base'][0]} blend={best['base'][1]:.3f}")
    print(f"best ft  : th={best['ft'][0]} blend={best['ft'][1]:.3f}")
    print(f"delta at each arm's best: {best['ft'][1]-best['base'][1]:+.3f}")
    json.dump(out, open(os.path.join(OUT, "eval_ft", "topo_sweep.json"), "w"), indent=1)
    print("TOPO SWEEP COMPLETE", flush=True)

if __name__ == "__main__":
    main()
