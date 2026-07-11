"""FINAL: official topometrics leaderboard score, baseline vs ft_ctrl.
Configs: base@0.4 (its best), ft@0.6 (its best), ft@0.4 (fixed-th reference).
Same 40 patches as topo_check/sweep (rng 7)."""
import os, glob, sys, json, time, csv
import numpy as np, tifffile, torch

sys.path.insert(0, "/mnt/vesuvius")
from loader059 import load_059, predict
import topometrics.leaderboard as lb

OUT = "/mnt/vesuvius/surf191_rand"
FT = "/mnt/vesuvius/experiments/FT191/ckpt_ft_ctrl.pth"
N = 40
CONFIGS = [("base", 0.4), ("ft", 0.6), ("ft", 0.4)]

def official(gt, pred):
    r = lb.compute_leaderboard_score(
        predictions=pred.astype(np.uint8), labels=gt.astype(np.uint8),
        dims=(0, 1, 2), spacing=(1., 1., 1.), surface_tolerance=2.0,
        voi_connectivity=26, voi_transform="one_over_one_plus", voi_alpha=0.3,
        combine_weights=(0.3, 0.35, 0.35), fg_threshold=None,
        ignore_label=2, ignore_mask=None)
    return (float(r.score), float(r.topo.toposcore),
            float(r.surface_dice), float(r.voi.voi_score))

def main():
    netA, norm, props = load_059()
    netF, _, _ = load_059()
    ck = torch.load(FT, map_location="cpu", weights_only=False)
    netF.load_state_dict(ck["model"], strict=True); netF.cuda().eval()
    nets = {"base": netA, "ft": netF}
    imgs = sorted(glob.glob(os.path.join(OUT, "imagesTr", "*_0000.tif")))
    rng = np.random.default_rng(7)
    picks = sorted(rng.choice(len(imgs), size=min(N, len(imgs)), replace=False))
    acc = {c: [] for c in CONFIGS}
    rows = []
    for j, i in enumerate(picks):
        f = imgs[i]
        lf = os.path.join(OUT, "labelsTr", os.path.basename(f).replace("_0000.tif", ".tif"))
        if not os.path.exists(lf): continue
        img = tifffile.imread(f); gt = tifffile.imread(lf) > 0
        if gt.sum() < 5000: continue
        t0 = time.time()
        probs = {arm: predict(nets[arm], img, 192, norm, props)
                 for arm in ("base", "ft")}
        line = [os.path.basename(f)]
        for c in CONFIGS:
            arm, th = c
            s = official(gt, probs[arm] > th)
            acc[c].append(s); line += list(s)
        rows.append(line)
        b = acc[("base", 0.4)][-1]; ftb = acc[("ft", 0.6)][-1]
        print(f"[{j+1}/{len(picks)}] {os.path.basename(f)[:22]:22s} "
              f"base {b[0]:.3f} -> ft@.6 {ftb[0]:.3f}  ({time.time()-t0:.0f}s)",
              flush=True)
    print("\n=== OFFICIAL LEADERBOARD SCORE (means over patches) ===")
    print("config     score   topo   sdice   voi")
    out = {}
    for c in CONFIGS:
        a = np.array(acc[c]).mean(0)
        out[f"{c[0]}@{c[1]}"] = [float(x) for x in a]
        print(f"{c[0]}@{c[1]:<4}  {a[0]:.3f}  {a[1]:.3f}  {a[2]:.3f}  {a[3]:.3f}")
    d = out["ft@0.6"][0] - out["base@0.4"][0]
    print(f"\nDELTA (ft@best vs base@best): {d:+.4f}")
    json.dump(out, open(os.path.join(OUT, "eval_ft", "official_score.json"), "w"), indent=1)
    with open(os.path.join(OUT, "eval_ft", "official_per_patch.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        hdr = ["patch"]
        for c in CONFIGS:
            hdr += [f"{c[0]}@{c[1]}_{k}" for k in ("score", "topo", "sdice", "voi")]
        w.writerow(hdr); w.writerows(rows)
    print("OFFICIAL SCORE COMPLETE", flush=True)

if __name__ == "__main__":
    main()
