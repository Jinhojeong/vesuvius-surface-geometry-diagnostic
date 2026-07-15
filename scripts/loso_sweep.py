"""Fair official-metric comparison on the unseen scroll (s4): sweep the
threshold for BOTH arms and compare each at its own best.

The first LOSO eval compared base at its in-distribution best (0.4) against
loso at 0.5/0.6, which is not a fair test: loso's AUC on s4 is much higher, so
its operating point may simply sit elsewhere. This decides whether the -0.064
official-metric drop is a calibration artifact or a real structural regression.
"""
import os, glob, sys, json, time, subprocess, tempfile
import numpy as np, tifffile, torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from loader059 import load_059, predict

RAND = "/mnt/vesuvius/surf191_rand"
FTP = "/mnt/vesuvius/surf191_ft"
CKPT = "/mnt/vesuvius/experiments/FT191/ckpt_ft_loso_s1.pth"
OUTD = "/mnt/vesuvius/experiments/FT191"
NP = 30
THR = [0.3, 0.4, 0.5, 0.6, 0.7]

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

def s4_files():
    a = [(f, os.path.join(RAND, "labelsTr",
          os.path.basename(f).replace("_0000.tif", ".tif")))
         for f in sorted(glob.glob(os.path.join(RAND, "imagesTr", "s4_*_0000.tif")))]
    b = [(f, os.path.join(FTP, "labelsTr",
          os.path.basename(f).replace("_0000.tif", ".tif")))
         for f in sorted(glob.glob(os.path.join(FTP, "imagesTr", "s4_*_0000.tif")))]
    rng = np.random.default_rng(7)
    b = [b[i] for i in rng.permutation(len(b))[:NP]]
    return [(f, l) for f, l in (a + b) if os.path.exists(l)][:NP]

def main():
    base, norm, props = load_059()
    loso, _, _ = load_059()
    ck = torch.load(CKPT, map_location="cpu", weights_only=False)
    assert ck.get("train_scrolls") == ["s1"]
    loso.load_state_dict(ck["model"], strict=True); loso.cuda().eval()
    arms = {"base": base, "loso_s1": loso}
    files = s4_files()
    print(f"threshold sweep on {len(files)} s4 patches (unseen by loso_s1)", flush=True)
    acc = {f"{a}@{t}": [] for a in arms for t in THR}
    for i, (f, lf) in enumerate(files):
        img = tifffile.imread(f); lbl = tifffile.imread(lf)
        if img.ndim != 3 or lbl.shape != img.shape: continue
        gt = lbl > 0
        if gt.sum() < 5000: continue
        t0 = time.time()
        for a, net in arms.items():
            p = predict(net, img, 192, norm, props)
            for t in THR:
                acc[f"{a}@{t}"].append(official(gt, p > t))
        print(f"[{i+1}/{len(files)}] ({time.time()-t0:.0f}s)", flush=True)
    print("\n=== OFFICIAL METRIC, THRESHOLD SWEEP, SCROLL 4 (unseen) ===")
    print(f"{'config':>14s}  {'score':>6s}  {'topo':>6s}  {'sdice':>6s}  {'voi':>6s}")
    out = {}
    for k, v in acc.items():
        if not v: continue
        m = np.nanmean(np.array(v), 0)
        out[k] = [float(x) for x in m]
        print(f"{k:>14s}  {m[0]:6.3f}  {m[1]:6.3f}  {m[2]:6.3f}  {m[3]:6.3f}")
    bb = max((k for k in out if k.startswith("base")), key=lambda k: out[k][0])
    bl = max((k for k in out if k.startswith("loso")), key=lambda k: out[k][0])
    d = out[bl][0] - out[bb][0]
    print(f"\nBEST base: {bb} = {out[bb][0]:.3f}")
    print(f"BEST loso: {bl} = {out[bl][0]:.3f}")
    print(f"FAIR DELTA (each at own best) = {d:+.3f}")
    print("VERDICT:", "no official-metric cost on the unseen scroll" if d > -0.01
          else "REAL regression on the unseen scroll (not a threshold artifact)")
    json.dump({"sweep": out, "best_base": bb, "best_loso": bl, "fair_delta": float(d),
               "n": len(files)}, open(os.path.join(OUTD, "loso_sweep.json"), "w"), indent=1)
    print("LOSO SWEEP COMPLETE", flush=True)

if __name__ == "__main__":
    main()
