"""Idea-4 v0: valley-carving sheet peel for fused compressed regions.

Where the predicted surface is thicker than one sheet, carve a cut along the
prob VALLEY (local min along the normal) — approximates the directional
min-cut without a graph solver. Score with the OFFICIAL topometrics blend on
the standard 40 patches: base@0.4, base+peel, full@0.5, full+peel.
Pre-registered: peel must not regress the blend and should raise topo.
"""
import os, glob, sys, json, time
import numpy as np, tifffile, torch
from scipy import ndimage as ndi

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from loader059 import load_059, predict
from official_score import official as _official_direct
import subprocess, tempfile
def official(gt, pred):
    # Betti-matching C++ can abort the whole process on some patches; isolate it
    # in a subprocess so one bad patch does not kill the run.
    with tempfile.TemporaryDirectory() as td:
        gp = f"{td}/gt.npy"; pp = f"{td}/pred.npy"
        np.save(gp, gt); np.save(pp, pred)
        try:
            r = subprocess.run([sys.executable, os.path.join(HERE, "official_one.py"), gp, pp],
                               capture_output=True, text=True, timeout=180)
            return [float(x) for x in r.stdout.split()]
        except Exception:
            return [float("nan")]*4

OUT = "/mnt/vesuvius/surf191_rand"
FULL = "/mnt/vesuvius/experiments/FT191/ckpt_ft_full.pth"
THICK_HALF = 2.6      # d_in above this = thicker than one sheet core
VALLEY_EPS = 0.0

def peel(pred_bin, prob):
    d_in = ndi.distance_transform_edt(pred_bin).astype(np.float32)
    d_out = ndi.distance_transform_edt(~pred_bin).astype(np.float32)
    sdf = ndi.gaussian_filter(d_out - d_in, 2.0)
    del d_out
    n = np.empty((3,) + sdf.shape, np.float32)
    for ax in range(3):
        n[ax] = np.gradient(sdf, axis=ax).astype(np.float32)
    del sdf
    mag = np.sqrt((n ** 2).sum(0)) + 1e-6
    n /= mag[None]
    pts = np.argwhere(pred_bin & (d_in >= 1.5))
    if len(pts) == 0: return pred_bin
    nv = n[:, pts[:, 0], pts[:, 1], pts[:, 2]].T
    p0 = prob[pts[:, 0], pts[:, 1], pts[:, 2]]
    pp = ndi.map_coordinates(prob, (pts + 1.5 * nv).T, order=1, mode="nearest")
    pm = ndi.map_coordinates(prob, (pts - 1.5 * nv).T, order=1, mode="nearest")
    dv = d_in[pts[:, 0], pts[:, 1], pts[:, 2]]
    cut = (p0 <= pp + VALLEY_EPS) & (p0 <= pm + VALLEY_EPS) & (dv >= THICK_HALF)
    out = pred_bin.copy()
    cp = pts[cut]
    cn = nv[cut]
    # v1: 3-vox cut along the normal so it actually disconnects (26-conn);
    # v0's 1-vox carve only made tunnels (topo crashed, VOI unchanged)
    for t in (-1.0, 0.0, 1.0):
        q = np.round(cp + t * cn).astype(int)
        q = q[((q >= 0) & (q < np.array(out.shape))).all(1)]
        out[q[:, 0], q[:, 1], q[:, 2]] = False
    return out

def main():
    netB, norm, props = load_059()
    netF, _, _ = load_059()
    ck = torch.load(FULL, map_location="cpu", weights_only=False)
    netF.load_state_dict(ck["model"], strict=True); netF.cuda().eval()
    imgs = sorted(glob.glob(os.path.join(OUT, "imagesTr", "*_0000.tif")))
    rng7 = np.random.default_rng(7)
    picks = sorted(rng7.choice(len(imgs), size=40, replace=False))
    acc = {k: [] for k in ("base", "base_peel", "full", "full_peel")}
    for j, i in enumerate(picks):
        f = imgs[i]
        lf = os.path.join(OUT, "labelsTr", os.path.basename(f).replace("_0000.tif", ".tif"))
        if not os.path.exists(lf): continue
        img = tifffile.imread(f); gt = tifffile.imread(lf) > 0
        if gt.sum() < 5000: continue
        t0 = time.time()
        pB = predict(netB, img, 192, norm, props)
        pF = predict(netF, img, 192, norm, props)
        cfg = {"base": pB > 0.4, "full": pF > 0.5}
        cfg["base_peel"] = peel(cfg["base"], pB)
        cfg["full_peel"] = peel(cfg["full"], pF)
        for k, pb in cfg.items():
            acc[k].append(official(gt, pb))
        cutB = int(cfg["base"].sum() - cfg["base_peel"].sum())
        print(f"[{j+1}/40] cutB={cutB:,} ({time.time()-t0:.0f}s)", flush=True)
    print("\n=== IDEA-4 PEEL: OFFICIAL (means over 40) ===")
    print("config      score   topo   sdice   voi")
    out = {}
    for k, v in acc.items():
        a = np.nanmean(np.array(v), 0)
        out[k] = [float(x) for x in a]
        print(f"{k:10s}  {a[0]:.3f}  {a[1]:.3f}  {a[2]:.3f}  {a[3]:.3f}")
    json.dump(out, open(os.path.join(OUT, "eval_ft", "peel_v1_official.json"), "w"), indent=1)
    print("PEEL COMPLETE", flush=True)

if __name__ == "__main__":
    main()
