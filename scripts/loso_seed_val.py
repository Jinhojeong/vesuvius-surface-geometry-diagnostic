"""Compact per-seed LOSO validation on unseen s4.
Reports <4vox AUC (transfer) and official metric across thresholds (regression),
for base vs one loso checkpoint (CKPT env). Seeds are compared to check the
AUC-up / official-down direction is not a single-seed artifact.
"""
import os, glob, sys, json, time, subprocess, tempfile
import numpy as np, tifffile, torch
from scipy import ndimage as ndi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from loader059 import load_059, predict
from diag2_191 import geometry, spacing_at
from diag4_auc import auc

RAND = "/mnt/vesuvius/surf191_rand"; FTP = "/mnt/vesuvius/surf191_ft"
OUTD = "/mnt/vesuvius/experiments/FT191"
CKPT = os.environ["CKPT"]; TAG = os.environ.get("TAG", "loso")
K = 30000; N_OFF = 30; THR = [0.4, 0.5, 0.6]
SB = np.array([0, 4, 6, 8, 11, 15, 22, 41])

def official(gt, pred):
    with tempfile.TemporaryDirectory() as td:
        gp = f"{td}/g.npy"; pp = f"{td}/p.npy"; np.save(gp, gt); np.save(pp, pred)
        try:
            r = subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), "official_one.py"), gp, pp],
                               capture_output=True, text=True, timeout=180)
            return [float(x) for x in r.stdout.split()]
        except Exception:
            return [float("nan")]*4

def files():
    a = [(f, os.path.join(RAND, "labelsTr", os.path.basename(f).replace("_0000.tif", ".tif")))
         for f in sorted(glob.glob(os.path.join(RAND, "imagesTr", "s4_*_0000.tif")))]
    b = [(f, os.path.join(FTP, "labelsTr", os.path.basename(f).replace("_0000.tif", ".tif")))
         for f in sorted(glob.glob(os.path.join(FTP, "imagesTr", "s4_*_0000.tif")))]
    rng = np.random.default_rng(7); b = [b[i] for i in rng.permutation(len(b))[:120]]
    return [(f, l) for f, l in a + b if os.path.exists(l)]

def main():
    base, norm, props = load_059()
    loso, _, _ = load_059()
    ck = torch.load(CKPT, map_location="cpu", weights_only=False)
    assert ck.get("train_scrolls") == ["s1"]
    loso.load_state_dict(ck["model"], strict=True); loso.cuda().eval()
    seed = ck.get("seed", "?")
    fs = files(); rng = np.random.default_rng(0)
    sp = {"b": [], "l": []}; isg = []
    off = {f"base@{t}": [] for t in THR}; off.update({f"loso@{t}": [] for t in THR})
    for i, (f, lf) in enumerate(fs):
        img = tifffile.imread(f); lbl = tifffile.imread(lf)
        if img.ndim != 3 or lbl.shape != img.shape: continue
        gt = lbl > 0
        if gt.sum() < 5000: continue
        pb = predict(base, img, 192, norm, props); pl = predict(loso, img, 192, norm, props)
        n, mom = geometry(gt)
        idx = np.argwhere(gt); pg = idx[rng.choice(len(idx), min(K, len(idx)), replace=False)]
        edt, inds = ndi.distance_transform_edt(~gt, return_indices=True)
        idxb = np.argwhere((edt >= 2) & (edt <= 10))
        pbg = idxb[rng.choice(len(idxb), min(K, len(idxb)), replace=False)]
        anc = np.stack([inds[c][pbg[:, 0], pbg[:, 1], pbg[:, 2]] for c in range(3)], 1)
        for pts, apts, ig in ((pg, pg, 1), (pbg, anc, 0)):
            sp["b"].append(pb[pts[:, 0], pts[:, 1], pts[:, 2]])
            sp["l"].append(pl[pts[:, 0], pts[:, 1], pts[:, 2]])
            sp.setdefault("spc", []).append(spacing_at(apts, gt, n))
            isg.append(np.full(len(pts), ig, np.uint8))
        if i < N_OFF:
            for t in THR:
                off[f"base@{t}"].append(official(gt, pb > t))
                off[f"loso@{t}"].append(official(gt, pl > t))
        if (i+1) % 20 == 0: print(f"  {i+1}/{len(fs)}", flush=True)
    ig = np.concatenate(isg) == 1
    spc = np.concatenate(sp["spc"]); Pb = np.concatenate(sp["b"]); Pl = np.concatenate(sp["l"])
    res = {"seed": seed, "auc": {}, "official": {}}
    print(f"\n=== SEED {seed}: transfer AUC on unseen s4 ===")
    for j in range(len(SB)-1):
        m = (spc >= SB[j]) & (spc < SB[j+1]); pos = m & ig; neg = m & ~ig
        if pos.sum() < 50 or neg.sum() < 50: continue
        ab = auc(Pb[pos], Pb[neg]); al = auc(Pl[pos], Pl[neg])
        res["auc"][f"{SB[j]}-{SB[j+1]}"] = [round(ab, 3), round(al, 3)]
        print(f"  [{SB[j]:2d},{SB[j+1]:2d})  base {ab:.3f}  loso {al:.3f}  ({al-ab:+.3f})")
    for k, v in off.items():
        if v: res["official"][k] = round(float(np.nanmean(np.array(v), 0)[0]), 3)
    bb = max((k for k in res["official"] if "base" in k), key=lambda k: res["official"][k])
    bl = max((k for k in res["official"] if "loso" in k), key=lambda k: res["official"][k])
    res["fair_delta"] = round(res["official"][bl] - res["official"][bb], 3)
    print(f"SEED {seed}: best base {bb}={res['official'][bb]}  best loso {bl}={res['official'][bl]}"
          f"  FAIR DELTA={res['fair_delta']:+.3f}")
    json.dump(res, open(os.path.join(OUTD, f"loso_val_s{seed}.json"), "w"), indent=1)
    print(f"SEED {seed} VAL COMPLETE", flush=True)

if __name__ == "__main__":
    main()
