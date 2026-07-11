"""Generate ignore-masked labels for the fine-tune set (arm: ft-ignore).

For each patch in surf191_ft: run model A, detect candidate unlabeled sheets
(same filters as qc192_labels), write labelsIg/{name}.tif with:
  0 = background, 1 = surface (GT), 2 = ignore (candidate unlabeled sheet).
"""
import os, glob, sys, time
import numpy as np, tifffile
from scipy import ndimage as ndi

FT = "/mnt/vesuvius/surf191_ft"
OUTL = os.path.join(FT, "labelsIg"); os.makedirs(OUTL, exist_ok=True)

sys.path.insert(0, "/mnt/vesuvius")
from loader059 import load_059, predict
from qc192_labels import sheetlike, TH_CONF, MIN_COMP, CT_PCT_MIN

def main():
    imgs = sorted(glob.glob(os.path.join(FT, "imagesTr", "*_0000.tif")))
    net, norm, props = load_059()
    tot_ign = 0; tot_gt = 0
    for i, f in enumerate(imgs):
        base = os.path.basename(f)
        lf = os.path.join(FT, "labelsTr", base.replace("_0000.tif", ".tif"))
        of = os.path.join(OUTL, base.replace("_0000.tif", ".tif"))
        if os.path.exists(of): continue
        try:
            img = tifffile.imread(f); lbl = tifffile.imread(lf)
        except Exception as e:
            print("skip", base, e); continue
        if img.ndim != 3 or lbl.shape != img.shape: continue
        gt = lbl > 0
        t0 = time.time()
        prob = predict(net, img, 192, norm, props)
        cand = (prob > TH_CONF) & (ndi.distance_transform_edt(~gt) > 3)
        qs = np.quantile(img, np.linspace(0, 1, 101))
        pct = np.searchsorted(qs, img).astype(np.float32)
        la, nc = ndi.label(cand)
        ign = np.zeros_like(gt)
        if nc:
            sizes = np.bincount(la.ravel())
            for cid in np.nonzero(sizes >= MIN_COMP)[0]:
                if cid == 0: continue
                m = la == cid
                pts = np.argwhere(m)
                ok, _ = sheetlike(pts[:: max(1, len(pts)//4000)])
                if ok and pct[m].mean() >= CT_PCT_MIN:
                    ign |= m
        out = gt.astype(np.uint8)
        out[ign] = 2
        tifffile.imwrite(of, out, compression="zlib")
        tot_ign += int(ign.sum()); tot_gt += int(gt.sum())
        print(f"[{i+1}/{len(imgs)}] {base[:30]:30s} ignore={ign.sum():,} "
              f"({time.time()-t0:.0f}s)", flush=True)
    print(f"\nTOTAL ignore/gt = {tot_ign:,}/{tot_gt:,} = {tot_ign/max(tot_gt,1):.3f}")
    print("IGNORE MASKS COMPLETE", flush=True)

if __name__ == "__main__":
    main()
