"""CT-brightness oracle: are pred+ voxels on brighter (sheet) CT than GT+ voxels?

For each patch: percentile-normalize CT, then compare mean CT percentile at
  - GT-positive voxels
  - pred-positive voxels (th=0.4)
  - pred-only (pred & ~GT): the voxels dice punishes as "false positives"
  - GT-only (GT & ~pred): the voxels counted as "misses"
If pred-only sits on bright papyrus while GT-only doesn't, labels are the weak
link (unlabeled real sheets), not the model. Correlate the gap with patch dice.
"""
import os, glob, sys, time
import numpy as np, tifffile
sys.path.insert(0, "/mnt/vesuvius")
from loader059 import load_059, predict

OUT = "/mnt/vesuvius/surf191_rand"
net, norm, props = load_059()
imgs = sorted(glob.glob(os.path.join(OUT, "imagesTr", "*_0000.tif")))
rows = []
for i, f in enumerate(imgs):
    lf = os.path.join(OUT, "labelsTr", os.path.basename(f).replace("_0000.tif", ".tif"))
    if not os.path.exists(lf): continue
    img = tifffile.imread(f); gt = tifffile.imread(lf) > 0
    if img.ndim != 3 or gt.shape != img.shape or gt.sum() < 5000: continue
    prob = predict(net, img, 192, norm, props)
    pred = prob > 0.4
    # CT percentile map (cheap: rank against sampled quantiles)
    qs = np.quantile(img, np.linspace(0, 1, 101))
    pct = np.searchsorted(qs, img).astype(np.float32)  # 0..100
    inter = (pred & gt).sum()
    dice = 2*inter/(pred.sum()+gt.sum()+1e-6)
    m_gt = pct[gt].mean()
    m_pr = pct[pred].mean() if pred.any() else np.nan
    po = pred & ~gt; go = gt & ~pred
    m_po = pct[po].mean() if po.any() else np.nan   # "false positives"
    m_go = pct[go].mean() if go.any() else np.nan   # "misses"
    m_bg = pct[~gt & ~pred].mean()
    rows.append((os.path.basename(f), float(dice), float(m_gt), float(m_pr),
                 float(m_po), float(m_go), float(m_bg)))
    if (i+1) % 25 == 0:
        print(f"[{i+1}/{len(imgs)}]", flush=True)
r = np.array([x[1:] for x in rows], np.float64)
dice, m_gt, m_pr, m_po, m_go, m_bg = r.T
print(f"\npatches: {len(rows)}")
print(f"mean CT percentile: GT+ {np.nanmean(m_gt):.1f} | pred+ {np.nanmean(m_pr):.1f} "
      f"| pred-only(FP) {np.nanmean(m_po):.1f} | GT-only(miss) {np.nanmean(m_go):.1f} "
      f"| background {np.nanmean(m_bg):.1f}")
lo = dice < 0.5; hi = dice >= 0.7
for name, m in (("dice<0.5", lo), ("dice>=0.7", hi)):
    print(f"{name} (n={m.sum()}): FP-voxel CT pct {np.nanmean(m_po[m]):.1f} vs "
          f"miss-voxel {np.nanmean(m_go[m]):.1f} vs GT+ {np.nanmean(m_gt[m]):.1f} "
          f"vs bg {np.nanmean(m_bg[m]):.1f}")
import csv
with open(f"{OUT}/diag2/oracle.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["patch","dice","ctpct_gt","ctpct_pred","ctpct_predonly","ctpct_gtonly","ctpct_bg"])
    w.writerows(rows)
print("ORACLE DONE", flush=True)
