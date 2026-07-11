"""H-GRID: causal test of the axis-alignment deficit.

Rotate each patch 30 deg about axis 0 (resampled, same content), run the model
on the rotated volume, and score the SAME physical points in both frames.
If the deficit is grid-interaction, points whose normals move
aligned->oblique should DISCRIMINATE BETTER in the rotated frame, and points
moving oblique->aligned should get WORSE. Content is identical; only the
grid relationship changes.
"""
import os, glob, sys, json, time
import numpy as np, tifffile, torch
from scipy import ndimage as ndi

OUT = "/mnt/vesuvius/surf191_rand"
FULL = "/mnt/vesuvius/experiments/FT191/ckpt_ft_full.pth"
K = 15000
ANG = 30.0
MARGIN = 45  # exclude points near border (rotation cval + receptive field)

sys.path.insert(0, "/mnt/vesuvius")
from diag2_191 import geometry
from loader059 import load_059, predict
from diag4_auc import auc

def rot_mats(theta_deg):
    t = np.deg2rad(theta_deg)
    # rotation in the (axis1, axis2) plane, about volume center (scipy axes=(1,2))
    R = np.array([[1, 0, 0],
                  [0, np.cos(t), -np.sin(t)],
                  [0, np.sin(t), np.cos(t)]])
    return R, R.T

def main():
    net, norm, props = load_059()
    ck = torch.load(FULL, map_location="cpu", weights_only=False)
    net.load_state_dict(ck["model"], strict=True); net.cuda().eval()
    imgs = sorted(glob.glob(os.path.join(OUT, "imagesTr", "*_0000.tif")))
    rng = np.random.default_rng(11)
    R, Rinv = rot_mats(ANG)
    rec = {k: [] for k in ("isg", "al0", "al1", "p0", "p1")}
    for i, f in enumerate(imgs[:120]):   # 120 patches is plenty
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
        c = (np.array(img.shape) - 1) / 2.0
        # original-frame prediction + geometry
        p_orig = predict(net, img, 192, norm, props)
        n, mom = geometry(gt)
        # rotated volume (content pulled: out[o] = in[R(o-c)+c])
        img_rot = ndi.rotate(img.astype(np.float32), ANG, axes=(1, 2),
                             reshape=False, order=1, mode="constant",
                             cval=float(props["mean"]))
        p_rot = predict(net, img_rot, 192, norm, props)
        # sample points away from borders
        idx = np.argwhere(gt)
        inb = ((idx > MARGIN) & (idx < np.array(img.shape) - MARGIN)).all(1)
        idx = idx[inb]
        if len(idx) < 2000: continue
        pts = idx[rng.choice(len(idx), size=min(K, len(idx)), replace=False)]
        edt, inds = ndi.distance_transform_edt(~gt, return_indices=True)
        shell = (edt >= 2) & (edt <= 10)
        idxb = np.argwhere(shell)
        inbb = ((idxb > MARGIN) & (idxb < np.array(img.shape) - MARGIN)).all(1)
        idxb = idxb[inbb]
        ptsb = idxb[rng.choice(len(idxb), size=min(K, len(idxb)), replace=False)]
        anchor = np.stack([inds[cc][ptsb[:, 0], ptsb[:, 1], ptsb[:, 2]]
                           for cc in range(3)], 1)
        for pts_, apts_, isg in ((pts, pts, 1), (ptsb, anchor, 0)):
            nn = np.stack([n[cc][apts_[:, 0], apts_[:, 1], apts_[:, 2]]
                           for cc in range(3)], 1)          # [K,3] normals
            al0 = np.max(np.abs(nn), axis=1)
            nn_rot = nn @ R.T                                # normal in rotated frame
            al1 = np.max(np.abs(nn_rot), axis=1)
            # content at input p lands at output o = R(p-c)+c
            # (verified empirically against scipy.ndimage.rotate)
            o = (pts_ - c) @ R.T + c
            v1 = ndi.map_coordinates(p_rot, o.T, order=1, mode="nearest")
            v0 = p_orig[pts_[:, 0], pts_[:, 1], pts_[:, 2]]
            rec["isg"].append(np.full(len(pts_), isg, np.uint8))
            rec["al0"].append(al0); rec["al1"].append(al1)
            rec["p0"].append(v0); rec["p1"].append(v1)
        print(f"[{i+1}/120] ({time.time()-t0:.0f}s)", flush=True)
    d = {k: np.concatenate(v) for k, v in rec.items()}
    isg = d["isg"] == 1
    # causal cells: aligned->oblique and oblique->aligned
    print("\n=== H-GRID: same content, two grid orientations (model=ft_full) ===")
    cells = [("aligned->oblique", (d["al0"] >= 0.94) & (d["al1"] < 0.80)),
             ("oblique->aligned", (d["al0"] < 0.80) & (d["al1"] >= 0.94)),
             ("stays oblique", (d["al0"] < 0.80) & (d["al1"] < 0.80)),
             ("stays aligned", (d["al0"] >= 0.94) & (d["al1"] >= 0.94))]
    out = {}
    for name, m in cells:
        a0 = auc(d["p0"][m & isg], d["p0"][m & ~isg])
        a1 = auc(d["p1"][m & isg], d["p1"][m & ~isg])
        out[name] = [a0, a1, int((m & isg).sum())]
        print(f"  {name:18s} n={int((m&isg).sum()):>9,}  AUC orig={a0:.3f} -> rot={a1:.3f}  "
              f"(delta {a1-a0:+.3f})")
    json.dump(out, open(os.path.join(OUT, "eval_ft", "hgrid_rotation.json"), "w"), indent=1)
    print("\nHGRID COMPLETE", flush=True)

if __name__ == "__main__":
    main()
