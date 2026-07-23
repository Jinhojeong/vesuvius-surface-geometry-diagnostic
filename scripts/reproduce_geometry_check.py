"""Reproduce the geometric-channel check: does local CT geometry carry any
weld information?

This is the cheapest of the measurements behind the #191 write-up and the one
most worth re-running, because it is the obvious next idea (feed the model
geometry instead of intensity) and it fails for a reason that is easy to miss.

The discrimination: a thick bright band in the CT is either
  WELD  - two sheets in physical contact (GT shows two runs under 4 voxels)
  THICK - one genuinely thick sheet (GT shows one run, above median thickness)
Both look like a single band to a binary model. If structure-tensor or curvature
features separate them, a geometric input channel could supply the boundary that
intensity lacks.

Features come from the raw CT only, never from GT; GT is used solely to label
the two classes. Three controls run alongside the headline number, because the
headline number is misleading on its own: the two classes differ in band
thickness by construction, and the profile statistics encode thickness.

Expected output (8 Dataset059 patches, n=674 as shipped):
  full feature set        AUC 0.79    <- looks like a pass
  thickness alone         AUC 0.98    <- the classes are separable by thickness
  structure tensor only   AUC 0.51    <- chance; geometry carries nothing
  thickness-matched       AUC 0.48    <- chance

Usage:
  python scripts/reproduce_geometry_check.py --data /path/to/patches
where the directory holds imagesTr/{id}_0000.tif (uint8 raw CT) and
labelsTr/{id}.tif (binary GT), in the layout of the patches-v1 release.

The numbers above need eight crops. patches-v1 ships two, which is enough to
check the code runs but not to reproduce the table: with two volumes the
per-class cap leaves n around 170, and the thickness-matched control falls
under the 60-pair minimum below, so that row prints nan while the other three
print normally. The matched row is the one the conclusion rests on.
"""
import argparse
import glob
import json
import os

import numpy as np
import tifffile
from scipy import ndimage as ndi

SPAN, STEP = 8.0, 0.5
K = int(2 * SPAN / STEP) + 1
PER_CLASS_PER_VOL = 60


def normals_of(mask):
    d_out = ndi.distance_transform_edt(~mask).astype(np.float32)
    d_in = ndi.distance_transform_edt(mask).astype(np.float32)
    sm = ndi.gaussian_filter(d_out - d_in, 1.0)
    n = np.stack(np.gradient(sm), 0).astype(np.float32)
    n /= (np.sqrt((n ** 2).sum(0)) + 1e-6)
    return n


def runs(profile):
    out, i = [], 0
    while i < len(profile):
        if profile[i]:
            j = i
            while j < len(profile) and profile[j]:
                j += 1
            out.append(((i + j - 1) / 2.0, i, j - 1))
            i = j
        else:
            i += 1
    return out


def st_features(ct, p0, sigma, R=5):
    """structure tensor and Hessian curvature of the CT around p0"""
    z, y, x = p0
    sl = (slice(max(z - R - 3, 0), z + R + 4), slice(max(y - R - 3, 0), y + R + 4),
          slice(max(x - R - 3, 0), x + R + 4))
    w = ct[sl].astype(np.float32)
    if min(w.shape) < 2 * R:
        return None
    g = np.gradient(ndi.gaussian_filter(w, sigma))
    T = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            T[a, b] = ndi.gaussian_filter(g[a] * g[b], sigma * 1.5).mean()
    l1, l2, l3 = [float(v) for v in np.sort(np.linalg.eigvalsh(T))[::-1]]
    s = l1 + 1e-9
    sm = ndi.gaussian_filter(w, sigma)
    H = np.zeros((3, 3))
    gg = np.gradient(sm)
    for a in range(3):
        gh = np.gradient(gg[a])
        for b in range(3):
            c = tuple(v // 2 for v in gh[b].shape)
            H[a, b] = gh[b][c]
    k1, k2, k3 = [float(v) for v in np.sort(np.linalg.eigvalsh((H + H.T) / 2))]
    return [(l1 - l2) / s, (l2 - l3) / s, l3 / s, (l1 - l3) / s, l1, l2, l3,
            k1, k2, k3, k1 * k3, (k1 + k3) / 2, abs(k1) / (abs(k3) + 1e-9)]


def profile_features(ct, p0, nv):
    offs = np.linspace(-SPAN, SPAN, K)
    cc = (np.array(p0, float)[None, :] + offs[:, None] * nv[None, :]).T
    pr = ndi.map_coordinates(ct.astype(np.float32), cc, order=1, mode="nearest")
    c = K // 2
    seg = pr[max(c - 14, 0):min(c + 15, K)]
    d1 = np.gradient(seg)
    d2 = np.gradient(d1)
    return [float(seg.mean()), float(seg.std()), float(seg.max() - seg.min()),
            float(seg[len(seg) // 2]), float(np.abs(d1).max()),
            float(d2.min()), float(d2.max()),
            float(np.percentile(seg, 25)), float(np.percentile(seg, 75))]


def collect(data_dir, limit=8):
    X, Y, W = [], [], []
    rng = np.random.default_rng(0)
    offs = np.linspace(-SPAN, SPAN, K)
    files = sorted(glob.glob(f"{data_dir}/imagesTr/*_0000.tif"))[:limit]
    if not files:
        raise SystemExit(f"no patches under {data_dir}/imagesTr/")
    for f in files:
        lf = f.replace("imagesTr", "labelsTr").replace("_0000.tif", ".tif")
        ct = tifffile.imread(f)
        gt = tifffile.imread(lf) > 0
        nrm = normals_of(gt)
        pts = np.argwhere(gt)
        pts = pts[rng.permutation(len(pts))[:12000]]
        m_ = int(SPAN) + 8
        pts = pts[((pts >= m_) & (pts < np.array(gt.shape) - m_)).all(1)]
        widths, cache = [], []
        for p0 in pts:
            nv = nrm[:, p0[0], p0[1], p0[2]]
            prof = ndi.map_coordinates(
                gt.astype(np.float32),
                (np.array(p0, float)[None, :] + offs[:, None] * nv[None, :]).T,
                order=0, mode="constant")
            rc = runs(prof > 0.5)
            here = [r for r in rc if r[1] <= K // 2 <= r[2]]
            if not rc or not here:
                continue
            w = (here[0][2] - here[0][1] + 1) * STEP
            widths.append(w)
            cache.append((tuple(int(v) for v in p0), nv, rc, w))
        if len(widths) < 200:
            continue
        med = float(np.median(widths))
        weld, thick = [], []
        for p0, nv, rc, w in cache:
            if len(weld) >= PER_CLASS_PER_VOL and len(thick) >= PER_CLASS_PER_VOL:
                break
            if len(rc) >= 2:
                sp = min(abs(rc[i + 1][0] - rc[i][0])
                         for i in range(len(rc) - 1)) * STEP
                if 0 < sp < 4.0 and len(weld) < PER_CLASS_PER_VOL:
                    bw = (max(r[2] for r in rc) - min(r[1] for r in rc) + 1) * STEP
                    weld.append((p0, nv, bw))
            elif len(rc) == 1 and w > 1.5 * med and len(thick) < PER_CLASS_PER_VOL:
                thick.append((p0, nv, w))
        for label, group in ((1, weld), (0, thick)):
            for p0, nv, bw in group:
                feats, ok = [], True
                for sg in (1.0, 2.0, 3.0):
                    fv = st_features(ct, p0, sg)
                    if fv is None:
                        ok = False
                        break
                    feats += fv
                if not ok:
                    continue
                feats += profile_features(ct, p0, nv)
                if not np.all(np.isfinite(feats)):
                    continue
                X.append(feats)
                Y.append(label)
                W.append(bw)
        print(f"{os.path.basename(f)[:26]}: weld {len(weld)} thick {len(thick)} "
              f"(median band {med:.1f} vox)", flush=True)
    return np.array(X), np.array(Y), np.array(W)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True,
                    help="directory with imagesTr/ and labelsTr/")
    ap.add_argument("--out", default=None, help="optional JSON output path")
    args = ap.parse_args()
    X, Y, W = collect(args.data)
    print(f"\ncollected n={len(Y)} ({int(Y.sum())} weld / {int((1-Y).sum())} thick)")
    if len(Y) < 60:
        raise SystemExit("not enough samples; point --data at more patches")

    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    def cv_auc(Xa, Ya):
        aucs = []
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=0).split(Xa, Ya):
            clf = make_pipeline(StandardScaler(),
                                LogisticRegression(max_iter=2000))
            clf.fit(Xa[tr], Ya[tr])
            aucs.append(roc_auc_score(Ya[te], clf.predict_proba(Xa[te])[:, 1]))
        return float(np.mean(aucs)), float(np.std(aucs))

    full = cv_auc(X, Y)
    geom = cv_auc(X[:, :39], Y)          # structure tensor + curvature only
    thick_only = roc_auc_score(Y, W)
    thick_only = max(thick_only, 1 - thick_only)

    # thickness-matched pairs: the decisive control
    wi, ti = np.where(Y == 1)[0], np.where(Y == 0)[0]
    used, pairs = set(), []
    for i in wi:
        cand = [j for j in ti if j not in used and abs(W[j] - W[i]) <= 0.5]
        if cand:
            j = min(cand, key=lambda j: abs(W[j] - W[i]))
            used.add(j)
            pairs += [i, j]
    matched = cv_auc(X[pairs], Y[pairs]) if len(pairs) >= 60 else (float("nan"),) * 2

    print(f"\nfull feature set          AUC {full[0]:.3f} +- {full[1]:.3f}")
    print(f"thickness alone           AUC {thick_only:.3f}")
    print(f"structure tensor only     AUC {geom[0]:.3f} +- {geom[1]:.3f}")
    print(f"thickness-matched (n={len(pairs)})  AUC {matched[0]:.3f}")
    if len(pairs) < 60:
        print(f"\nNOT REPRODUCED: the thickness-matched control needs 60 pairs and "
              f"got {len(pairs)}, so the row above is nan. That control is what "
              f"separates a real geometric signal from a thickness artifact, and "
              f"without it the other three rows do not support any conclusion. "
              f"Re-run on eight crops; two are not enough.")
    else:
        print("\nReading: the full-set number is a thickness artifact. Geometry alone "
              "is at chance, so a geometric input channel has no weld signal to give "
              "a model.")
    if args.out:
        json.dump({"n": int(len(Y)), "full": full, "thickness_only": thick_only,
                   "geometry_only": geom, "thickness_matched": matched,
                   "matched_n": len(pairs)}, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
