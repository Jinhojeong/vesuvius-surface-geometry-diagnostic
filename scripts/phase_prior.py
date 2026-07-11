"""H2 test: phase-demodulation prior for compressed sheet regions (CPU-only).

Hypothesis: where sheets are tightly packed (<8 vox), voxel-contrast detection
fails (models: A below chance 0.41, B 0.77) but the stack's quasi-PERIODICITY
survives. Score each point by band-passed (at the locally estimated period)
signal along the sheet normal; if per-spacing-bin AUC of this model-free score
beats the learned models in tight bins, the phase prior carries information the
networks miss -> build it in as a feature channel / canonicalization.

Protocol mirrors diag4 (GT vs near-shell background points, same spacing bins).
"""
import os, glob, sys, json, time
import numpy as np, tifffile
from scipy import ndimage as ndi

OUT = "/mnt/vesuvius/surf191_rand"
DIAG = os.path.join(OUT, "phase_prior"); os.makedirs(DIAG, exist_ok=True)
K = 30000          # match diag4 exactly -> rng sequence aligns point-for-point
T = 24             # half-length of the normal profile
PERIODS = np.arange(4, 21)  # candidate periods (vox)

sys.path.insert(0, "/mnt/vesuvius")
from diag2_191 import geometry, spacing_at

def profiles_along_normal(img, pts, n):
    """[K, 2T+1] intensity profiles along +/- normal at pts."""
    ts = np.arange(-T, T + 1, dtype=np.float32)
    pz, py, px = pts[:, 0], pts[:, 1], pts[:, 2]
    nz = n[0][pz, py, px]; ny = n[1][pz, py, px]; nx = n[2][pz, py, px]
    cz = pz[:, None] + ts[None] * nz[:, None]
    cy = py[:, None] + ts[None] * ny[:, None]
    cx = px[:, None] + ts[None] * nx[:, None]
    coords = np.stack([cz.ravel(), cy.ravel(), cx.ravel()])
    v = ndi.map_coordinates(img.astype(np.float32), coords, order=1, mode="nearest")
    return v.reshape(len(pts), 2 * T + 1)

# precompute windowed cosine (morlet-like real) kernels per period
_KER = {}
for p in PERIODS:
    t = np.arange(-T, T + 1, dtype=np.float32)
    env = np.exp(-(t**2) / (2 * (1.2 * p) ** 2))
    k = env * np.cos(2 * np.pi * t / p)
    k -= k.mean()
    _KER[p] = (k / (np.linalg.norm(k) + 1e-9)).astype(np.float32)
KMAT = np.stack([_KER[p] for p in PERIODS])       # [P, 2T+1]

def phase_scores(prof):
    """Per point: detrend, estimate the local period by AUTOCORRELATION
    (phase-invariant), then return the PHASE-LOCKED matched-filter response at
    the center for that single period. Discriminative signal = is there a peak
    exactly at t=0, not merely 'is the texture periodic'."""
    sm = ndi.gaussian_filter1d(prof, 8.0, axis=1)
    x = prof - sm
    x = x / (x.std(axis=1, keepdims=True) + 1e-6)
    L = x.shape[1]
    # autocorrelation via FFT, per point
    F = np.fft.rfft(x, n=2*L, axis=1)
    ac = np.fft.irfft(F * np.conj(F), axis=1)[:, :L]
    ac = ac / (ac[:, :1] + 1e-9)
    pstar = PERIODS[np.argmax(ac[:, PERIODS], axis=1)]        # [K]
    pidx = np.searchsorted(PERIODS, pstar)
    return np.einsum("kt,kt->k", x, KMAT[pidx])               # phase-locked @ t=0

def main():
    imgs = sorted(glob.glob(os.path.join(OUT, "imagesTr", "*_0000.tif")))
    rng = np.random.default_rng(0)
    allsp = []; allgt = []; allsc = []
    for i, f in enumerate(imgs):
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
        n, mom = geometry(gt)
        idx = np.argwhere(gt)
        take = rng.choice(len(idx), size=min(K, len(idx)), replace=False)
        pts_gt = idx[take]
        edt, inds = ndi.distance_transform_edt(~gt, return_indices=True)
        shell = (edt >= 2) & (edt <= 10)
        idxb = np.argwhere(shell)
        takeb = rng.choice(len(idxb), size=min(K, len(idxb)), replace=False)
        pts_bg = idxb[takeb]
        anchor = np.stack([inds[c][pts_bg[:, 0], pts_bg[:, 1], pts_bg[:, 2]]
                           for c in range(3)], 1)
        for pts, apts, isg in ((pts_gt, pts_gt, 1), (pts_bg, anchor, 0)):
            prof = profiles_along_normal(img, pts, n)
            sc = phase_scores(prof)
            sp = spacing_at(apts, gt, n)
            allsp.append(sp); allsc.append(sc)
            allgt.append(np.full(len(pts), isg, np.uint8))
        print(f"[{i+1}/{len(imgs)}] {os.path.basename(f)[:26]:26s} ({time.time()-t0:.0f}s)",
              flush=True)
    sp = np.concatenate(allsp); sc = np.concatenate(allsc)
    isg = np.concatenate(allgt) == 1
    from diag4_auc import auc
    sbins = np.array([0, 4, 6, 8, 11, 15, 22, 39, 41])
    print("\n=== PHASE-PRIOR per-spacing-bin AUC (compare: A=0.41/0.56/0.67, B=0.77/0.76/0.79 in first 3 bins) ===")
    tab = []
    for b in range(len(sbins)-1):
        m = (sp >= sbins[b]) & (sp < sbins[b+1])
        a = auc(sc[m & isg], sc[m & ~isg])
        tab.append([float(sbins[b]), float(sbins[b+1]), int((m & isg).sum()),
                    int((m & ~isg).sum()), a])
        print(f"  [{sbins[b]:5.1f},{sbins[b+1]:5.1f})  nP={int((m&isg).sum()):>8,} "
              f"nN={int((m&~isg).sum()):>8,}  AUC_phase={a:.3f}")
    json.dump(tab, open(os.path.join(DIAG, "phase_auc.json"), "w"), indent=1)
    # complementarity vs learned models (points align with diag4 via same rng)
    d4 = np.load("/mnt/vesuvius/surf191_rand/diag4/points.npz")
    if len(d4["is_gt"]) == len(isg) and (d4["is_gt"] == isg.astype(np.uint8)).all():
        z = lambda v: (v - v.mean()) / (v.std() + 1e-9)
        combos = {"pA": d4["pA"], "pB": d4["pB"], "phase": sc,
                  "pA+phase": z(d4["pA"]) + z(sc), "pB+phase": z(d4["pB"]) + z(sc)}
        print("\n=== COMPLEMENTARITY: per-spacing-bin AUC ===")
        hdr = "  bin           " + "".join(f"{k:>10s}" for k in combos)
        print(hdr); ctab = {k: [] for k in combos}
        for b in range(len(sbins)-1):
            m = (sp >= sbins[b]) & (sp < sbins[b+1])
            cells = ""
            for k, v in combos.items():
                a = auc(v[m & isg], v[m & ~isg]); ctab[k].append(a)
                cells += f"{a:>10.3f}"
            print(f"  [{sbins[b]:4.0f},{sbins[b+1]:4.0f})" + cells)
        json.dump(ctab, open(os.path.join(DIAG, "combo_auc.json"), "w"), indent=1)
    else:
        print("WARN: diag4 point alignment failed — combo analysis skipped")
    np.savez_compressed(os.path.join(DIAG, "points_phase.npz"),
                        sp=sp, sc=sc, is_gt=isg.astype(np.uint8))
    print("\nPHASE PRIOR TEST COMPLETE ->", DIAG, flush=True)

if __name__ == "__main__":
    main()
