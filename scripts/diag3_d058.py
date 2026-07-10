"""#191 diag3 — same 200 patches, same strata, SECOND model (nnUNet surface_recto).

Reuses geometry/curvature/spacing from diag2_191. Inference via nnUNetPredictor
(fold_0, no TTA for speed). If the custom trainer class is missing from vanilla
nnunetv2, fall back to trainer_name="nnUNetTrainer" (loss-only trainers don't
change architecture).
"""
import os, glob, sys, json, time, csv
import numpy as np, tifffile, torch

OUT = "/mnt/vesuvius/surf191_rand"
DIAG = os.path.join(OUT, "diag3"); os.makedirs(DIAG, exist_ok=True)
MODEL = "/mnt/vesuvius/models/surface_recto_nnu"
TH = 0.4

sys.path.insert(0, "/mnt/vesuvius")
from diag2_191 import geometry, curvature_at, spacing_at

def make_predictor():
    from nnunetv2.inference.predict_from_raw_data import nnUNetPredictor
    import nnunetv2.inference.predict_from_raw_data as prd
    pred = nnUNetPredictor(tile_step_size=0.75, use_mirroring=False,
                           use_gaussian=True, verbose=False,
                           device=torch.device("cuda"),
                           allow_tqdm=False)
    try:
        pred.initialize_from_trained_model_folder(MODEL, use_folds=(0,),
                                                  checkpoint_name="checkpoint_final.pth")
    except Exception as e:
        print("first init failed:", type(e).__name__, str(e)[:200])
        # patch trainer_name inside checkpoint and retry
        ck = torch.load(f"{MODEL}/fold_0/checkpoint_final.pth",
                        map_location="cpu", weights_only=False)
        print("ckpt trainer_name:", ck.get("trainer_name"))
        ck["trainer_name"] = "nnUNetTrainer"
        torch.save(ck, f"{MODEL}/fold_0/checkpoint_final.pth")
        pred.initialize_from_trained_model_folder(MODEL, use_folds=(0,),
                                                  checkpoint_name="checkpoint_final.pth")
    return pred

def main():
    imgs = sorted(glob.glob(os.path.join(OUT, "imagesTr", "*_0000.tif")))
    print(f"{len(imgs)} patches", flush=True)
    pred = make_predictor()
    props = {"spacing": [1.0, 1.0, 1.0]}
    all_cv = []; all_sp = []; all_hit = []
    rows = []
    rng = np.random.default_rng(0)   # same seed as diag2 -> same sample points
    for i, f in enumerate(imgs):
        lf = os.path.join(OUT, "labelsTr", os.path.basename(f).replace("_0000.tif", ".tif"))
        if not os.path.exists(lf): continue
        try:
            img = tifffile.imread(f); lbl = tifffile.imread(lf)
        except Exception as e:
            print("skip unreadable", os.path.basename(f), e); continue
        if img.ndim != 3 or lbl.shape != img.shape: continue
        gt = lbl > 0
        if gt.sum() < 5000: continue
        t0 = time.time()
        prob = pred.predict_single_npy_array(img[None].astype(np.float32), props,
                                             None, None, True)[1][1]  # P(class1)
        p = prob > TH
        inter = (p & gt).sum()
        dice = 2*inter/(p.sum()+gt.sum()+1e-6); rec = inter/gt.sum()
        n, mom = geometry(gt)
        idx = np.argwhere(gt)
        take = rng.choice(len(idx), size=min(30000, len(idx)), replace=False)
        pts = idx[take]
        cv = curvature_at(pts, mom)
        sp = spacing_at(pts, gt, n)
        hit = p[pts[:, 0], pts[:, 1], pts[:, 2]]
        all_cv.append(cv); all_sp.append(sp); all_hit.append(hit)
        scroll = os.path.basename(f).split("_")[0]
        rows.append((os.path.basename(f), scroll, float(dice), float(rec)))
        print(f"[{i+1}/{len(imgs)}] {os.path.basename(f)[:30]:30s} {scroll} "
              f"dice={dice:.3f} rec={rec:.3f} ({time.time()-t0:.0f}s)", flush=True)
    with open(os.path.join(DIAG, "per_patch.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["patch","scroll","dice","recall"]); w.writerows(rows)
    cv = np.concatenate(all_cv); sp = np.concatenate(all_sp)
    hit = np.concatenate(all_hit)
    def table(vals, bins, name):
        print(f"\n=== recall by {name} ==="); out = []
        for b in range(len(bins)-1):
            m = (vals >= bins[b]) & (vals < bins[b+1])
            nn = int(m.sum()); r = float(hit[m].mean()) if nn else float("nan")
            out.append([float(bins[b]), float(bins[b+1]), nn, r])
            print(f"  [{bins[b]:.3f},{bins[b+1]:.3f})  n={nn:>9,}  recall={r:.3f}")
        return out
    cbins = np.unique(np.quantile(cv, np.linspace(0, 1, 11))); cbins[-1] += 1e-6
    sbins = np.array([0, 4, 6, 8, 11, 15, 22, 39, 41])
    ct = table(cv, cbins, "curvature (decile bins)")
    st = table(sp, sbins, "spacing to next sheet (vox)")
    import collections
    g = collections.defaultdict(list)
    for _, s, d, *_ in rows: g[s].append(d)
    scr = {s: [float(np.mean(v)), float(np.median(v)), len(v)] for s, v in g.items()}
    print("\n=== per-scroll dice ===")
    for s, (m, md, k) in sorted(scr.items()): print(f"  {s}: mean={m:.3f} median={md:.3f} n={k}")
    json.dump({"curvature": ct, "spacing": st, "per_scroll": scr, "th": TH},
              open(os.path.join(DIAG, "strata3.json"), "w"), indent=1)
    print("\nDIAG3 COMPLETE ->", DIAG, flush=True)

if __name__ == "__main__":
    main()
