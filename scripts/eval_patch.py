"""Patch-mode eval entry point for external splitters (issue: IyanDopico/
vesuvius-sheet-tools#1).

Scores a prediction against a GT surface patch with the same harness used in
the #191 diagnostics, and reports the metric family separately so topology
changes are visible even when the blend is not.

Inputs
  --gt      GT surface label patch: .tif (binary uint8) or .npy
  --pred    prediction: .tif/.npy binary mask, or .npz with an int 'labels'
            array (instance ids; mask = labels > 0)
  --pred-b  optional second prediction (e.g. splitter ON) -> prints an
            A -> B delta table; A is --pred (e.g. splitter OFF)

Reported per prediction
  blend / topo / surface-dice / voi   (official blend, subprocess-isolated)
  components (26-conn) for pred and GT, largest-component share
  if instance labels were given: instance count and instance-level VOI
  against GT 26-conn components (informative only for segment-style
  instances: expect over-segmentation by construction; read the A -> B
  delta, not the absolute value)

Patch layout this consumes (see repo README):
  imagesTr/{scroll}_z{Z}_y{Y}_x{X}_0000.tif  uint8 raw CT, 300^3, zyx
  labelsTr/{scroll}_z{Z}_y{Y}_x{X}.tif       uint8 binary GT, same grid
"""
import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def load_any(path):
    """Returns (mask bool array, instance labels or None)."""
    if path.endswith(".npz"):
        with np.load(path) as d:
            key = "labels" if "labels" in d else list(d.keys())[0]
            lab = d[key]
        return lab > 0, lab.astype(np.int64)
    if path.endswith(".npy"):
        a = np.load(path)
    else:
        import tifffile
        a = tifffile.imread(path)
    if a.dtype.kind in "iu" and a.max() > 1:
        return a > 0, a.astype(np.int64)
    return a > 0.5 if a.dtype.kind == "f" else a > 0, None


def official_scores(gt, pred):
    """Subprocess-isolated official blend (a Betti abort cannot kill us)."""
    here = os.path.dirname(os.path.abspath(__file__))
    with tempfile.TemporaryDirectory() as td:
        gp, pp = os.path.join(td, "g.npy"), os.path.join(td, "p.npy")
        np.save(gp, gt)
        np.save(pp, pred)
        r = subprocess.run([sys.executable, os.path.join(here, "official_one.py"),
                            gp, pp], capture_output=True, text=True, timeout=600)
        v = [float(x) for x in r.stdout.split()]
    return v if len(v) == 4 else [float("nan")] * 4


def component_stats(mask):
    from scipy import ndimage as ndi
    lab, n = ndi.label(mask, structure=np.ones((3, 3, 3)))
    if n == 0:
        return 0, 0.0
    sizes = np.bincount(lab.ravel())[1:]
    return int(n), float(sizes.max() / max(mask.sum(), 1))


def instance_voi(gt_mask, inst):
    """VOI between GT 26-conn components and predicted instances, on the
    intersection of both foregrounds. Lower is better."""
    from scipy import ndimage as ndi
    gl, _ = ndi.label(gt_mask, structure=np.ones((3, 3, 3)))
    m = (gl > 0) & (inst > 0)
    if m.sum() < 100:
        return float("nan")
    a = gl[m].astype(np.int64)
    b = inst[m].astype(np.int64)
    n = float(len(a))

    def ent(x):
        c = np.bincount(x)
        p = c[c > 0] / n
        return float(-(p * np.log2(p)).sum())

    joint = a.astype(np.int64) * (b.max() + 1) + b
    _, joint = np.unique(joint, return_inverse=True)
    return 2 * ent(joint.astype(np.int64)) - ent(a) - ent(b)


def report(tag, gt, mask, inst):
    blend, topo, sdice, voi = official_scores(gt, mask)
    ncomp, share = component_stats(mask)
    gcomp, _ = component_stats(gt)
    row = {"blend": blend, "topo": topo, "sdice": sdice, "voi": voi,
           "components": ncomp, "gt_components": gcomp,
           "largest_share": share}
    if inst is not None:
        row["instances"] = int(len(np.unique(inst[inst > 0])))
        row["instance_voi"] = instance_voi(gt, inst)
    print(f"\n[{tag}]")
    for k, v in row.items():
        print(f"  {k:15s} {v:.4f}" if isinstance(v, float) else f"  {k:15s} {v}")
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gt", required=True)
    ap.add_argument("--pred", required=True)
    ap.add_argument("--pred-b", default=None,
                    help="second prediction for an A->B delta (e.g. ON)")
    args = ap.parse_args()
    gt, _ = load_any(args.gt)
    mask_a, inst_a = load_any(args.pred)
    if gt.shape != mask_a.shape:
        raise SystemExit(f"shape mismatch: gt {gt.shape} vs pred {mask_a.shape}")
    ra = report("pred" if not args.pred_b else "A (off)", gt, mask_a, inst_a)
    if args.pred_b:
        mask_b, inst_b = load_any(args.pred_b)
        rb = report("B (on)", gt, mask_b, inst_b)
        print("\n[A -> B delta]")
        for k in ra:
            if k in rb and isinstance(ra[k], float) and isinstance(rb[k], float):
                print(f"  {k:15s} {rb[k]-ra[k]:+.4f}")
            elif k in rb:
                print(f"  {k:15s} {rb[k]-ra[k]:+d}")


if __name__ == "__main__":
    main()
