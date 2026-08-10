"""Prereg section 9 secondary: instance VOI with VOI_merge reported separately.

For each of the six models (seeds 40-42, both arms): predict the first 20
evaluation tiles in hash order, threshold at that model's own frozen-run t*,
form connected components of (foreground minus predicted boundary), and score
VOI against both reference segmentations, the components of foreground minus
each arm's boundary target. voi_split = H(PR|GT), voi_merge = H(GT|PR), the
convention fixed in villa PR #1301. One JSON, medians per (model, reference).
"""
import hashlib, json, sys
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

sys.path.insert(0, "/mnt/vesuvius/experiments/retrain_ab")
from instrument import predict_tile, boundary_target, V1, V2, AB

import torch
sys.path.insert(0, "/mnt/vesuvius")
from loader059 import load_059

SEEDS = (40, 41, 42)


def components(fg, boundary):
    lab, n = ndi.label(fg & ~boundary)
    return lab


def voi_terms(gt, pr):
    """H(PR|GT) and H(GT|PR) over the joint foreground, nats."""
    m = (gt > 0) & (pr > 0)
    if not m.any():
        return None, None
    g, p = gt[m].ravel(), pr[m].ravel()
    joint = {}
    for gi, pi in zip(g[::17], p[::17]):          # 1/17 stride subsample
        joint[(gi, pi)] = joint.get((gi, pi), 0) + 1
    n = sum(joint.values())
    pg, pp = {}, {}
    for (gi, pi), c in joint.items():
        pg[gi] = pg.get(gi, 0) + c
        pp[pi] = pp.get(pi, 0) + c
    h_pr_gt = -sum(c / n * np.log(c / pg[gi]) for (gi, pi), c in joint.items())
    h_gt_pr = -sum(c / n * np.log(c / pp[pi]) for (gi, pi), c in joint.items())
    return float(h_pr_gt), float(h_gt_pr)


def main() -> None:
    split = json.loads((AB / "frozen/tile_split.json").read_text())["split"]
    tiles = sorted((n for n, g in split.items() if g == "eval"),
                   key=lambda n: hashlib.md5(n.encode()).hexdigest())[:20]
    net, norm, props = load_059()
    mean, std = props["mean"], props["std"]
    out = {}
    for arm in ("v1", "v2"):
        for seed in SEEDS:
            tag = f"{arm}_s{seed}"
            sc = json.loads((AB / f"scores/ckpts_ckpt_{tag}.json").read_text())
            t_star = sc["t_star"]
            ck = torch.load(AB / f"ckpts/ckpt_{tag}.pth", map_location="cpu",
                            weights_only=False)
            net.load_state_dict(ck["model"], strict=True)
            net.cuda().eval()
            rows = []
            for name in tiles:
                slab, tile = name.split("_", 1)
                with np.load(V1 / slab / f"{tile}.npz") as d:
                    l1 = d["labels"].astype(np.int32)
                with np.load(V2 / slab / f"{tile}.npz") as d:
                    l2 = d["labels"].astype(np.int32)
                fg = l1 > 0
                ct = np.load(AB / "evalct" / f"{name}.npy")
                ct = ct[:l1.shape[0], :l1.shape[1], :l1.shape[2]]
                p = predict_tile(net, ct, mean, std)
                pred = components(fg, p > t_star)
                row = {"tile": name}
                for rname, rlab in (("v1", l1), ("v2", l2)):
                    ref = components(rlab > 0, boundary_target(rlab))
                    sp, mg = voi_terms(ref, pred)
                    row[f"voi_split_vs_{rname}"] = sp
                    row[f"voi_merge_vs_{rname}"] = mg
                rows.append(row)
                del p, pred
            med = {k: float(np.median([r[k] for r in rows if r[k] is not None]))
                   for k in rows[0] if k != "tile"}
            out[tag] = {"t_star": t_star, "medians": med, "rows": rows}
            print(tag, json.dumps(med), flush=True)
    (AB / "scores/VOI_SECONDARY.json").write_text(json.dumps(out, indent=1))
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
