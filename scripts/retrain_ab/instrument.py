"""Prereg build, day 4b: the read-out instrument.

For one model (or one oracle arm), over the frozen site lists: predict the
boundary probability for every eval tile that carries sites, set ONE global
threshold by the matched-budget rule, then score each site by whether the
thresholded boundary intersects its flagging ray's single-id thick run.

Ray geometry is validate_v4.py's, verbatim: normal from the v1 mask's SDF in a
33-vox window, run bounds walked on the mask along the normal, sample points at
STEP spacing. The mask is byte-identical across arms so the geometry is
arm-independent by construction.

Budget rule (frozen): threshold t* is the smallest 1/2000 grid value whose
pooled predicted-positive fraction over ALL scored voxels of the site-carrying
eval tiles is at or below b*, where b* is the v2 boundary-target fraction over
the same voxels. Ties inside a bin are not re-ordered; the bin boundary is the
threshold, which is conservative in the direction of no separation.

Modes:
  oracle-v1 / oracle-v2   the arm's own boundary target as the "prediction",
                          threshold ignored (targets are binary). Unit test:
                          oracle-v1 must separate ~0 of the primary sites,
                          oracle-v2 must land near the recast pass rate.
  model:<ckpt path>       GPU, sliding 160^3 window, 80-vox step, uniform blend.

Memory shape: the first version held every tile's probability map for a second
pass and died at 46 GB on a 31 GB box. The verdict only needs the MAXIMUM
probability along each site's run points, since thresholding then intersecting
equals comparing that maximum to the threshold. So pass 1 keeps the pooled
histogram and one float per site and frees each map immediately.

Output: one JSON per (mode) under experiments/retrain_ab/scores/.
"""
from __future__ import annotations
import json, os, sys, time
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

AB = Path("/mnt/vesuvius/experiments/retrain_ab")
V1 = Path("/mnt/vesuvius/p1218_full/blocks")
V2 = Path("/mnt/vesuvius/kaggle_p1218_repair_v2/blocks_repaired")
SC = AB / "scores"
SC.mkdir(exist_ok=True)

K, STEP = 33, 0.75
SPAN = (K // 2) * STEP
OFFS = np.arange(K, dtype=np.float32) * STEP - SPAN
W = 16
NBIN = 2000
PS = 160


def sdf(mask):
    return (ndi.distance_transform_edt(mask)
            - ndi.distance_transform_edt(~mask)).astype(np.float32)


def crop_of(shape, p0, w):
    sl = tuple(slice(max(0, p0[i] - w), min(shape[i], p0[i] + w + 1))
               for i in range(3))
    c = tuple(p0[i] - sl[i].start for i in range(3))
    return sl, c


def normal_at(mask, p0, w=W):
    sl, c = crop_of(mask.shape, p0, w)
    sm = sdf(mask[sl])
    n = np.array([np.gradient(sm, axis=a)[c] for a in range(3)], np.float32)
    return n / (np.linalg.norm(n) + 1e-6), sl, c


def run_bounds(maskf, p0, nv):
    cc = (np.asarray(p0, float)[None, :] + OFFS[:, None] * nv[None, :]).T
    on = ndi.map_coordinates(maskf, cc, order=0, mode="constant") > 0.5
    c = K // 2
    if not on[c]:
        return None
    a = c
    while a > 0 and on[a - 1]:
        a -= 1
    b = c
    while b < K - 1 and on[b + 1]:
        b += 1
    return a, b


def run_points(p0, nv, a, b, shape):
    kk = np.arange(a, b + 1)
    q = np.round(np.asarray(p0, float)[None, :]
                 + (kk * STEP - SPAN)[:, None] * nv[None, :]).astype(int)
    ok = ((q >= 0) & (q < np.array(shape)).all() if False else
          ((q >= 0) & (q < np.array(shape))).all(1))
    return q[ok]


def boundary_target(lab):
    out = np.zeros(lab.shape, bool)
    fg = lab > 0
    for ax in range(3):
        a = [slice(None)] * 3
        b = [slice(None)] * 3
        a[ax] = slice(None, -1)
        b[ax] = slice(1, None)
        a, b = tuple(a), tuple(b)
        la, lb = lab[a], lab[b]
        d = (la > 0) & (lb > 0) & (la != lb)
        out[a] |= d
        out[b] |= d
    return out & fg


def predict_tile(net, ct, mean, std):
    import torch
    a = ((ct.astype(np.float32) - mean) / std)
    Z, Y, X = a.shape
    prob = np.zeros(a.shape, np.float32)
    wsum = np.zeros(a.shape, np.float32)
    stp = PS // 2
    zs = sorted(set(list(range(0, Z - PS + 1, stp)) + [Z - PS]))
    ys = sorted(set(list(range(0, Y - PS + 1, stp)) + [Y - PS]))
    xs = sorted(set(list(range(0, X - PS + 1, stp)) + [X - PS]))
    with torch.no_grad():
        for z0 in zs:
            for y0 in ys:
                for x0 in xs:
                    t = torch.from_numpy(
                        a[z0:z0+PS, y0:y0+PS, x0:x0+PS][None, None]).cuda()
                    with torch.autocast("cuda", torch.float16):
                        o = net(t)
                    if isinstance(o, dict):
                        o = o["surface"]
                    if isinstance(o, (list, tuple)):
                        o = o[0]
                    p = torch.softmax(o.float(), 1)[0, 1].cpu().numpy()
                    prob[z0:z0+PS, y0:y0+PS, x0:x0+PS] += p
                    wsum[z0:z0+PS, y0:y0+PS, x0:x0+PS] += 1.0
    return prob / np.maximum(wsum, 1e-6)


def main() -> None:
    mode = sys.argv[1]
    tag = mode.replace("model:", "").replace("/", "_").replace(".pth", "") \
        if mode.startswith("model:") else mode
    names = os.environ.get("LISTS", "primary,onesided,offset,background").split(",")
    lists = {n: json.loads((AB / "frozen" / f"sites_{n}.json").read_text())["sites"]
             for n in names}
    tag = tag + os.environ.get("TAG_SUFFIX", "")
    by_tile: dict[str, dict[str, list]] = {}
    for n, ss in lists.items():
        for s in ss:
            s["_id"] = f"{n}:{s['key']}"      # offset rows reuse the primary key
            by_tile.setdefault(s["tile"], {}).setdefault(n, []).append(s)

    net = mean = std = None
    if mode.startswith("model:"):
        import torch
        sys.path.insert(0, "/mnt/vesuvius")
        from loader059 import load_059
        net, norm, props = load_059()
        ck = torch.load(mode[6:], map_location="cpu", weights_only=False)
        net.load_state_dict(ck["model"], strict=True)
        net.cuda().eval()
        mean, std = props["mean"], props["std"]

    # single pass: pooled histogram + per-site max probability along the run
    hists = np.zeros(NBIN, np.int64)
    n_scored = 0
    b_star_num = 0
    site_max: dict[str, tuple[str, float]] = {}   # _id -> (cause, max prob)
    tiles = sorted(by_tile)
    t0 = time.time()
    for i, name in enumerate(tiles):
        slab, tile = name.split("_", 1)
        with np.load(V2 / slab / f"{tile}.npz") as d:
            l2 = d["labels"].astype(np.int32)
        scored = l2 > 0                      # foreground = scored region
        if mode == "oracle-v2":
            p = boundary_target(l2).astype(np.float32)
        elif mode == "oracle-v1":
            with np.load(V1 / slab / f"{tile}.npz") as d:
                p = boundary_target(d["labels"].astype(np.int32)).astype(np.float32)
        else:
            ct = np.load(AB / "evalct" / f"{name}.npy")
            # edge tiles are smaller than the standard block and the CT cache
            # is zero-padded to it; crop back to the label shape or the scored
            # boolean index breaks (post-freeze amendment, recorded in hashes)
            ct = ct[:l2.shape[0], :l2.shape[1], :l2.shape[2]]
            p = predict_tile(net, ct, mean, std)
        hists += np.histogram(p[scored], bins=NBIN, range=(0.0, 1.0))[0]
        n_scored += int(scored.sum())
        b_star_num += int(boundary_target(l2).sum())
        del scored
        with np.load(V1 / slab / f"{tile}.npz") as d:
            mask = d["labels"] > 0
        maskf = mask.astype(np.float32)
        for n, ss in by_tile[name].items():
            for s in ss:
                p0 = (int(s["z"]), int(s["y"]), int(s["x"]))
                if not all(0 <= p0[k] < mask.shape[k] for k in range(3)) \
                        or not mask[p0]:
                    site_max[s["_id"]] = ("off_mask", -1.0)
                    continue
                nv, sl, c = normal_at(mask, p0)
                rb = run_bounds(maskf[sl], c, nv)
                if rb is None:
                    site_max[s["_id"]] = ("no_run", -1.0)
                    continue
                q = run_points(c, nv, *rb, mask[sl].shape)
                off = np.array([sl[k].start for k in range(3)])
                qq = q + off
                site_max[s["_id"]] = ("ok", float(p[qq[:, 0], qq[:, 1], qq[:, 2]].max()))
        del p, mask, maskf, l2
        if (i + 1) % 25 == 0:
            print(f"pass [{i+1}/{len(tiles)}] {(time.time()-t0)/60:.1f}min",
                  flush=True)

    b_star = b_star_num / max(n_scored, 1)
    tail = np.concatenate([np.cumsum(hists[::-1])[::-1], [0]])
    idx = int(np.argmax(tail / max(n_scored, 1) <= b_star))
    t_star = idx / NBIN
    print(f"b*={b_star:.6f} t*={t_star:.5f}", flush=True)

    thr = 0.5 if mode.startswith("oracle") else t_star
    out = {n: [] for n in lists}
    causes = {n: {"off_mask": 0, "no_run": 0} for n in lists}
    for n, ss in lists.items():
        for s in ss:
            cause, mx = site_max[s["_id"]]
            v = -1 if mx < 0 else int(mx > thr)
            if v < 0:
                causes[n][cause] += 1
            out[n].append({"key": s["key"], "verdict": v, "max_run_prob": mx})

    res = {"mode": mode, "b_star": b_star, "t_star": t_star,
           "n_tiles": len(tiles)}
    for n in out:
        ok = [r for r in out[n] if r["verdict"] >= 0]
        res[n] = {"n": len(out[n]), "n_evaluable": len(ok),
                  "not_evaluable": causes[n],
                  "separated": (sum(r["verdict"] for r in ok) / len(ok))
                  if ok else None}
    (SC / f"{tag}.json").write_text(json.dumps(
        {**res, "rows": out}, indent=0))
    print(json.dumps(res, indent=1), flush=True)


if __name__ == "__main__":
    main()
