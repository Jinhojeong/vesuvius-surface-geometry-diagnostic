"""Background negative control: on-mask points far from every census site.

The +30y offset control measured badly: 69 percent of displaced sites leave the
mask, and the survivors sit near other repair material, moving +0.062 under the
v2 oracle. This control replaces it in the primary gate. Points are sampled on
the v1 mask of eval tiles, at least 40 vox from every census site of that tile
(any decision), hash-ordered, capped at 5,000; the eligible population is smaller and 2,265 are kept. Same schema as the other site lists.
"""
import hashlib, json
import numpy as np
from pathlib import Path

AB = Path("/mnt/vesuvius/experiments/retrain_ab")
SRC = Path("/mnt/vesuvius/kaggle_p1218_repair_v2")
V1 = Path("/mnt/vesuvius/p1218_full/blocks")
CAP, DMIN = 5000, 40

split = json.loads((AB / "frozen/tile_split.json").read_text())["split"]
eval_tiles = sorted(n for n, g in split.items() if g == "eval")
out = []
for name in eval_tiles:
    slab, tile = name.split("_", 1)
    rp = SRC / "records" / f"{name}.json"
    sites = json.loads(rp.read_text())["sites"] if rp.exists() else []
    sxyz = np.array([[s["z"], s["y"], s["x"]] for s in sites]) if sites else None
    try:
        with np.load(V1 / slab / f"{tile}.npz") as d:
            lab = d["labels"]
    except FileNotFoundError:
        continue
    mask = lab > 0
    rng = np.random.default_rng(int(hashlib.md5(name.encode()).hexdigest()[:8], 16))
    fg = np.argwhere(mask[32:-32:4, 32:-32:4, 32:-32:4])
    if not len(fg):
        continue
    picks = fg[rng.permutation(len(fg))[:200]] * 4 + 32
    kept = 0
    for p in picks:
        if kept >= 40:
            break
        if sxyz is not None and len(sxyz) and np.min(
                np.linalg.norm(sxyz - p, axis=1)) < DMIN:
            continue
        z, y, x = int(p[0]), int(p[1]), int(p[2])
        out.append({"tile": name, "slab": slab, "z": z, "y": y, "x": x,
                    "inst": int(lab[z, y, x]), "th": None, "ratio": None,
                    "n_sites": None,
                    "key": hashlib.md5(f"bg:{name}:{z}:{y}:{x}".encode()).hexdigest()})
        kept += 1
out.sort(key=lambda s: s["key"])
out = out[:CAP]
(AB / "frozen/sites_background.json").write_text(json.dumps(
    {"what": f"on-mask, >= {DMIN} vox from every census site, hash order, cap {CAP}",
     "n": len(out), "sites": out}, indent=0))
print("background sites:", len(out))
