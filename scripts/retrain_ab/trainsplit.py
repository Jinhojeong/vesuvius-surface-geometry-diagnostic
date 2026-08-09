"""Gate list: SPLIT sites on TRAIN tiles, 60 hash-ordered tiles, cap 2000.
The pre-freeze pilot is scored here so the frozen eval sites stay untouched
until the prereg is public."""
import hashlib, json
from pathlib import Path
AB = Path("/mnt/vesuvius/experiments/retrain_ab")
SRC = Path("/mnt/vesuvius/kaggle_p1218_repair_v2")
split = json.loads((AB / "frozen/tile_split.json").read_text())["split"]
train = sorted((n for n, g in split.items() if g == "train"),
               key=lambda n: hashlib.md5(n.encode()).hexdigest())[:60]
out = []
for name in train:
    rp = SRC / "records" / f"{name}.json"
    if not rp.exists():
        continue
    r = json.loads(rp.read_text())
    for s in r["sites"]:
        if s["decision"] == "SPLIT":
            out.append({"tile": name, "slab": r["slab"],
                        "z": s["z"], "y": s["y"], "x": s["x"],
                        "inst": s["inst"], "th": s["th"], "ratio": s["ratio"],
                        "n_sites": s.get("n_sites"),
                        "key": hashlib.md5(
                            f"ts:{name}:{s['i']}:{s['z']}:{s['y']}:{s['x']}"
                            .encode()).hexdigest()})
out.sort(key=lambda s: s["key"])
out = out[:2000]
(AB / "frozen/sites_trainsplit.json").write_text(json.dumps(
    {"what": "SPLIT sites on 60 hash-ordered train tiles, gate use only",
     "n": len(out), "sites": out}, indent=0))
print("trainsplit sites:", len(out), "tiles:", len({s["tile"] for s in out}))
