"""Prereg build, day 1: freeze the tile split and the site lists.

Emits, under /mnt/vesuvius/experiments/retrain_ab/frozen/:
  tile_split.json     train/eval assignment for every tile. The eval-list names
                      (200 in eval_tiles_197.json, 15 without records) are eval
                      by fiat; the
                      remaining tiles split 60/40 train/heldout by a
                      deterministic md5 of the tile name, no RNG anywhere.
  sites_primary.json  held-out SPLIT sites on eval tiles, recast-pass filter
                      deferred to the recast join (TODO noted in file), capped
                      at 5,000 by deterministic hash order.
  sites_onesided.json 2,000 ONE_SIDED sites on eval tiles, same ordering.
  sites_offset.json   the primary list displaced +30 vox in y, flagged where
                      the displacement leaves the 128^3 tile and dropped there.
  counts.json         every population count this freeze saw, for the prereg.

Nothing here reads labels or CT. It is a pure join over records/ and the
shipped eval-tile list, so it is reproducible from the public v2.0 artifacts.
"""
from __future__ import annotations
import hashlib, json
from pathlib import Path

SRC = Path("/mnt/vesuvius/kaggle_p1218_repair_v2")
OUT = Path("/mnt/vesuvius/experiments/retrain_ab/frozen")
OUT.mkdir(parents=True, exist_ok=True)

Y_BOUND = 512   # blocks are (256, 512, 512) int32, measured
OFFSET_Y = 30
CAP_PRIMARY, CAP_ONESIDED = 5000, 2000


def h(name: str) -> str:
    return hashlib.md5(name.encode()).hexdigest()


def main() -> None:
    ev = json.loads((SRC / "validation/eval_tiles_197.json").read_text())
    eval_names = set(ev["tiles"])
    # tiles the repair never touched have no record file; eval membership is
    # still by fiat, they simply contribute no sites
    split_pre = {n: "eval" for n in eval_names}

    split, counts = dict(split_pre), {"eval": len(split_pre), "train": 0,
                                       "heldout_nontrain": 0}
    sites = {"SPLIT": [], "ONE_SIDED": []}
    dec_totals: dict[str, int] = {}

    for rp in sorted((SRC / "records").glob("*.json")):
        r = json.loads(rp.read_text())
        name = f"{r['slab']}_{r['tile']}"
        if name in eval_names:
            grp = "eval"
            counts[grp] -= 1        # pre-seeded above; do not double count
        else:
            grp = "train" if int(h(name), 16) % 100 < 60 else "heldout_nontrain"
        split[name] = grp
        counts[grp] += 1
        for s in r["sites"]:
            d = s.get("decision", "?")
            dec_totals[d] = dec_totals.get(d, 0) + 1
            if grp == "eval" and d in sites:
                sites[d].append({
                    "tile": name, "slab": r["slab"],
                    "z": s["z"], "y": s["y"], "x": s["x"],
                    "inst": s["inst"], "th": s["th"], "ratio": s["ratio"],
                    "n_sites": s.get("n_sites"),
                    "key": h(f"{name}:{s['i']}:{s['z']}:{s['y']}:{s['x']}"),
                })

    for k in sites:
        sites[k].sort(key=lambda s: s["key"])   # deterministic, name-hash order

    primary = sites["SPLIT"][:CAP_PRIMARY]
    onesided = sites["ONE_SIDED"][:CAP_ONESIDED]
    offset, off_dropped = [], 0
    for s in primary:
        y = s["y"] + OFFSET_Y
        if y >= Y_BOUND:
            off_dropped += 1
            continue
        offset.append({**s, "y": y, "offset_of": s["key"]})

    (OUT / "tile_split.json").write_text(json.dumps(
        {"rule": "eval = eval_tiles_197.json by fiat; others md5(name) mod 100 "
                 "< 60 -> train else heldout_nontrain",
         "split": split}, indent=0))
    (OUT / "sites_primary.json").write_text(json.dumps(
        {"what": "SPLIT sites on eval tiles, deterministic-hash order, "
                 f"cap {CAP_PRIMARY}. TODO-MEASURE: recast-pass filter joins "
                 "repairs/ before freeze.",
         "n": len(primary), "sites": primary}, indent=0))
    (OUT / "sites_onesided.json").write_text(json.dumps(
        {"n": len(onesided), "sites": onesided}, indent=0))
    (OUT / "sites_offset.json").write_text(json.dumps(
        {"offset_y": OFFSET_Y, "n": len(offset), "n_dropped_out_of_tile":
         off_dropped, "sites": offset}, indent=0))
    counts_all = {
        "tiles": counts,
        "decisions_all_tiles": dec_totals,
        "eval_SPLIT_available": len(sites["SPLIT"]),
        "eval_ONE_SIDED_available": len(sites["ONE_SIDED"]),
        "primary_frozen": len(primary),
        "onesided_frozen": len(onesided),
        "offset_frozen": len(offset),
    }
    (OUT / "counts.json").write_text(json.dumps(counts_all, indent=1))
    print(json.dumps(counts_all, indent=1))


if __name__ == "__main__":
    main()
