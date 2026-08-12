#!/usr/bin/env python3
"""placement_score.py - score every placement run with score_ab.score_run
verbatim (2000-ray census scorer + crossing metric, near_r=6), run = unit.
Identical machinery to powered_score.py / powered16_score.py; only the run
directory, site list (off_sites.txt / offseeds_placement.json), and site
names (o0..o7) differ, so rows are directly comparable to the archived
powered rows. score_ab.py is NOT modified."""
import json
import os
import sys
import time

sys.path.insert(0, "/mnt/vesuvius/vcbuild")
import score_ab as S

PLACE = "/mnt/vesuvius/vcbuild/demo_out/placement"
OFFSEEDS = "/mnt/vesuvius/hazard_zarr_smoke/offseeds_placement.json"
NEAR_R, RAYS = 6.0, 2000

offs = json.load(open(OFFSEEDS))
tiles = S.TileCache(S.BLOCKS_DIR, S.CENSUS_DIR)
rows, missing = [], []
t_start = time.time()
for i, o in enumerate(offs):
    site = "o%d_%s_%s" % (i, o["slab"], o["tile"])
    for arm in ("A", "B"):
        for rep in range(1, 6):
            d = os.path.join(PLACE, site, "%s_r%d" % (arm, rep), "out")
            run = S.load_tifxyz(d)
            if run is None:
                missing.append(d)
                print("MISSING", d, flush=True)
                continue
            m = S.score_run(run, tiles, NEAR_R, RAYS)
            m.update({"site": site, "site_idx": i, "arm": arm, "rep": rep})
            rows.append(m)
            ry, cr = m["ray"], m["crossing"]
            print("%-34s %s r%d q=%5d dbl=%-6s on=%4d/%4d near=%-6s area=%.0f"
                  % (site, arm, rep, m["n_valid_quads"],
                     ry.get("frac_double"), ry.get("n_on_sheet", 0),
                     ry.get("n_sampled", 0), cr.get("frac_quads_near"),
                     m.get("area_vx2") or 0), flush=True)

out = {"generated": time.strftime("%F %T"),
       "params": {"near_r": NEAR_R, "ray_samples": RAYS,
                  "scorer": "score_ab.score_run (powered-identical)",
                  "design": "8 off-seeds (offseeds_placement.json, "
                            "prereg_placement_ab @ 8531276, Amendment 1 @ "
                            "083f424) x 2 arms (A=field only, B=field+"
                            "p1218_conf_v2_amp4.zarr) x 5 fresh reps"},
       "n_runs": len(rows), "missing": missing, "runs": rows}
json.dump(out, open(PLACE + "/placement_scores.json", "w"), indent=1)
print("wrote %s (%d runs, %d missing, %.0fs)"
      % (PLACE + "/placement_scores.json", len(rows), len(missing),
         time.time() - t_start), flush=True)
