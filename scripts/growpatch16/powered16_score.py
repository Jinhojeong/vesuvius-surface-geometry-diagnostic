#!/usr/bin/env python3
"""powered16_score.py - score every powered16 run (new sites 8-15, arms
A/B/C) with score_ab.score_run verbatim (2000-ray census scorer + crossing
metric, near_r=6), run = unit. Identical machinery to powered_score.py /
randctl_score.py; only the run directory, site indices, and arm list
differ, so rows are directly comparable to the archived powered/randctl
rows."""
import json
import os
import sys
import time

sys.path.insert(0, "/mnt/vesuvius/vcbuild")
import score_ab as S

PW16 = "/mnt/vesuvius/vcbuild/demo_out/powered16"
SITES16 = "/mnt/vesuvius/hazard_zarr_smoke/demo_sites16.json"
NEAR_R, RAYS = 6.0, 2000

sites = json.load(open(SITES16))
tiles = S.TileCache(S.BLOCKS_DIR, S.CENSUS_DIR)
rows, missing = [], []
t_start = time.time()
for i in range(8, 16):
    s = sites[i]
    site = "s%d_%s_%s" % (i, s["slab"], s["tile"])
    for arm in ("A", "B", "C"):
        for rep in range(1, 6):
            d = os.path.join(PW16, site, "%s_r%d" % (arm, rep), "out")
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
                  "design": "8 NEW sites (idx 8-15, demo_sites16.json) x "
                            "3 arms (A=field only, B=field+"
                            "p1218_conf_v2_amp4.zarr, C=field+"
                            "p1218_conf_rand_amp4_ext8.zarr) x 5 fresh reps",
                  "control": "random field, matched count/sigma/amp/"
                             "amp_scale 4, positions redrawn uniformly in "
                             "the NEW region cubes (radius 288 L1, "
                             "rng 20260812)"},
       "n_runs": len(rows), "missing": missing, "runs": rows}
json.dump(out, open(PW16 + "/powered16_scores_new8.json", "w"), indent=1)
print("wrote %s (%d runs, %d missing, %.0fs)"
      % (PW16 + "/powered16_scores_new8.json", len(rows), len(missing),
         time.time() - t_start), flush=True)
