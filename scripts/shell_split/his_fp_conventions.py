"""His published shell enrichments, converted to false-positive mass shares.

His `results/m7_margin_fp.json` publishes per-shell ENRICHMENT, which is a ratio of two
shares and carries no mass on its own. Multiplying a shell's enrichment by that shell's
share of the scored non-sheet volume recovers the share of his false-positive mass that
sits in the shell. The shell volumes are a label-only quantity, so they are identical in
his run and ours and can be read straight off our control rows.

Two conventions are reported because they are not the same number. A per-volume median is
the median over volumes of a quantity formed on each volume first. A pooled figure weights
each volume by its false-positive count. Sums of medians are not medians, so every
aggregate here is formed per volume and only then reduced.

What it needs
  --his    his published results/m7_margin_fp.json
           (github.com/TAUIL-Abd-Elilah/vesuvius-repro)
  --ours   control_his60.jsonl written by run_shells.py over the same 60 volumes
  --out    destination json

    python his_fp_conventions.py --his m7_margin_fp.json \
        --ours control_his60.jsonl --out his_fp_mass_conventions.json
"""
from __future__ import annotations
import argparse, json
from pathlib import Path

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--his", required=True)
ap.add_argument("--ours", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

his = {r["sample"]: r for r in json.loads(Path(a.his).read_text())["rows"]
       if r.get("status") == "ok"}
ours = {json.loads(l)["sample"]: json.loads(l)
        for l in Path(a.ours).read_text().splitlines() if l.strip()}
common = sorted(set(his) & set(ours))

imp, w = [], []
for s in common:
    d = ours[s]["desc"]                      # shell_vox is label-only, identical to his
    sv = [v / d["n_nonsheet_scored"] for v in d["shell_vox"]]
    hs = [his[s]["shells"][f"shell_{k}"] for k in range(1, 6)]
    imp.append([h * v for h, v in zip(hs, sv)])
    w.append(his[s]["n_fp"])
imp = np.array(imp) * 100.0
w = np.array(w, float)

per_volume_shells = np.median(imp, axis=0)
per_volume_45 = np.median(imp[:, 3] + imp[:, 4])
beyond3 = 100.0 - imp[:, :3].sum(axis=1)     # shells 4, 5 and everything past 5
pooled = (imp * w[:, None]).sum(axis=0) / w.sum()

out = {
    "n_volumes": len(common),
    "his_implied_fp_mass_share_pct": {
        "per_volume_median": [round(float(x), 2) for x in per_volume_shells],
        "per_volume_median_shells_4_plus_5": round(float(per_volume_45), 2),
        "per_volume_median_beyond_shell3": round(float(np.median(beyond3)), 2),
        "per_volume_q1_q3_beyond_shell3": [round(float(np.percentile(beyond3, 25)), 2),
                                           round(float(np.percentile(beyond3, 75)), 2)],
        "fp_weighted_pooled": [round(float(x), 2) for x in pooled],
        "fp_weighted_pooled_shells_4_plus_5": round(float(pooled[3] + pooled[4]), 2),
        "fp_weighted_pooled_beyond_shell3": round(float(100.0 - pooled[:3].sum()), 2),
        "note": ("His quoted 11 percent for the region past shell 3 sits between the two "
                 "readings of his own rows, 10.32 for shells 4 and 5 together and 12.04 for "
                 "everything past shell 3, so it holds. Both are per-volume medians, formed "
                 "on each volume and reduced afterwards. Subtracting the first three shell "
                 "medians from 100 instead gives 13.76, which is not a median of anything, "
                 "so it is not reported here."),
    },
    "his_quoted": [42, 27, 16, 11],
}
Path(a.out).write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
