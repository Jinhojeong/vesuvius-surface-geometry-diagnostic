"""The control gate: pass criteria written down before the run, then the verdict.

Tolerance and why it is what it is
----------------------------------
We are not reproducing his inference stack bit for bit and never could: he drives villa's
`vesuvius.predict` wrapper under a Windows conda env, we drive nnUNetv2 2.8.1's own
sliding window. Patch placement, blend weights and fp16 rounding all differ a little, so
an exact match is not the bar. The bar is that the statistic his control rests on comes
out the same to a precision that still separates the two hypotheses that control was
built to separate: a STEP at shell 1 versus a SMOOTH DECAY.

  G1  every one of the five per-volume median shell enrichments within 10% relative
  G2  per-volume Pearson and Spearman correlation of shell-1 enrichment >= 0.95
  G3  per-volume median obs/pred within 0.05 absolute of his 0.7554
  G4  per-volume median margin enrichment within 10% relative of his 3.104
  G5  ordering shell_1 > shell_2 > ... > shell_5 preserved, and median obs/pred on the
      same side of 1 as his

10% on shell 1 is +-0.31 on a value of 3.10. The gap his control turns on is shell 1
(3.10) against shell 2 (1.85), which is 1.25. A 10% band is a fifth of that gap, so it
still discriminates; anything looser stops discriminating and the gate stops being a gate.
"""
import json, sys
from pathlib import Path
import numpy as np

HIS = {"shells": {"shell_1": 3.103, "shell_2": 1.851, "shell_3": 1.012,
                  "shell_4": 0.588, "shell_5": 0.306},
       "obs_over_pred_median": 0.7554, "margin_enrichment_median": 3.104,
       "frac_obs_over_pred_above_1": 0.10}

s = json.loads(Path(sys.argv[1]).read_text())["summary"]
ours_sh = s["shell_enrichment_median"]
v = s["vs_his"]

g1 = {}
for k in range(1, 6):
    h, o = HIS["shells"][f"shell_{k}"], ours_sh[f"shell_{k}"]
    g1[f"shell_{k}"] = {"his": h, "ours": o, "rel_err": round((o - h) / h, 4),
                        "pass": abs(o - h) / h <= 0.10}
g2 = {"pearson": v["shell_1"]["pearson_r"], "spearman": v["shell_1"]["spearman_r"],
      "pass": v["shell_1"]["pearson_r"] >= 0.95 and v["shell_1"]["spearman_r"] >= 0.95}
op = s["obs_over_pred"]["median"]
g3 = {"his": HIS["obs_over_pred_median"], "ours": op,
      "abs_diff": round(abs(op - HIS["obs_over_pred_median"]), 4),
      "pass": abs(op - HIS["obs_over_pred_median"]) <= 0.05}
me = s["margin_enrichment"]["median"]
g4 = {"his": HIS["margin_enrichment_median"], "ours": me,
      "rel_err": round((me - HIS["margin_enrichment_median"]) / HIS["margin_enrichment_median"], 4),
      "pass": abs(me - HIS["margin_enrichment_median"]) / HIS["margin_enrichment_median"] <= 0.10}
order = all(ours_sh[f"shell_{k}"] > ours_sh[f"shell_{k+1}"] for k in range(1, 5))
g5 = {"monotone_decay": order,
      "his_obs_over_pred_below_1": HIS["obs_over_pred_median"] < 1,
      "ours_obs_over_pred_below_1": op < 1,
      "same_side_of_1": (op < 1) == (HIS["obs_over_pred_median"] < 1),
      "pass": order and (op < 1) == (HIS["obs_over_pred_median"] < 1)}

out = {"criteria_doc": __doc__,
       "G1_shell_medians_within_10pct": g1,
       "G2_shell1_per_volume_correlation_ge_0.95": g2,
       "G3_obs_over_pred_within_0.05": g3,
       "G4_margin_enrichment_within_10pct": g4,
       "G5_shape_and_sign_preserved": g5}
out["G1_pass"] = all(x["pass"] for x in g1.values())
out["verdict"] = "PASS" if all([out["G1_pass"], g2["pass"], g3["pass"], g4["pass"], g5["pass"]]) else "FAIL"
out["failed"] = [n for n, p in [("G1", out["G1_pass"]), ("G2", g2["pass"]), ("G3", g3["pass"]),
                                ("G4", g4["pass"]), ("G5", g5["pass"])] if not p]
print(json.dumps(out, indent=1))
Path(sys.argv[2]).write_text(json.dumps(out, indent=1))
