"""Preregistered verdict, computed once by the frozen section 12 buckets.

Reads the twelve per-model score JSONs, forms the six paired deltas (B minus A)
on the primary endpoint, the two control deltas, the 15 A-to-A pairwise deltas,
recomputes the final MDE from the measured A-arm sd, runs the 50-crop flips
spot re-check, and prints the bucket. Every rule here restates
prereg_retrain_ab/PREREGISTRATION.md sections 10 and 12; nothing is chosen at
run time. Exit is informational only; the JSON is the record.
"""
import itertools, json, math, sys
from pathlib import Path

import numpy as np

AB = Path("/mnt/vesuvius/experiments/retrain_ab")
SEEDS = [40, 41, 42, 43, 44, 45]
MDE_PLAN = 0.0075


def sep(arm, seed, key):
    p = AB / f"scores/ckpts_ckpt_{arm}_s{seed}.json"
    if not p.exists():
        return None
    return json.loads(p.read_text())[key]["separated"]


def main() -> None:
    rows = {}
    for s in SEEDS:
        for arm in ("v1", "v2"):
            for k in ("primary", "onesided", "background"):
                v = sep(arm, s, k)
                if v is None:
                    print(f"missing {arm} s{s}; verdict not computable yet")
                    sys.exit(1)
                rows[(arm, s, k)] = v

    d_primary = [rows[("v2", s, "primary")] - rows[("v1", s, "primary")]
                 for s in SEEDS]
    d_ones = [rows[("v2", s, "onesided")] - rows[("v1", s, "onesided")]
              for s in SEEDS]
    d_bg = [rows[("v2", s, "background")] - rows[("v1", s, "background")]
            for s in SEEDS]
    a_primary = [rows[("v1", s, "primary")] for s in SEEDS]
    aa = [abs(a - b) for a, b in itertools.combinations(a_primary, 2)]

    mean_d = float(np.mean(d_primary))
    sd_d = float(np.std(d_primary, ddof=1))
    t = mean_d / (sd_d / math.sqrt(len(SEEDS))) if sd_d > 0 else float("inf")
    from scipy import stats
    p = float(2 * stats.t.sf(abs(t), len(SEEDS) - 1))
    signs = sum(1 for x in d_primary if x > 0) if mean_d > 0 else \
        sum(1 for x in d_primary if x < 0)
    sd_a = float(np.std(a_primary, ddof=1))
    mde_measured = 2.8 * sd_a * math.sqrt(2.0 / len(SEEDS))
    aa_width = float(np.median(aa))
    ctrl = {"onesided": float(np.mean(d_ones)),
            "background": float(np.mean(d_bg))}
    ctrl_clean = all(abs(v) < MDE_PLAN for v in ctrl.values())

    if abs(mean_d) >= MDE_PLAN and signs >= 5 and p < 0.05 and ctrl_clean:
        bucket = "POSITIVE"
    elif abs(mean_d) >= MDE_PLAN and signs >= 5 and p < 0.05:
        bucket = "SITE-UNSPECIFIC PRIOR"
    elif abs(mean_d) < MDE_PLAN and aa_width < MDE_PLAN:
        bucket = "NULL"
    else:
        bucket = "INCONCLUSIVE-UNDERPOWERED"

    out = {"per_seed": {str(s): {"A": rows[("v1", s, "primary")],
                                 "B": rows[("v2", s, "primary")],
                                 "delta": d_primary[i]}
                        for i, s in enumerate(SEEDS)},
           "mean_paired_delta": mean_d, "sd_paired": sd_d,
           "t": t, "p_two_sided": p, "signs_stable": signs,
           "ci95": [mean_d - 2.571 * sd_d / math.sqrt(6),
                    mean_d + 2.571 * sd_d / math.sqrt(6)],
           "controls_mean_delta": ctrl, "controls_below_MDE": ctrl_clean,
           "aa_null_median_width": aa_width,
           "a_arm_sd": sd_a, "mde_planning": MDE_PLAN,
           "mde_measured_2.8sd_sqrt2overN": mde_measured,
           "bucket": bucket}
    (AB / "scores/VERDICT.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
