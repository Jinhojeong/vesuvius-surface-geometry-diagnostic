#!/usr/bin/env python3
"""Does recall track the intensity axes that define the located/non-located split?

Joins the per-volume intensity endpoints (intensity_split_892.json, the four
prespecified measures from the 08-07 grey-level analysis) to the six-arm recall
table (per_volume.csv) over the located60 cohort, and reports tie-corrected
Spearman correlations of recall against each covariate under two normalizations.

The point is not the size of any one coefficient. It is that recall is
uncorrelated with these axes under the checkpoint's own plans normalization and
correlated with them under per-volume z-scoring, so the split reads as a
property of the model on one path and of the intensity on the other.

Range restriction: the located60 cohort spans a narrower contrast range than the
full located population, so these coefficients are attenuated relative to what a
fair draw would give. The direction and the contrast between arms are the
readable part; the magnitudes are not.

Usage: covariates.py [intensity_split_892.json] [per_volume.csv] [out.json]
Pure stdlib, deterministic.
"""
import csv
import json
import sys

INT = sys.argv[1] if len(sys.argv) > 1 else "results/intensity_split_892.json"
CSV = sys.argv[2] if len(sys.argv) > 2 else "results/histmatch_confound/per_volume.csv"
OUT = sys.argv[3] if len(sys.argv) > 3 else "results/histmatch_confound/covariates.json"

COVARIATES = ("bg_median", "sheet_median", "contrast", "iqr_all")
ARMS = {"A_orig_ctplans": "plans CTNormalization",
        "B_orig_instzs": "instance z-score"}


def ranks(xs):
    """Average ranks, ties shared."""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    out = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        r = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            out[order[k]] = r
        i = j + 1
    return out


def spearman(xs, ys):
    rx, ry = ranks(xs), ranks(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    return num / (dx * dy) if dx and dy else float("nan")


def main():
    intens = {r["sample"]: r for r in json.load(open(INT))["rows"]}
    rows = []
    with open(CSV) as f:
        for r in csv.DictReader(f):
            if r.get("cohort") != "located60":
                continue
            s = r["sample"]
            if s not in intens:
                continue
            rec = {"sample": s}
            ok = True
            for arm in ARMS:
                v = r.get(arm + "_recall") or r.get(arm)
                if v in (None, "", "None"):
                    ok = False
                    break
                rec[arm] = float(v)
            for c in COVARIATES:
                v = intens[s].get(c)
                if v is None:
                    ok = False
                    break
                rec[c] = float(v)
            if ok:
                rows.append(rec)

    out = {"n": len(rows),
           "cohort": "located60",
           "note": __doc__.strip().split("\n\n")[2],
           "spearman": {}}
    for arm, label in ARMS.items():
        out["spearman"][arm] = {"normalization": label}
        for c in COVARIATES:
            out["spearman"][arm][c] = round(
                spearman([r[c] for r in rows], [r[arm] for r in rows]), 4)
    # range restriction, stated rather than implied
    def q(v, f):
        v = sorted(v)
        i = f * (len(v) - 1)
        lo = int(i)
        return v[lo] + (i - lo) * (v[min(lo + 1, len(v) - 1)] - v[lo])
    out["covariate_spread_in_cohort"] = {
        c: {"min": min(r[c] for r in rows), "q1": q([r[c] for r in rows], 0.25),
            "q3": q([r[c] for r in rows], 0.75), "max": max(r[c] for r in rows)}
        for c in COVARIATES}
    json.dump(out, open(OUT, "w"), indent=1)
    for arm in ARMS:
        print(arm, {c: out["spearman"][arm][c] for c in COVARIATES})
    print("n =", len(rows), "->", OUT)


if __name__ == "__main__":
    main()
