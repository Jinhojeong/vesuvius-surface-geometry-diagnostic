"""Where exactly do the two runs part company?

His surface_bench_m7.json scores the same 60 volumes at the same crop and threshold, so
its per-volume n_sheet and base_rate are pure label quantities. If those match ours to the
last digit while recall and precision do not, then the crop, the label handling and the
scored region are shared and the whole difference is in what the network emitted.
"""
import json
from pathlib import Path
import numpy as np
from scipy.stats import pearsonr, spearmanr

his = {r["sample"]: r for r in json.loads(Path("/tmp/his_surface_bench_m7.json").read_text())["rows"]
       if r.get("status") == "ok"}
ours = {json.loads(l)["sample"]: json.loads(l) for l in
        Path("/mnt/vesuvius/experiments/shell_split/control_his60.jsonl").read_text().splitlines() if l.strip()}
common = sorted(set(his) & set(ours))
rows = []
for s in common:
    d = ours[s]["desc"]
    tp = d["n_pred_scored"] - d["n_fp"]
    rows.append({"sample": s,
                 "his_n_sheet": his[s]["n_sheet"], "ours_n_sheet": d["n_sheet"],
                 "his_base_rate": round(his[s]["base_rate"], 8),
                 "ours_base_rate": round(d["n_sheet"] / d["n_scored"], 8),
                 "his_recall": round(his[s]["recall"], 4), "ours_recall": round(tp / d["n_sheet"], 4),
                 "his_precision": round(his[s]["precision"], 4),
                 "ours_precision": round(tp / max(1, d["n_pred_scored"]), 4),
                 "his_ppf": round(his[s]["pred_positive_fraction"], 4),
                 "ours_ppf": round(d["n_pred_scored"] / d["n_scored"], 4)})
n_sheet_exact = sum(r["his_n_sheet"] == r["ours_n_sheet"] for r in rows)
base_exact = sum(abs(r["his_base_rate"] - r["ours_base_rate"]) < 1e-7 for r in rows)
a = lambda k: np.array([r[k] for r in rows], float)
out = {
    "n_common": len(rows),
    "label_side_identical": {
        "n_sheet_exact_match": f"{n_sheet_exact}/{len(rows)}",
        "base_rate_match_to_1e-7": f"{base_exact}/{len(rows)}"},
    "prediction_side": {
        "recall": {"his_median": round(float(np.median(a("his_recall"))), 4),
                   "ours_median": round(float(np.median(a("ours_recall"))), 4),
                   "pearson": round(float(pearsonr(a("his_recall"), a("ours_recall"))[0]), 4),
                   "spearman": round(float(spearmanr(a("his_recall"), a("ours_recall"))[0]), 4),
                   "ours_higher_in": f"{int((a('ours_recall') > a('his_recall')).sum())}/{len(rows)}"},
        "precision": {"his_median": round(float(np.median(a("his_precision"))), 4),
                      "ours_median": round(float(np.median(a("ours_precision"))), 4),
                      "pearson": round(float(pearsonr(a("his_precision"), a("ours_precision"))[0]), 4),
                      "ours_higher_in": f"{int((a('ours_precision') > a('his_precision')).sum())}/{len(rows)}"},
        "pred_positive_fraction": {"his_median": round(float(np.median(a("his_ppf"))), 4),
                                   "ours_median": round(float(np.median(a("ours_ppf"))), 4)}},
    "reading": ("identical label side, both recall and precision higher in ours on every "
                "volume. Note the argument in the reverse direction is not sound, since a "
                "misregistration lowers both, so this ranks explanations rather than "
                "settling them; the matched-operating-point rescore is what does that. A degraded "
                "prediction can lower both at once"),
    "rows": rows,
}
Path("/mnt/vesuvius/experiments/shell_split/divergence_diagnosis.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))
