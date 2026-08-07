"""Four-way comparison: villa's wrapper against our two normalization arms and his published rows.

  W  villa's `vesuvius.predict`, default normalization, run here      wrapper_check/W_default.jsonl
  C  villa's wrapper fed an already-CT-normalized block, --normalization none
                                                                     wrapper_check/W_ctnorm.jsonl
  P  our nnUNetPredictor path, plans CTNormalization                  norm_ablation16.json ct_plans
                                                                     (+ control_his60.jsonl, indep.)
  Z  our nnUNetPredictor path, instance z-score                       norm_ablation16.json instance_zscore
  H  TAUIL-Abd-Elilah's published rows                                surface_bench_m7.json / m7_margin_fp.json

Endpoints are [recall, precision, predicted-positive fraction] over the inner 128^3 of the centred
256^3 block, threshold 0.2, class 2 excluded. Identical arithmetic in every arm.

The margin block is the sharper test. `enrichment` and the one-voxel Euclidean shell profile are
what the public claim is actually about, and P and H disagree on the SHAPE there, not just the
level, so an arm cannot land between them by accident.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path

import numpy as np

KEYS = ("recall", "precision", "pred_pos_frac")


def from_desc(d: dict) -> list:
    """Endpoints out of geom.descriptives counts. Exact, no re-thresholding."""
    tp = d["n_pred_scored"] - d["n_fp"]
    return [round(tp / d["n_sheet"], 4), round(tp / max(d["n_pred_scored"], 1), 4),
            round(d["n_pred_scored"] / d["n_scored"], 4)]


def load_jsonl(p: Path) -> dict:
    out = {}
    if p.exists():
        for ln in p.read_text().splitlines():
            if ln.strip():
                r = json.loads(ln)
                out[r["sample"]] = r
    return out


def mad(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.median(np.abs(a - b)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--wrapper", required=True)
    ap.add_argument("--wrapper-ctnorm", default=None)
    ap.add_argument("--norm-ablation", required=True)
    ap.add_argument("--control-his60", required=True)
    ap.add_argument("--his-bench", required=True)
    ap.add_argument("--his-margin", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    W = load_jsonl(Path(a.wrapper))
    C = load_jsonl(Path(a.wrapper_ctnorm)) if a.wrapper_ctnorm else {}
    na = {r["sample"]: r for r in json.loads(Path(a.norm_ablation).read_text())["rows"]}
    p60 = load_jsonl(Path(a.control_his60))
    hb = {r["sample"]: r for r in json.loads(Path(a.his_bench).read_text())["rows"]}
    hm = {r["sample"]: r for r in json.loads(Path(a.his_margin).read_text())["rows"]}

    names = [n for n in sorted(na) if n in W and W[n].get("status") == "ok"]
    rows = []
    for n in names:
        r = {"sample": n,
             "blk_mean": na[n]["blk_mean"], "blk_std": na[n]["blk_std"],
             "W": W[n]["endpoints"], "P": na[n]["ct_plans"], "Z": na[n]["instance_zscore"],
             "H": na[n]["his"],
             "P_indep": from_desc(p60[n]["desc"]) if n in p60 else None,
             "W_from_desc": from_desc(W[n]["desc"]),
             "margin": {
                 "W_enrichment": W[n].get("enrichment"),
                 "P_enrichment": p60[n].get("enrichment") if n in p60 else None,
                 "H_enrichment": hm[n].get("enrichment") if n in hm else None,
                 "W_shells": W[n].get("shells"),
                 "P_shells": p60[n].get("shells") if n in p60 else None,
                 "H_shells": hm[n].get("shells") if n in hm else None,
                 "W_margin_share": W[n].get("margin_share_of_nonsheet"),
                 "H_margin_share": hm[n].get("margin_share_of_nonsheet") if n in hm else None,
             }}
        if n in C and C[n].get("status") == "ok":
            r["C"] = C[n]["endpoints"]
        rows.append(r)

    arms = ["W"] + (["C"] if any("C" in r for r in rows) else []) + ["P", "Z", "H"]
    g = {k: np.array([r[k] for r in rows if k in r], float) for k in arms}

    summary = {"n_volumes": len(rows), "arms": arms,
               "median": {k: [round(float(np.median(g[k][:, i])), 4) for i in range(3)] for k in arms},
               "pairwise_median_abs_dev": {}}
    for i, x in enumerate(arms):
        for y in arms[i + 1:]:
            if len(g[x]) != len(g[y]):
                continue
            summary["pairwise_median_abs_dev"][f"{x}-{y}"] = {
                KEYS[j]: round(mad(g[x][:, j], g[y][:, j]), 4) for j in range(3)}

    # margin block: median enrichment and median shell profile per arm
    def med_shell(key: str, k: int):
        v = [r["margin"][key][f"shell_{k}"] for r in rows
             if r["margin"].get(key) and r["margin"][key].get(f"shell_{k}") is not None]
        return round(float(np.median(v)), 3) if v else None

    summary["margin"] = {}
    for arm, ek, sk in (("W", "W_enrichment", "W_shells"), ("P", "P_enrichment", "P_shells"),
                        ("H", "H_enrichment", "H_shells")):
        e = [r["margin"][ek] for r in rows if r["margin"].get(ek) is not None]
        summary["margin"][arm] = {
            "median_enrichment": round(float(np.median(e)), 3) if e else None,
            "median_shells": [med_shell(sk, k) for k in range(1, 6)]}

    Path(a.out).write_text(json.dumps({"rows": rows, "summary": summary}, indent=1))

    w = 9
    print(f"\nn = {len(rows)} volumes\n")
    print("per volume, [recall / precision / pred-pos-frac]")
    hdr = "sample        " + "".join(f"{k:>26}" for k in arms)
    print(hdr)
    for r in rows:
        line = f"{r['sample']:<14}"
        for k in arms:
            v = r.get(k)
            line += "  " + (" ".join(f"{x:{w}.4f}" for x in v) if v else " " * 26)
        print(line)
    print("\nmedians")
    for k in arms:
        m = summary["median"][k]
        print(f"  {k:<3} {m[0]:.4f}  {m[1]:.4f}  {m[2]:.4f}")
    print("\npairwise median absolute deviation")
    for kk, vv in summary["pairwise_median_abs_dev"].items():
        print(f"  {kk:<6} recall {vv['recall']:.4f}   precision {vv['precision']:.4f}   ppf {vv['pred_pos_frac']:.4f}")
    print("\nmargin / shells (median over volumes)")
    for arm, d in summary["margin"].items():
        print(f"  {arm:<3} enrichment {d['median_enrichment']}   shells {d['median_shells']}")


if __name__ == "__main__":
    main()
