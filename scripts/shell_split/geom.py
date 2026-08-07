"""Label-side geometry: exact re-implementation of TAUIL-Abd-Elilah's margin and shell
definitions (thin_labels.across_sheet_dirs + margin_relabel.relabel_margin +
m7_margin_fp.distance_profile / analyse), read out of recon_repro @ 9afa412.

Nothing here touches the model. Everything is a function of (CT, labels) only.
"""
from __future__ import annotations
import numpy as np
from scipy.ndimage import gaussian_filter, distance_transform_edt

SIGMA = 1.0
TRIM = 64
SIZE = 256
THRESH = 0.2


def across_sheet_dirs(ct: np.ndarray, pts: np.ndarray) -> np.ndarray:
    """thin_labels.py:43 verbatim."""
    sm = gaussian_filter(ct.astype(np.float32), SIGMA)
    g = np.gradient(sm)
    H = np.empty((len(pts), 3, 3), dtype=np.float32)
    for i in range(3):
        gi = np.gradient(g[i])
        for j in range(3):
            H[:, i, j] = gi[j][pts[:, 0], pts[:, 1], pts[:, 2]]
    H = 0.5 * (H + np.transpose(H, (0, 2, 1)))
    _, v = np.linalg.eigh(H)
    d = v[:, :, 0]
    n = np.linalg.norm(d, axis=1, keepdims=True)
    return d / np.maximum(n, 1e-6)


def relabel_margin(ct: np.ndarray, lab: np.ndarray, chunk: int = 400_000):
    """margin_relabel.py:44 verbatim. Returns (labels_margin, stats)."""
    out = lab.copy()
    pts = np.argwhere(lab == 1)
    if len(pts) == 0:
        return out, {"n_sheet": 0, "n_relabelled": 0}
    shape = np.array(lab.shape)
    touched = np.zeros(lab.shape, dtype=bool)
    for s in range(0, len(pts), chunk):
        p = pts[s:s + chunk].astype(np.int32)
        d = across_sheet_dirs(ct, p)
        for sign in (1.0, -1.0):
            q = np.rint(p + sign * d).astype(np.int64)
            np.clip(q, 0, shape - 1, out=q)
            touched[q[:, 0], q[:, 1], q[:, 2]] = True
    sel = touched & (lab == 0)
    out[sel] = 2
    return out, {"n_sheet": int((lab == 1).sum()), "n_relabelled": int(sel.sum())}


def margin_mask(ct: np.ndarray, lab: np.ndarray) -> np.ndarray:
    """m7_margin_fp.margin_mask: (labels_margin == 2) & (labels == 0)."""
    lm, _ = relabel_margin(ct, lab)
    return (lm == 2) & (lab == 0)


def distance_profile(lab: np.ndarray, pred: np.ndarray, max_d: int = 5) -> dict:
    """m7_margin_fp.py:97 verbatim."""
    sheet, ignore = lab == 1, lab == 2
    scored = ~ignore
    d = distance_transform_edt(~sheet)
    fp = pred & scored & ~sheet
    n_fp, n_ns = float(fp.sum()), float((scored & ~sheet).sum())
    if n_fp == 0 or n_ns == 0:
        return {}
    out = {}
    for k in range(1, max_d + 1):
        shell = (d > k - 1) & (d <= k) & scored & ~sheet
        n_shell = float(shell.sum())
        if n_shell < 100:
            out[f"shell_{k}"] = None
            continue
        out[f"shell_{k}"] = round((float((fp & shell).sum()) / n_fp) / (n_shell / n_ns), 3)
    return out


def analyse(lab: np.ndarray, margin: np.ndarray, pred: np.ndarray) -> dict:
    """m7_margin_fp.py:136 verbatim."""
    sheet, ignore = lab == 1, lab == 2
    scored = ~ignore
    fp = pred & scored & ~sheet
    non_sheet = scored & ~sheet
    n_fp = float(fp.sum())
    n_ns = float(non_sheet.sum())
    n_margin = float((margin & scored).sum())
    if n_fp == 0 or n_ns == 0 or n_margin == 0:
        return {"status": "degenerate", "n_fp": n_fp, "n_margin": n_margin}
    share_fp = float((fp & margin).sum()) / n_fp
    share_vol = n_margin / n_ns
    return {
        "status": "ok",
        "fp_share_in_margin": round(share_fp, 6),
        "margin_share_of_nonsheet": round(share_vol, 6),
        "enrichment": round(share_fp / share_vol, 3),
        "margin_hit_rate": round(float((pred & margin & scored).sum()) / n_margin, 5),
        "nonmargin_fp_rate": round((n_fp - float((fp & margin).sum()))
                                   / max(n_ns - n_margin, 1.0), 5),
        "n_fp": int(n_fp), "n_margin_scored": int(n_margin),
    }


def descriptives(lab: np.ndarray, pred: np.ndarray, max_d: int = 5) -> dict:
    """NOT in his code. FP mass share by shell, and the beyond-k share of scored
    non-sheet volume. Raw counts kept so both pooled and per-volume conventions
    can be formed downstream."""
    sheet, ignore = lab == 1, lab == 2
    scored = ~ignore
    d = distance_transform_edt(~sheet)
    ns = scored & ~sheet
    fp = pred & ns
    n_fp, n_ns = int(fp.sum()), int(ns.sum())
    row = {"n_fp": n_fp, "n_nonsheet_scored": n_ns,
           "n_scored": int(scored.sum()), "n_sheet": int(sheet.sum()),
           "n_pred_scored": int((pred & scored).sum()),
           "shell_fp": [], "shell_vox": []}
    for k in range(1, max_d + 1):
        sh = (d > k - 1) & (d <= k) & ns
        row["shell_fp"].append(int((fp & sh).sum()))
        row["shell_vox"].append(int(sh.sum()))
    for k in (2, 3, 4):
        beyond = (d > k) & ns
        row[f"beyond{k}_vox"] = int(beyond.sum())
        row[f"beyond{k}_fp"] = int((fp & beyond).sum())
    return row
