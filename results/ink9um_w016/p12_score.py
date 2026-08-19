"""Score an ink probability map inside the w016 validation mask.

Endpoints fixed by PREREGISTRATION.md (8ceda1e7a0536613):
  primary   F1 at the threshold that maximises F1, with the threshold reported
  secondary AUC, and F1 at 0.5 after the tutorial rescaling (p-0.25)/0.5
"""
import json, sys
import numpy as np

LBL = "/mnt/vesuvius/ink9um_w016/w016_val_labels.npz"


def load_labels():
    d = np.load(LBL)
    z = int(np.argmax((d["validation"] > 0).reshape(d["validation"].shape[0], -1).sum(1)))
    m = d["validation"][z] > 0
    y = (d["ink"][z] > 0)
    return m, y


def auc(scores, y):
    o = np.argsort(scores, kind="mergesort")
    s, yy = scores[o], y[o]
    n1 = int(yy.sum()); n0 = len(yy) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = np.empty(len(s), float)
    i = 0
    r = 1
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        ranks[i:j + 1] = (r + (r + (j - i))) / 2.0
        r += (j - i + 1)
        i = j + 1
    return float((ranks[yy].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def score(prob_map, tag):
    m, y = load_labels()
    if prob_map.shape != m.shape:
        raise SystemExit("shape mismatch %s vs mask %s" % (prob_map.shape, m.shape))
    p = prob_map[m].astype(np.float32)
    t = y[m]
    best = (-1, None)
    for thr in np.arange(0.05, 0.96, 0.01):
        pred = p >= thr
        tp = float((pred & t).sum()); fp = float((pred & ~t).sum()); fn = float((~pred & t).sum())
        f1 = 0.0 if (2 * tp + fp + fn) == 0 else 2 * tp / (2 * tp + fp + fn)
        if f1 > best[0]:
            best = (f1, float(thr))
    r = np.clip((p - 0.25) / 0.5, 0, 1)
    pred05 = r >= 0.5
    tp = float((pred05 & t).sum()); fp = float((pred05 & ~t).sum()); fn = float((~pred05 & t).sum())
    f105 = 0.0 if (2 * tp + fp + fn) == 0 else 2 * tp / (2 * tp + fp + fn)
    return dict(tag=tag, n=int(m.sum()), pos_rate=round(float(t.mean()), 5),
                f1_best=round(best[0], 5), thr_best=round(best[1], 3),
                auc=round(auc(p, t), 5), f1_at_05_rescaled=round(f105, 5),
                calling_rate_at_best=round(float((p >= best[1]).mean()), 5))


if __name__ == "__main__":
    import tifffile
    out = []
    for path in sys.argv[1:]:
        img = tifffile.imread(path).astype(np.float32) / 255.0
        out.append(score(img, path.split("/")[-1]))
        print(json.dumps(out[-1]))
