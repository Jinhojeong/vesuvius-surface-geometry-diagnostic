"""Full-corpus m7 shell/margin pass over the public 892 volumes.

Reuses, unmodified, the shipped inference path and scorer in
/mnt/vesuvius/experiments/shell_split/gen (run_shells.py + geom.py):
  build_predictor / predict_prob  -> nnUNetPredictor, checkpoint_best.pth, fold_0,
                                     3d_fullres, patch 192^3, tile_step 0.5, gaussian,
                                     mirroring off, plans CTNormalization,
                                     centred 256^3 block, TRIM=64 -> inner 128^3,
                                     mask = sigmoid(l1 - l0) > 0.2
  margin_mask / analyse / distance_profile / descriptives  -> label-side geometry

Nothing under shell_split/ is written to. This file and its outputs live in
/mnt/vesuvius/experiments/shell892/.

Shape of the run: one GPU job at a time in the parent, a pool of CPU workers for the
scoring. margin_mask is a Hessian eigendecomposition over the full 320^3 and costs
9-18 s/volume against 2.1 s of GPU, so leaving it in the GPU loop wastes ~85% of the
wall clock. Rows are appended as they finish, so a drop costs at most the volumes in
flight and a restart resumes.
"""
from __future__ import annotations
import argparse, json, os, sys, time
from collections import deque
from pathlib import Path

GEN = "/mnt/vesuvius/experiments/shell_split/gen"
sys.path.insert(0, GEN)

import numpy as np
import tifffile

from geom import margin_mask, distance_profile, analyse, descriptives, TRIM, SIZE, THRESH

IM = Path("/mnt/vesuvius/kaggle892/images")
LB = Path("/mnt/vesuvius/kaggle892/labels")


def _init_worker() -> None:
    for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
              "NUMEXPR_NUM_THREADS"):
        os.environ[v] = "2"


# His per-volume files carry recall at these thresholds, so carry exact counts at the
# same six. NBIN is the resolution of the retained probability histogram, which is what
# lets any budget-calibrated threshold be recovered after the fact without a rerun.
THRESH_GRID = (0.05, 0.10, 0.20, 0.30, 0.50, 0.70)
NBIN = 2000


def score_one(nm: str, p: np.ndarray, lo: int, hi: int) -> dict:
    """CPU worker. Re-reads the volume (warm in page cache from the GPU pass) and runs
    the shipped scorer verbatim on the inner 128^3 crop.

    `p` arrives as float32, the same dtype run_shells.py thresholds on, so the
    `p > THRESH` decision is bit-identical to the shipped path."""
    t0 = time.time()
    ct = np.asarray(tifffile.imread(str(IM / f"{nm}.tif")))
    lab = np.asarray(tifffile.imread(str(LB / f"{nm}.tif")))
    sl = (slice(lo, hi),) * 3
    labc = lab[sl]
    mg = margin_mask(ct, lab)[sl]
    pred = p > THRESH
    r = analyse(labc, mg, pred)
    r["shells"] = distance_profile(labc, pred)
    r["desc"] = descriptives(labc, pred)
    r["sample"] = nm
    r["crop"] = [lo, hi]

    # Threshold sweep and retained histograms. Scored region only (class 2 excluded),
    # which is the whole point of the correction.
    scored = labc != 2
    sheet = labc == 1
    ps = p[scored]
    psh = p[sheet]
    r["threshold_sweep"] = {
        f"{t:.2f}": {"n_pred_scored": int((ps > t).sum()), "n_tp": int((psh > t).sum())}
        for t in THRESH_GRID}
    hist = {"sample": nm, "nbin": NBIN, "range": [0.0, 1.0],
            "n_scored": int(scored.sum()), "n_sheet": int(sheet.sum()),
            "hist_scored": np.histogram(ps, bins=NBIN, range=(0.0, 1.0))[0].tolist(),
            "hist_sheet": np.histogram(psh, bins=NBIN, range=(0.0, 1.0))[0].tolist()}
    r["_hist"] = hist
    r["cpu_seconds"] = round(time.time() - t0, 1)
    return r


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--device", default="cuda")
    a = ap.parse_args()

    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor

    names = json.loads(Path(a.names).read_text())
    outp = Path(a.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if outp.exists():
        for ln in outp.read_text().splitlines():
            if ln.strip():
                done.add(json.loads(ln)["sample"])
    todo = [n for n in names if n not in done]

    # Pool is created with 'spawn' and before CUDA is touched: the workers are pure
    # numpy/scipy and must not inherit a CUDA context.
    ctx = mp.get_context("spawn")
    pool = ProcessPoolExecutor(max_workers=a.workers, mp_context=ctx,
                               initializer=_init_worker)

    import run_shells as RS  # imports torch; created after the pool exists
    pred_obj, ip = RS.build_predictor(a.device)

    print(f"MARKER init ok  workers={a.workers}  todo={len(todo)}  already={len(done)}",
          flush=True)
    print(f"MARKER geom TRIM={TRIM} SIZE={SIZE} THRESH={THRESH}  intensity={ip}", flush=True)
    print(f"MARKER num_seg_heads={pred_obj.label_manager.num_segmentation_heads}", flush=True)

    # his tables carry pred_on_empty_ct; it is a property of the model, not of a volume,
    # so measure it once instead of 892 times.
    zero = np.zeros((320, 320, 320), dtype=np.uint8)
    pz, _, _ = RS.predict_prob(zero, pred_obj, ip)
    print(f"MARKER pred_on_empty_ct={float((pz > THRESH).mean()):.8f} "
          f"max_p={float(pz.max()):.6f}", flush=True)
    del zero, pz

    histp = outp.parent / "pred_histograms.jsonl"
    hh = histp.open("a")
    fh = outp.open("a")
    t0 = time.time()
    n_done = 0
    pending: deque = deque()
    max_pending = max(2, a.workers * 2)

    def drain(target: int) -> None:
        nonlocal n_done
        while len(pending) > target:
            fut, nm = pending.popleft()
            try:
                r = fut.result()
            except Exception as e:  # keep the run alive, record the failure
                r = {"sample": nm, "status": "error", "error": repr(e)}
            h = r.pop("_hist", None)
            if h is not None:
                hh.write(json.dumps(h) + "\n")
                hh.flush()
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            n_done += 1
            el = time.time() - t0
            rate = el / max(n_done, 1)
            eta = rate * (len(todo) - n_done)
            s1 = (r.get("shells") or {}).get("shell_1")
            print(f"MARKER [{n_done}/{len(todo)}] {nm} status={r.get('status')} "
                  f"enr={r.get('enrichment')} shell1={s1} "
                  f"n_sheet={(r.get('desc') or {}).get('n_sheet')} "
                  f"cpu={r.get('cpu_seconds')}s elapsed={el:.0f}s "
                  f"rate={rate:.2f}s/vol eta={eta/60:.0f}min", flush=True)

    for nm in todo:
        ct = np.asarray(tifffile.imread(str(IM / f"{nm}.tif")))
        tg = time.time()
        p, lo, hi = RS.predict_prob(ct, pred_obj, ip)
        gpu_s = time.time() - tg
        del ct
        pending.append((pool.submit(score_one, nm, p, lo, hi), nm))
        print(f"MARKER gpu {nm} {gpu_s:.1f}s queued={len(pending)}", flush=True)
        drain(max_pending - 1)

    drain(0)
    fh.close()
    hh.close()
    pool.shutdown(wait=True)
    print(f"MARKER done  {n_done} rows in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
