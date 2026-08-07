"""Run villa's own `vesuvius.predict` wrapper on the m7 checkpoint and score it on our geometry.

Why this exists. Our claim is that `vesuvius.predict` normalizes an nnU-Net checkpoint with
per-volume instance z-score rather than the CTNormalization the checkpoint's plans.json declares.
We read the source and we ran the two normalizations inside our OWN predictor path
(norm_ablation.py). We had never run villa's wrapper itself. This does.

Path, mirroring TAUIL-Abd-Elilah's bench_m7_recall.py / m7_margin_fp.py @ 9afa412:

  * write the CT block to a uint8 zarr, chunks 128
  * `vesuvius.predict --model_path ... --input_dir ct.zarr --output_dir logits --device cuda
    --disable_tta --batch_size 1 --num_workers 2`, with NO --normalization, so the default applies
  * `blending.main()` over the logits directory into merged.zarr
  * p = sigmoid(logit1 - logit0), threshold 0.2, class 2 excluded, scored over the inner 128^3

Two geometry modes, because the installed vesuvius 0.2.4 has no --bbox argument (his version did):

  crop256  the zarr holds only the centred 256^3 block, so villa tiles exactly the voxels his
           `--bbox 32:288,32:288,32:288` asked for. Confirmed on the wire: villa reports patch
           origins {0,64}^3, the same 8 positions a 256-long bbox at step 96... i.e. patch 192,
           step 0.5, last position clamped to 64. This is the arm directly comparable to
           run_shells.py / norm_ablation.py, which also crop first.
  full320  the zarr holds the whole 320^3 volume and villa tiles all of it; the inner 128^3 is
           then sliced out. Different patch grid, kept as a control on the mode choice.

--prenorm ct writes a float32 zarr already CT-normalized with the plans constants and passes
`--normalization none`, which is the fidelity control: if villa's own machinery is faithful, that
arm has to land on our CTNormalization rows. It isolates normalization as the only free variable.

Both stdout and stderr of every child process are kept, so a normalization warning (or its
absence) is on the record.
"""
from __future__ import annotations
import argparse, json, os, shutil, subprocess, sys, time
from pathlib import Path

import numpy as np
import tifffile
import zarr

sys.path.insert(0, "/mnt/vesuvius/experiments/shell_split/gen")
from geom import margin_mask, distance_profile, analyse, descriptives, TRIM, SIZE, THRESH  # noqa: E402

PY = os.environ.get("VENV_PY", "/home/jinhojeong/Vesuvius/.venv/bin/python")
HERE = Path(__file__).resolve().parent
MODEL = "/mnt/vesuvius/models/surface_m7_nnunet"
IM = Path("/mnt/vesuvius/kaggle892/images")
LB = Path("/mnt/vesuvius/kaggle892/labels")

ENV = {**os.environ, "nnUNet_compile": "0", "TORCHDYNAMO_DISABLE": "1",
       "PYTHONIOENCODING": "utf-8"}


def run_wrapper(zpath: Path, outdir: Path, norm: str | None, extra: list[str]) -> subprocess.CompletedProcess:
    cmd = [PY, str(HERE / "villa_predict.py"),
           "--model_path", MODEL,
           "--input_dir", str(zpath),
           "--output_dir", str(outdir),
           "--device", "cuda", "--disable_tta",
           "--batch_size", "1", "--num_workers", "2"]
    if norm:
        cmd += ["--normalization", norm]
    cmd += extra
    return subprocess.run(cmd, env=ENV, capture_output=True, text=True, errors="replace")


def run_blend(logits: Path, merged: Path) -> subprocess.CompletedProcess:
    code = ("import sys\n"
            "from vesuvius.models.run import blending\n"
            f"sys.argv=['b', r'{logits}', r'{merged}']\n"
            "blending.main()\n")
    return subprocess.run([PY, "-c", code], env=ENV, capture_output=True, text=True, errors="replace")


def endpoints(p: np.ndarray, labc: np.ndarray) -> list:
    """Identical arithmetic to norm_ablation.endpoints, so W is on the same scale as P/Z."""
    sheet, scored = labc == 1, labc != 2
    m = p > THRESH
    tp = float((m & sheet).sum())
    npp = float((m & scored).sum())
    return [round(tp / float(sheet.sum()), 4), round(tp / max(npp, 1.0), 4),
            round(npp / float(scored.sum()), 4)]


def ct_normalize(x: np.ndarray, ip: dict) -> np.ndarray:
    """nnunetv2 CTNormalization, verbatim maths; same function as run_shells.ct_normalize."""
    x = x.astype(np.float32, copy=True)
    np.clip(x, ip["percentile_00_5"], ip["percentile_99_5"], out=x)
    x -= float(ip["mean"])
    x /= max(float(ip["std"]), 1e-8)
    return x


def one(nm: str, work: Path, mode: str, norm: str | None, logdir: Path,
        prenorm: str = "none") -> dict:
    t0 = time.time()
    ct = np.asarray(tifffile.imread(str(IM / f"{nm}.tif")))
    lab = np.asarray(tifffile.imread(str(LB / f"{nm}.tif")))
    off = (ct.shape[0] - SIZE) // 2
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    if mode == "crop256":
        blk = np.ascontiguousarray(ct[off:off + SIZE, off:off + SIZE, off:off + SIZE])
        lo_local, hi_local = TRIM, SIZE - TRIM
    else:
        blk = ct
        lo_local, hi_local = off + TRIM, off + SIZE - TRIM
    if prenorm == "ct":
        ip = json.loads((Path(MODEL) / "plans.json").read_text())[
            "foreground_intensity_properties_per_channel"]["0"]
        blk = ct_normalize(blk, ip)
        dt = "float32"
    else:
        dt = "uint8"
    zpath = work / "ct.zarr"
    z = zarr.open(str(zpath), mode="w", shape=blk.shape, chunks=(128, 128, 128), dtype=dt)
    z[:] = blk

    logits = work / "logits"
    r1 = run_wrapper(zpath, logits, norm, [])
    (logdir / f"{nm}.predict.log").write_text(
        "CMD-RC %s\n--- STDOUT ---\n%s\n--- STDERR ---\n%s\n" % (r1.returncode, r1.stdout, r1.stderr))
    if not logits.exists() or not any(logits.iterdir()):
        return {"sample": nm, "status": "predict_failed", "rc": r1.returncode,
                "tail": ((r1.stdout or "") + (r1.stderr or ""))[-1200:]}

    merged = work / "merged.zarr"
    r2 = run_blend(logits, merged)
    (logdir / f"{nm}.blend.log").write_text(
        "CMD-RC %s\n--- STDOUT ---\n%s\n--- STDERR ---\n%s\n" % (r2.returncode, r2.stdout, r2.stderr))
    if not merged.exists():
        return {"sample": nm, "status": "blend_failed", "rc": r2.returncode,
                "tail": ((r2.stdout or "") + (r2.stderr or ""))[-1200:]}

    a = zarr.open(str(merged), mode="r")
    arr = a if not hasattr(a, "keys") or hasattr(a, "shape") else a[list(a.keys())[0]]
    sl = slice(lo_local, hi_local)
    l0 = np.asarray(arr[0, sl, sl, sl]).astype(np.float32)
    l1 = np.asarray(arr[1, sl, sl, sl]).astype(np.float32)
    p = 1.0 / (1.0 + np.exp(-(l1 - l0)))

    gl = slice(off + TRIM, off + SIZE - TRIM)          # absolute coords for labels
    labc = lab[gl, gl, gl]
    mg = margin_mask(ct, lab)[gl, gl, gl]
    pred = p > THRESH
    out = analyse(labc, mg, pred)
    out["shells"] = distance_profile(labc, pred)
    out["desc"] = descriptives(labc, pred)
    out["endpoints"] = endpoints(p, labc)
    out["sample"] = nm
    out["mode"] = mode
    out["normalization_flag"] = norm or "DEFAULT(none passed)"
    out["prenorm"] = prenorm
    out["merged_shape"] = list(arr.shape)
    out["logit_stats"] = {"l0_mean": round(float(l0.mean()), 4), "l1_mean": round(float(l1.mean()), 4),
                          "p_mean": round(float(p.mean()), 4)}
    out["blk_mean"] = round(float(ct[off:off + SIZE, off:off + SIZE, off:off + SIZE].mean()), 2)
    out["blk_std"] = round(float(ct[off:off + SIZE, off:off + SIZE, off:off + SIZE].std()), 2)
    out["seconds"] = round(time.time() - t0, 1)
    out["status"] = out.get("status", "ok")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", default="crop256", choices=["crop256", "full320"])
    ap.add_argument("--norm", default=None, help="value for --normalization; omit for villa's default")
    ap.add_argument("--prenorm", default="none", choices=["none", "ct"],
                    help="'ct' writes an already-CT-normalized float32 zarr (pair with --norm none)")
    ap.add_argument("--work", default="/mnt/vesuvius/experiments/wrapper_check/work")
    a = ap.parse_args()

    names = json.loads(Path(a.names).read_text())
    outp = Path(a.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    logdir = outp.parent / (outp.stem + "_logs")
    logdir.mkdir(parents=True, exist_ok=True)
    done = set()
    if outp.exists():
        for ln in outp.read_text().splitlines():
            if ln.strip():
                done.add(json.loads(ln)["sample"])
    fh = outp.open("a")
    print(f"MARKER start n={len(names)} done={len(done)} mode={a.mode} norm={a.norm} "
          f"prenorm={a.prenorm}", flush=True)
    for k, nm in enumerate(names):
        if nm in done:
            continue
        try:
            r = one(nm, Path(a.work), a.mode, a.norm, logdir, a.prenorm)
        except Exception as exc:  # noqa: BLE001
            import traceback
            r = {"sample": nm, "status": "error", "error": f"{type(exc).__name__}: {exc}",
                 "tb": traceback.format_exc()[-1500:]}
        fh.write(json.dumps(r) + "\n")
        fh.flush()
        print(f"MARKER [{k+1}/{len(names)}] {nm} {r['status']} {r.get('endpoints')} "
              f"{r.get('seconds')}s", flush=True)
    fh.close()
    print("MARKER done", flush=True)


if __name__ == "__main__":
    main()
