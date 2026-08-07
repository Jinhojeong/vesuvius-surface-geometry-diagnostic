"""run_gpu_roi.py equivalent: the entry TAUIL-Abd-Elilah's scripts use to call vesuvius.predict.

Taken from his run_gpu_roi_cached.py @ 9afa412 with the chunk-cache half deleted, which is what
that file describes itself as ("run_gpu_roi.py plus an on-disk chunk cache ... with the variable
unset it is byte-for-byte the behaviour of run_gpu_roi.py").

Two of the three things below are his and change nothing. The allocator cap is environment work,
and the torch.compiler.disable shim is a no-op on torch 2.6, which accepts `reason`.

The third is mine and is opt-in via VESUVIUS_ARCH_SHIM=1. Without it, villa 0.2.4 cannot load an
nnU-Net checkpoint at all against nnunetv2 2.8.1, and dies before any data is read:

  * load_nnunet_model.load_model first tries `trainer_class.build_network_architecture(
    arch_class_name, arch_kwargs, req_import, n_in, n_out, enable_deep_supervision=False)`.
    nnunetv2 2.8.1 changed that signature to (plans_manager, configuration_manager, n_in, n_out,
    enable_deep_supervision), so the call raises TypeError. villa swallows it.
  * the fallback, load_nnunet_model.initialize_network, does `exec(f"import {i}")` over
    _kw_requires_import. Those entries are KWARG NAMES ('conv_op', 'norm_op', ...), not module
    paths, so it raises ModuleNotFoundError: No module named 'conv_op'. nnU-Net's own code
    pydoc.locate()s the kwarg VALUES instead. This is the same in villa main as in 0.2.4.
  * Inferer.infer() catches the exception, prints, and returns None; main() then dies unpacking
    it. Exit code stays 0 with an empty output directory, exactly the failure mode his scripts
    warn about.

The shim replaces initialize_network with nnU-Net's own construction: pydoc.locate the class,
pydoc.locate the kwargs named in _kw_requires_import, instantiate. That is what an older
nnunetv2 would have handed villa through the trainer path, so it restores his environment rather
than changing villa's behaviour. It touches network construction ONLY. Nothing on the
normalization, dataloading, sliding-window or blending path is modified, and that path is the
thing under test.
"""
from __future__ import annotations
import os
import sys

os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")

import torch  # noqa: E402

if torch.cuda.is_available():
    torch.cuda.set_per_process_memory_fraction(
        float(os.environ.get("VESUVIUS_GPU_MEM_FRACTION", "0.88"))
    )

_orig_disable = torch.compiler.disable


def _disable(fn=None, *, recursive=True, reason=None):  # noqa: ARG001
    return _orig_disable(fn, recursive=recursive)


torch.compiler.disable = _disable

if os.environ.get("VESUVIUS_ARCH_SHIM") == "1":
    import pydoc
    import vesuvius.utils.models.load_nnunet_model as _lnm

    def _initialize_network(architecture_class_name, arch_init_kwargs,
                            arch_init_kwargs_req_import, num_input_channels,
                            num_output_channels, enable_deep_supervision=False):
        kw = dict(arch_init_kwargs)
        for name in arch_init_kwargs_req_import:
            if kw.get(name) is not None:
                kw[name] = pydoc.locate(kw[name])
        cls = pydoc.locate(architecture_class_name)
        if cls is None:
            raise RuntimeError(f"could not locate {architecture_class_name}")
        net = cls(input_channels=num_input_channels, num_classes=num_output_channels,
                  deep_supervision=enable_deep_supervision, **kw)
        print(f"[arch shim] built {architecture_class_name} via pydoc.locate", file=sys.stderr)
        return net

    _lnm.initialize_network = _initialize_network

from vesuvius.models.run import inference  # noqa: E402

if __name__ == "__main__":
    sys.argv = ["vesuvius.predict"] + sys.argv[1:]
    inference.main()
