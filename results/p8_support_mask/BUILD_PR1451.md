# Building the PR #1451 tracer from an extracted container image

The VC3D binaries here come from the published container image, unpacked to a
rootfs rather than run under docker. Its CMake build tree hard-codes the
container's `/src`, which does not exist on the host and cannot be symlinked
without root, so an in-place rebuild looks impossible at first. `bwrap`
reproduces the path without privileges:

    R=/path/to/vc3d-rootfs
    bwrap --bind $R / --dev-bind /dev /dev --proc /proc \
          --bind /mnt/vesuvius /mnt/vesuvius --bind /tmp /tmp \
          --setenv PATH /usr/local/bin:/usr/bin:/bin --setenv HOME /root \
          /usr/bin/ninja -C /src/build vc_grow_seg_from_seed

Four of the five hunks in the PR applied to this tree; the first needed hand
placement because the image predates the upstream context around
`volume.chunkedCache()`. The rebuild is incremental and finishes with zero
errors; the resulting binary is at `/src/build/bin/vc_grow_seg_from_seed` and
carries the new `volume support test enabled` string, while the image's own
`/usr/local/bin` copy is left untouched so existing runs are unaffected.
