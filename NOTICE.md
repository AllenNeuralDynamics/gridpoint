# Provenance and attribution

This project is a **derived work** of
[AllenInstitute/stpt_registration](https://github.com/AllenInstitute/stpt_registration),
copyright © Allen Institute for Brain Science, distributed under the Allen Institute
Software License (2-clause BSD **plus a non-commercial clause**). That license applies
to this repository in full — see [LICENSE](LICENSE). In particular, **redistribution and
use for commercial purposes require the Allen Institute's written permission**
(terms@alleninstitute.org).

## What derives from what

| File here | Derived from (upstream, under `lims2_modules/common/register/`) |
|---|---|
| `src/StptLinearRegistration.{h,cxx}` | `Registration/idpMutualInfoVersorRigidVolumeRegistration.{h,txx}`, `Registration/idpMutualInfoCenteredAffineRegularGradientDescentRegistration.{h,txx}` |
| `src/gridpoint.cxx` | New. Replaces the `idpProjection*Module` entry points with a two-volume CLI. |

## What changed, and why

The upstream code is a LIMS pipeline: it reads an Allen image-series XML plus a
structured `model_directory`, and links Kakadu (JPEG2000), TinyXML, Boost and VTK. It
also **cannot register two 3D volumes** — `ImageSeriesUtilities::PopulateVolume`
assembles its volume slice-by-slice from one 2D file per `<sub-image>`, pinning the
z-extraction to slice 0.

This project keeps the *algorithm* and discards the pipeline:

- **ITK is the only dependency.** No Kakadu/JPEG2000, TinyXML, Boost or VTK; no XML and
  no `model_directory`. Inputs are two ordinary 3D volumes in any format ITK reads.
- **Ported from ITK 3.20 to ITK 5.** The v3-style registration classes used upstream
  (`MattesMutualInformationImageToImageMetric`, `MultiResolutionImageRegistrationMethod`,
  `VersorRigid3DTransformOptimizer`, `RegularStepGradientDescentOptimizer`,
  `CenteredVersorTransformInitializer`) are still shipped by ITK 5, so the numerics carry
  over directly.
- **Affine transform.** Upstream names its class `CenteredAffine` but scales only
  parameters 0..11, which matches ITK's 12-parameter `AffineTransform`; that is what is
  used here.

### Two deliberate deviations, both evidence-backed

1. **Optimizer scaling is derived from the image by default.** Upstream hard-codes
   `TranslationScale = 1e7` with `MaximumStepLength = 60`, values calibrated to the Allen
   data's physical units. On synthetic 48³/96³/128³ pairs those constants diverge on two
   of the three sizes (the optimizer leaves the overlap entirely, and the metric aborts
   with "too many samples map outside moving image buffer" or "joint PDF summed to zero").
   By default the translation scale is therefore set to the fixed image's physical extent,
   with step lengths 0.2 → 1e-4; this recovers the exact transform at every size tested,
   in both millimetre- and micron-spaced data. Pass `--upstream-optimizer` for the
   original constants, or set `--translation-scale` / `--max-step` / `--min-step`
   explicitly.

2. **The multi-resolution schedule is capped to the image size.** The upstream schedule
   (4 levels, coarsest shrink 8) assumes a full-size brain volume. On a small image the
   coarsest level can collapse to a handful of voxels and make the Mattes metric
   degenerate, so the level count is reduced to keep at least 16 voxels on the shortest
   axis. A warning is printed whenever this happens.

3. **The deformable control-grid spacing is capped to the image size**, for the same
   reason as the pyramid: the upstream 30-voxel coarsest grid assumes a full-size volume,
   and on a small image it leaves a 1-node grid with nothing to optimize. A warning is
   printed when it applies.

Deviations 1-3 only affect inputs on which the upstream constants fail outright.

### Bugs fixed while porting

These are defects in the original, not consequences of the port:

- **`delete x;` on array allocations** (11 sites) is undefined behaviour; changed to
  `delete[] x;`.
- **The random sample index could run one past the end.**
  `i1 = (int)(rand() * pow(gridstep,3) / float(RAND_MAX))` reaches `gridstep^3` exactly
  when `rand()` returns `RAND_MAX`, and the result is used unguarded to offset into the
  image. On glibc (`RAND_MAX` = 2^31-1) that is a 1-in-2-billion event and effectively
  invisible; with the Microsoft CRT (`RAND_MAX` = 32767) it is near-certain across the
  ~65k draws per level, and it access-violated immediately. The index is now clamped.
- **A progress divisor could be zero.** `frac = sz1 / 25` is 0 whenever there are fewer
  than 25 grid nodes, and `i % frac` then divides by zero. Guarded at both sites.

### Observation

`fixedInd` / `movingiInd` and the `NumBin` / `MaxFixed` / `MaxMoving` settings are
allocated, filled and freed but never read: the deformable similarity term is an
uncentered cross-correlation, not a histogram metric. They appear to be vestigial from an
earlier mutual-information implementation. They are kept here so the port stays faithful,
but they do not affect the result.
