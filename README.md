# gridpoint

Register **two arbitrary 3D volumes** with the Allen STPT/CCF registration algorithm —
the method used to build the Allen Mouse Common Coordinate Framework average template.

A single executable that depends on **ITK alone**. The upstream
[AllenInstitute/stpt_registration](https://github.com/AllenInstitute/stpt_registration)
code is a LIMS pipeline that needs Kakadu (JPEG2000), TinyXML, Boost, VTK, an Allen
image-series XML and a structured `model_directory` — and cannot take a 3D volume at all.
This keeps the algorithm and drops the pipeline. See [NOTICE.md](NOTICE.md) for exactly
what derives from upstream and what changed.

> **License:** Allen Institute Software License = 2-clause BSD **plus a non-commercial
> clause**. Commercial use requires the Allen Institute's written permission. See
> [LICENSE](LICENSE).

## Status

| Stage | State |
|---|---|
| Rigid (versor) — Mattes mutual information | **working** |
| Affine — Mattes mutual information | **working** |
| Deformable — discrete symmetric MRF, cross-correlation | **working** |

All three stages run. The deformable stage is the distinctive one: displacements live on
a control grid and are chosen by discrete optimization (minimum-spanning-tree message
passing over a label space) rather than continuous gradient descent, and it is run
forward and backward each level and composed, giving a symmetric field.

## Build

Needs ITK 5 (any build with the default modules) and a C++17 compiler.

```bash
cmake -S . -B build -DITK_DIR=/path/to/ITK-build
cmake --build build --config Release
```

On Windows with Visual Studio:

```bash
cmake -S . -B build -G "Visual Studio 16 2019" -A x64 -DITK_DIR=C:/path/to/ITK-build
cmake --build build --config Release
```

Only the ITK modules actually used are requested, so an ITK build missing optional image
format backends (JPEG2000 and friends) still works.

## Usage

```
gridpoint -f fixed.nii.gz -m moving.nii.gz -o out_dir [options]
```

Run `gridpoint --help` for the full list. Common options:

| Option | Meaning |
|---|---|
| `--stages rigid,affine` | which stages to run |
| `--fixed-mask` / `--moving-mask` | metric masks |
| `--threads N` | thread cap |
| `--levels N` | multi-resolution levels (default 4 ⇒ shrink 8,4,2,1) |
| `--bins`, `--linear-samples`, `--linear-iterations` | linear mutual-information knobs |
| `--grid 30,15,8,4` | deformable control-grid schedule, voxels (coarsest first) |
| `--lambda 1e-5` | deformable smoothness weight |
| `--samples 1024` | deformable random samples per node |
| `--write-deformation-field` | also write the dense displacement field |
| `--translation-scale`, `--max-step`, `--min-step` | optimizer scaling (auto by default) |
| `--upstream-optimizer` | use the original Allen constants verbatim |
| `--warped-ext .nii` | warped-image extension (default `.nii.gz`) |

Outputs, written into `out_dir`:

| File | Contents |
|---|---|
| `stpt_Affine.txt` | collapsed rigid+affine transform (ITK format, maps fixed → moving) |
| `stpt_Warped<ext>` | moving image resampled onto the fixed grid |
| `stpt_Warp.nii.gz` | dense deformation field — **only with `--write-deformation-field`** |
| `stpt.log` | settings and final parameters |

The deformation field is a 3-component float vector image on the fixed grid, in physical
units (voxel displacements are converted using the fixed image's spacing and direction
cosines). That makes it an ordinary ITK displacement field, so ANTs can consume it:

```
antsApplyTransforms -d 3 -i moving.nii.gz -r fixed.nii.gz -o redone.nii.gz \
    -t out_dir/stpt_Warp.nii.gz -t out_dir/stpt_Affine.txt
```

Order matters — `antsApplyTransforms` applies `-t` right to left, so the affine must be
listed *last*. Verified on a 96³ textured pair: that ordering reproduces gridpoint's own
`stpt_Warped.nii.gz` to NCC 0.9990 (the residual is interpolation, not convention),
while the reverse order gives 0.8992.

## Validation

**Linear stages.** On synthetic pairs with a known translation the recovered transform
matches exactly and the warped image is identical to the fixed image (normalized
cross-correlation 1.0000):

| Case | Recovered translation | Expected | NCC after |
|---|---|---|---|
| 48³, shift 5 vox | `+5.000, 0, 0` | `+5, 0, 0` | 1.0000 |
| 96³, shift 7 vox | `+7.000, 0, 0` | `+7, 0, 0` | 1.0000 |
| 128³, shift 10 vox | `+10.000, 0, 0` | `+10, 0, 0` | 1.0000 |
| 96³, spacing 0.025 | `+0.175, 0, 0` | `+0.175, 0, 0` | 1.0000 |

**Deformable stage.** On a textured 96³ phantom warped by a known sinusoidal field, the
deformable stage adds exactly what it should: little when the deformation is gentle
enough for an affine to absorb, a lot when it is local.

| Applied warp | before | linear only | linear + deformable | gain |
|---|---|---|---|---|
| amplitude 4 vox, 1 cycle | 0.8395 | 0.9862 | 0.9882 | +0.0019 |
| amplitude 8 vox, 2 cycles | 0.7260 | 0.8021 | **0.8577** | **+0.0556** |
| amplitude 6 vox, 3 cycles | 0.7888 | 0.7734 | **0.8556** | **+0.0822** |
| amplitude 10 vox, 3 cycles | 0.6492 | 0.6766 | **0.7193** | **+0.0427** |

A caveat worth knowing: the similarity term is an *uncentered* cross-correlation over
sampled voxels, so it is nearly flat inside large uniform regions and label choice there
is arbitrary. On a phantom of constant-valued blocks the deformable stage can therefore
make an already-good alignment slightly worse. Real anatomy is textured everywhere and
does not show this.

## Using it from Python

[`commandants`](https://github.com/AllenNeuralDynamics/commandants) wraps this binary as
its `stpt` backend, giving the same interface it provides for ANTs and elastix:

```python
from commandants import stpt

reg = stpt.presets.affine("fixed.nii.gz", "moving.nii.gz", "out_dir")
print(reg.to_shell())
result = reg.run()
warped = result.load("warped")     # -> SimpleITK.Image
```

Put `gridpoint` on `PATH` or set `STPTPATH` to its directory. commandants only
invokes the binary, so that MIT-licensed package stays independent of this one's license.
