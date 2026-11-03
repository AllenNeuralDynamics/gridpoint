# stpt_register

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
| Deformable — discrete symmetric B-spline MRF | not yet implemented |

Requesting `--stages ... ,deformable` currently exits with an error rather than silently
producing a linear-only result.

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
stpt_register -f fixed.nii.gz -m moving.nii.gz -o out_dir [options]
```

Run `stpt_register --help` for the full list. Common options:

| Option | Meaning |
|---|---|
| `--stages rigid,affine` | which stages to run |
| `--fixed-mask` / `--moving-mask` | metric masks |
| `--threads N` | thread cap |
| `--levels N` | multi-resolution levels (default 4 ⇒ shrink 8,4,2,1) |
| `--bins`, `--linear-samples`, `--linear-iterations` | mutual-information knobs |
| `--translation-scale`, `--max-step`, `--min-step` | optimizer scaling (auto by default) |
| `--upstream-optimizer` | use the original Allen constants verbatim |
| `--warped-ext .nii` | warped-image extension (default `.nii.gz`) |

Outputs, written into `out_dir`:

| File | Contents |
|---|---|
| `stpt_Affine.txt` | collapsed rigid+affine transform (ITK format, maps fixed → moving) |
| `stpt_Warped<ext>` | moving image resampled onto the fixed grid |
| `stpt.log` | settings and final parameters |

## Validation

On synthetic pairs with a known translation, the recovered transform matches exactly and
the warped image is identical to the fixed image (normalized cross-correlation 1.0000):

| Case | Recovered translation | Expected | NCC after |
|---|---|---|---|
| 48³, shift 5 vox | `+5.000, 0, 0` | `+5, 0, 0` | 1.0000 |
| 96³, shift 7 vox | `+7.000, 0, 0` | `+7, 0, 0` | 1.0000 |
| 128³, shift 10 vox | `+10.000, 0, 0` | `+10, 0, 0` | 1.0000 |
| 96³, spacing 0.025 | `+0.175, 0, 0` | `+0.175, 0, 0` | 1.0000 |

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

Put `stpt_register` on `PATH` or set `STPTPATH` to its directory. commandants only
invokes the binary, so that MIT-licensed package stays independent of this one's license.
