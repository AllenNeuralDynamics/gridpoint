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

ITK 5 and a C++17 compiler. That is the whole dependency list — no Kakadu/JPEG2000,
TinyXML, Boost or VTK, which the original code required.

### Linux (easiest: ITK from apt)

Ubuntu 22.04 and 24.04 package ITK 5 directly. Earlier releases ship ITK 4 and will not
work.

```bash
sudo apt-get install -y build-essential cmake libinsighttoolkit5-dev
```

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j"$(nproc)"
```

No `ITK_DIR` is needed when ITK is installed system-wide. The binary lands at
`build/gridpoint`; `sudo cmake --install build` puts it in `/usr/local/bin`.

### From your own ITK build tree

Point `ITK_DIR` at the directory containing `ITKConfig.cmake`:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DITK_DIR=/path/to/ITK-build
```

### Windows (Visual Studio)

```bash
cmake -S . -B build -G "Visual Studio 16 2019" -A x64 -DITK_DIR=C:/path/to/ITK-build
```

```bash
cmake --build build --config Release
```

### Verify

```bash
./build/gridpoint --help
```

### Notes and troubleshooting

- Only the ITK modules actually used are requested, so an ITK build missing optional
  image-format backends (JPEG2000 and friends) still links.
- The project enables both C and C++ even though no C source exists here. ITK's HDF5
  module reaches CMake's `FindHDF5`, which validates the toolchain by compiling a `.c`
  probe; a C++-only project fails that probe with *"Unknown extension .c"*.
- Configuring against Ubuntu's ITK prints a wall of warnings about missing `vtkgdcm*`
  and `gdcm*` files. That is a known defect in the `libgdcm` package, whose
  `GDCMTargets.cmake` references files the binary package does not install. It is noise;
  the build succeeds.

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
| `stpt_Warp.nii.gz` | forward deformation field, fixed → moving — **only with `--write-deformation-field`** |
| `stpt_InverseWarp.nii.gz` | backward field, moving → fixed — same flag |
| `stpt.log` | settings and final parameters |

Each deformation field is a 3-component float vector image on the fixed grid, in physical
units (voxel displacements are converted using the fixed image's spacing and direction
cosines). That makes them ordinary ITK displacement fields, so ANTs can consume them:

```
antsApplyTransforms -d 3 -i moving.nii.gz -r fixed.nii.gz -o redone.nii.gz \
    -t out_dir/stpt_Warp.nii.gz -t out_dir/stpt_Affine.txt
```

Order matters — `antsApplyTransforms` applies `-t` right to left, so the affine must be
listed *last*. Verified on a 96³ textured pair: that ordering reproduces gridpoint's own
`stpt_Warped.nii.gz` to NCC 0.9990 (the residual is interpolation, not convention),
while the reverse order gives 0.8992.

### Going the other way

The deformable stage is symmetric: it optimizes forward and backward at every level and
composes them, so the backward field is produced alongside the forward one rather than
recovered by inverting it afterwards. That is what `stpt_InverseWarp.nii.gz` holds, and
it is what you need to carry fixed-space data — atlas labels, annotations, point sets —
onto the moving image:

```
antsApplyTransforms -d 3 -i ccf_annotation.nii.gz -r moving.nii.gz -o labels_in_moving.nii.gz \
    -n NearestNeighbor -t [out_dir/stpt_Affine.txt,1] -t out_dir/stpt_InverseWarp.nii.gz
```

Note the order is the mirror of the forward recipe — the inverted affine comes *first* on
the command line now — and that `,1` is what inverts it.

Both fields are written on the fixed grid, because the deformable stage runs after linear
alignment has already resampled the moving image onto that grid. ANTs interpolates a
displacement field wherever it needs it, so using the inverse field with `-r moving.nii.gz`
works; it simply has nothing to say outside the fixed image's bounding box.

Measured on the same 96³ pair, composing the two fields returns points to where they
started — median residual **0.006 voxels**, p95 **0.26**, against forward displacements
averaging 1.76 voxels (p95 4.82). End to end, mapping fixed into moving space:

| Transform chain | NCC vs. the real moving image |
|---|---|
| `[Affine,1]` then `InverseWarp` | **0.8931** |
| `InverseWarp` then `[Affine,1]` (wrong order) | 0.8225 |
| `[Affine,1]` alone | 0.8197 |
| `[Affine,1]` with the *forward* field | 0.7416 |

The inverse direction scores about what the forward direction does (0.8889), which is the
symmetry showing up where it should. Both wrong variants land at or below the affine-only
baseline, so a convention mistake here is loud rather than subtle.

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
