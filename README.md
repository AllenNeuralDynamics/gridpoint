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
| `--fixed-mask` / `--moving-mask` | metric masks, see [Masking](#masking) |
| `--threads N` | thread cap |
| `--levels N` | multi-resolution levels (default 4 ⇒ shrink 8,4,2,1) |
| `--bins`, `--linear-samples`, `--linear-iterations` | linear mutual-information knobs (linear stages only) |
| `--grid 30,15,8,4` | deformable control-grid schedule, voxels (coarsest first) |
| `--lambda 1e-5` | deformable smoothness weight |
| `--samples 1024` | deformable random samples per node |
| `--write-deformation-field` | also write the displacement fields |
| `--warp-resolution grid\|full` | which grid to write them on, see [Deformation fields](#deformation-fields) |
| `--translation-scale`, `--max-step`, `--min-step` | optimizer scaling (auto by default) |
| `--upstream-optimizer` | use the original Allen constants verbatim |
| `--warped-ext .nii` | warped-image extension (default `.nii.gz`) |

## Masking

`--fixed-mask` and `--moving-mask` take `unsigned char` images, non-zero meaning inside,
interpreted in physical space -- they need correct origin/spacing/direction but do not
have to share the image's grid. They restrict **both** the linear and the deformable
stages, but they do so in different ways, and the difference matters.

**Linear stages.** The masks are handed to the Mattes mutual-information metric
(`SetFixedImageMask` / `SetMovingImageMask`), restricting where it is evaluated. This is
the same concept as the ANTs `-x [fixedMask,movingMask]` option.

**Deformable stage.** Masking is applied to the *control grid*. A node whose cell holds
no in-mask voxel is given a flat zero data cost, so the regularizer alone places it.

Two consequences worth understanding:

- **Masked-out regions are not held still.** Their displacement is a smooth
  extrapolation from the masked-in neighbourhood. That is usually what is wanted, but
  it is not zero, and it is not "the mask is a hard boundary".
- **It is per node, not per sample.** A cell that straddles the mask edge is included
  whole, so the effective boundary is quantized to the control-grid spacing of each
  level -- 30 voxels at the coarsest default level, 4 at the finest. Sample-level
  rejection would be sharper, but would leave different candidate displacements with
  different numbers of valid samples, biasing the correlation between them.

The moving mask is carried through the linear transform before use, because the
deformable stage sees the moving image only after that transform has resampled it onto
the fixed grid. It is then re-displaced by the deformation field at every pyramid level,
so the mask follows the image as the alignment improves -- and in the backward half of
the symmetric optimization the roles swap, the fixed mask being displaced by the inverse
field.

A masked run is also **faster**: a skipped node skips the sampling loop, which is the
most expensive step in the stage. Expect the saving to concentrate at the fine levels;
a coarse cell is 30^3 voxels and usually clips the mask somewhere.

Run with `-v` to see how much of the grid is actually live:

```
  deformable: fixed mask covers 26.8% of the volume
    active control nodes: forward 2304/4096, backward 2048/4096
```

Note that `--bins` has no effect on the deformable stage; it is a linear-stage knob
only. The deformable stage uses normalized cross-correlation, not mutual information.

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

### Deformation fields

`--warp-resolution` chooses the grid the two fields are written on. **`grid` is the
default.**

Every level of the deformable stage accumulates at *control-grid* resolution -- `genDfm`,
`smoothDfm` and `symmetricMapping` all operate on the control array -- and the dense
field is produced by a single trilinear upsample at the end of the level. So the control
grid is not a reduction of the field, it *is* the field; writing one vector per voxel
writes every number `spacing^3` times over. At the default schedule the finest spacing is
4, so that is 64x.

| | 192³ phantom | 10 µm CCF (1320×800×1140) |
|---|---|---|
| `full` | 58 MB | ~12 GB |
| `grid` | 0.87 MB | ~195 MB |

ITK and ANTs read a displacement field through its origin, spacing and direction cosines
and have never required it to share the reference image's grid, so the smaller file is a
drop-in -- the `antsApplyTransforms` commands above are unchanged.

Verified on a 192³ pair with a local sinusoidal warp:

* the written grid field equals a stride-4 decimation of the dense field exactly
  (`max|diff| = 0`), which is what makes it lossless rather than merely small;
* `antsApplyTransforms` produces the same warped image from either, to `1.4e-06` on
  intensities spanning 0 to 2.57 -- float rounding, not a difference in the transform;
* `--warp-resolution full` is bit-identical to the pre-0.3.0 build at 1, 4 and 24 threads.

The grid field carries one extra sample per axis, replicating the last. The control
grid's own extent stops short of the fixed grid's by up to `spacing - 1` voxels, and ANTs
pads with zero outside a field rather than clamping to its edge; without that sample the
outer shell loses its displacement. The spacing is derived from the grid dimensions
rather than assumed equal to the control spacing, because the core computes
`nodes = size / spacing` with integer division -- when an axis length is not a multiple
of the spacing, the true ratio is not the spacing.

Use `full` when something downstream wants the dense array itself rather than a transform
ANTs will resample. Note that anything which makes ANTs *write* a field --
`antsApplyTransforms -o [composite.nii.gz,1]`, or a subsequent `antsRegistration` --
emits it on the reference grid, so the compact form does not propagate.

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

## Memory

The deformable stage holds the whole volume in RAM -- the discrete optimization is a
global minimum-spanning-tree message pass over the control grid, so there is nothing to
tile or stream. Peak usage is therefore set almost entirely by voxel count:

```
peak  ~=  85 bytes/voxel  +  the largest single level's cost volumes  +  ~20 MiB
```

Measured against peak working set, with the default grid schedule:

| Volume | Voxels | Measured peak | Before |
|---|---|---|---|
| 64^3 | 262,144 | 40 MiB | 64 MiB |
| 96^3 | 884,736 | 121 MiB | 202 MiB |
| 128^3 | 2,097,152 | 226 MiB | 418 MiB |
| 160^3 | 4,096,000 | 433 MiB | 808 MiB |
| 192^3 | 7,077,888 | 744 MiB | 1.36 GiB |

The size-linear part of that fell from 201 to 106 bytes per voxel, with
bit-identical output. Two things went:

* **Sampling-coordinate arrays.** `interpolate3d` was handed three full-size float
  arrays of where to sample. At the finest level that interpolation runs on a grid
  twice as fine on every axis, so the three arrays are 96 bytes per voxel of the input
  -- the single largest live allocation in the stage, and measurement put the peak
  inside it. Every entry was the output index over a constant, so they are now derived
  in the loop instead of stored.
* **Mutual-information scaffolding.** The deformable stage scores with NCC (see
  the metric note below); `--bins` only reaches the linear stages. It nonetheless built a
  histogram-bin index per voxel for each image and a private copy of each image, and
  read none of it. That was 16 bytes per voxel of allocation, though only about half of
  it was ever written and so only half ever became resident.

The cost-volume term is the only part the parameters move: each level allocates
`nodes x labels x 4` bytes twice (forward and backward). Those two factors nearly cancel
across the schedule -- coarse levels have few nodes over a huge label space, fine levels
the reverse -- so the per-level figures stay within about 2x of each other and peak at
the *finest* level.

What that means at atlas scale:

| Case | Shape | Voxels | Estimated peak |
|---|---|---|---|
| CCF 100 um | 132x80x114 | 1,203,840 | 256 MiB |
| CCF 50 um | 264x160x228 | 9,630,720 | 1.9 GiB |
| CCF 25 um | 528x320x456 | 77,045,760 | 14.8 GiB |
| STPT 20 um | 600x400x500 | 120,000,000 | 23.0 GiB |
| CCF 10 um | 1320x800x1140 | 1,203,840,000 | 230.6 GiB |

gridpoint prints this estimate before the deformable stage starts, and warns if it
exceeds the machine's physical memory. The estimate uses the schedule you asked for, not
the capped one, so on small images where the level or grid caps fire it reads slightly
high.

### Fitting a large registration

Memory scales with voxel count, so halving the resolution costs about **8x less**.
Register coarse and apply the result at full resolution:

```
gridpoint -f fixed_25um.nii.gz -m moving_25um.nii.gz -o out --write-deformation-field
antsApplyTransforms -d 3 -i moving_10um.nii.gz -r fixed_10um.nii.gz -o warped_10um.nii.gz \
    -t out/stpt_Warp.nii.gz -t out/stpt_Affine.txt
```

**`--grid` is specified in voxels, not microns.** If you change resolution, scale the
schedule by the same factor or you are optimizing a physically different deformation --
the default `30,15,8,4` means 750/375/200/100 um at 25 um, but 300/150/80/40 um at 10 um.

Little is lost doing this. The finest control grid is 4 voxels, so the field was never
going to carry detail below that; registering finer buys a denser sampling of the metric,
not a finer warp. On a 96^3 textured phantom, halving the resolution and applying the
result at full resolution scored NCC 0.9110 against 0.8862 for registering at full
resolution directly -- though that phantom's warp is a smooth global sinusoid, which
downsampling preserves almost perfectly, so treat it as evidence that the loss is small
rather than that it is negative.

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
