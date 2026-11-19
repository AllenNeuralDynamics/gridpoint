/*=========================================================================
 *
 *  StptDeformableCuda.h
 *
 *  Optional CUDA backend for the deformable stage's cost volume -- the
 *  optimization() loop, which is ~99% of the stage's arithmetic.
 *
 *  Copyright (c) Allen Institute for Brain Science. All rights reserved.
 *  Licensed under the Allen Institute Software License; see LICENSE.
 *
 *  THE POINT OF THIS BACKEND IS THAT IT IS BIT-IDENTICAL TO THE CPU.
 *
 *  That is not an aspiration, it is a structural property of the algorithm and the
 *  reason a GPU port is worth doing here at all. Every label owns a private
 *  accumulator and the sum over samples is the outer loop, so assigning one label to
 *  one GPU thread performs exactly the same sequence of float operations in exactly
 *  the same order as the scalar CPU loop. Nothing is reassociated, because nothing
 *  is summed across threads.
 *
 *  What that buys: a GPU build still reproduces every template built with a CPU
 *  build, --device can be switched mid-project without invalidating anything, and
 *  the backend is validated with the same byte-for-byte output comparison as the CPU
 *  path rather than with a tolerance.
 *
 *  What it requires, and these are not optional:
 *    - `-fmad=false`. nvcc contracts a*b+c into an FMA by default, which rounds once
 *      instead of twice. This is the same trap as MSVC's /fp:precise, which was
 *      measured to break 9 of 23 outputs on the CPU side.
 *    - `-prec-div=true` (nvcc's default; -use_fast_math would disable it). The node
 *      and sample index decomposition divides by a float grid spacing, so division
 *      must be correctly rounded.
 *    - The final energy expression evaluated in double, as the CPU does.
 *
 *  Both are set in CMakeLists.txt. If you are tempted to add -use_fast_math, that
 *  single flag silently discards every guarantee above.
 *
 *=========================================================================*/
#ifndef StptDeformableCuda_h
#define StptDeformableCuda_h

#include <cstdint>
#include <string>

namespace stpt
{

/** Which backend the deformable cost volume should run on. */
enum class Device
{
  Auto, // GPU when this build has CUDA and a usable device is present, else CPU
  Cpu,  // always the CPU path
  Gpu   // GPU required; startup fails if there is not one
};

namespace cuda
{

/** True when this binary was compiled with the CUDA backend. Independent of whether
 *  a device is present -- a CUDA build on a machine with no GPU still returns true. */
bool
BuiltWithCuda();

/** True when a usable CUDA device is present and initialised right now. Safe to call
 *  on a build without CUDA (returns false) and on a machine with no driver. */
bool
DeviceAvailable();

/** One line naming the device, e.g. "NVIDIA GeForce RTX 3080 (10.0 GiB, cc 8.6)".
 *  Empty when no device is available. */
std::string
DeviceDescription();

/** Why no device is available, phrased for an error message -- "this build has no
 *  CUDA support", "no CUDA driver found", and so on. Empty when one IS available. */
std::string
UnavailableReason();

/** Everything optimization() needs to hand the cost volume to the GPU.
 *
 *  Pointers are HOST pointers; the backend owns all device allocation and transfer.
 *  Names and semantics mirror the locals in optimization() one for one, deliberately:
 *  the kernel has to reproduce that arithmetic exactly, and anything renamed in
 *  transit is somewhere a difference can hide. */
struct OptimizationInputs
{
  // Fixed image, indexed fixed[yy + xx*m + zz*m*n].
  const float * fixed{ nullptr };
  std::int64_t  m{ 0 };
  std::int64_t  n{ 0 };
  std::int64_t  o{ 0 };

  // Moving image as optimization() sees it: the caller's buffer, or the 2x
  // interpolated one in subpixel mode. Indexed movingi[y + x*mi + z*mi*ni].
  const float * movingi{ nullptr };
  std::int64_t  mi{ 0 };
  std::int64_t  ni{ 0 };
  std::int64_t  oi{ 0 };

  // Control grid.
  int          gridx{ 0 };
  int          gridy{ 0 };
  int          gridz{ 0 };
  float        gridstep{ 0.0f }; // float on purpose: the index arithmetic divides by it
  std::int64_t sz1{ 0 };         // node count = gridx*gridy*gridz

  // Label space. Offsets are per-axis arrays of `len` entries, already carrying the
  // subpixel x2 scaling where it applies.
  int         len{ 0 };
  int         lengthfull{ 0 }; // len^3
  const int * yOff{ nullptr };
  const int * xOff{ nullptr };
  const int * zOff{ nullptr };

  // Sampling.
  int         maxsamp{ 0 };
  bool        randommode{ false };
  const int * sampleIdx{ nullptr }; // sz1*maxsamp entries when randommode, else null

  // Boundary test inputs.
  bool subpixel{ false };
  int  hw_local{ 0 };

  // Per-node mask; null means every node is active.
  const unsigned char * activeNodes{ nullptr };

  float epsilon{ 0.0f };

  // Output: sz1*lengthfull floats, label-major within each node.
  float * energy{ nullptr };
};

/** Compute the cost volume on the GPU.
 *
 *  Returns false when this particular call cannot be served (no device, or the
 *  problem shape does not fit the kernel's assumptions) -- the caller must then run
 *  the CPU path, which is always correct. Returns true only when `energy` has been
 *  fully written.
 *
 *  Throws std::runtime_error on an actual CUDA failure, which is a real error rather
 *  than a reason to silently produce a slower correct answer. */
bool
RunOptimization(const OptimizationInputs & in);

} // namespace cuda
} // namespace stpt

#endif
