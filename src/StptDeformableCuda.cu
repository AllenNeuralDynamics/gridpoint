/*=========================================================================
 *
 *  StptDeformableCuda.cu
 *
 *  CUDA backend for optimization(), the deformable stage's NCC cost volume.
 *  Compiled only when GRIDPOINT_ENABLE_CUDA is on; otherwise
 *  StptDeformableCudaStub.cxx supplies the same symbols.
 *
 *  Copyright (c) Allen Institute for Brain Science. All rights reserved.
 *  Licensed under the Allen Institute Software License; see LICENSE.
 *
 *  Read the bit-exactness contract in StptDeformableCuda.h before changing anything
 *  in the kernel. Every arithmetic expression below is transcribed from
 *  StptDeformableMRF.cxx deliberately verbatim, INCLUDING the places where it looks
 *  odd -- the float grid-step divisions, the int64/float mixing in the boundary
 *  test, the bitwise `|` between comparisons. Those are not accidents to tidy up
 *  here; they decide which voxel is read, and "cleaning them up" is how the GPU
 *  silently stops matching the CPU.
 *
 *=========================================================================*/
#include "StptDeformableCuda.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace stpt
{
namespace cuda
{
namespace
{

void
Check(cudaError_t err, const char * what)
{
  if (err != cudaSuccess)
  {
    std::ostringstream msg;
    msg << "CUDA error in " << what << ": " << cudaGetErrorString(err);
    throw std::runtime_error(msg.str());
  }
}

/** Frees a device pointer however the scope is left. */
struct DeviceBuffer
{
  void * p{ nullptr };
  DeviceBuffer() = default;
  explicit DeviceBuffer(std::size_t bytes)
  {
    if (bytes > 0) { Check(cudaMalloc(&p, bytes), "cudaMalloc"); }
  }
  ~DeviceBuffer()
  {
    if (p) { cudaFree(p); }
  }
  DeviceBuffer(const DeviceBuffer &) = delete;
  DeviceBuffer & operator=(const DeviceBuffer &) = delete;
  DeviceBuffer(DeviceBuffer && other) noexcept : p(other.p) { other.p = nullptr; }
  template <typename T>
  T * as() const
  {
    return static_cast<T *>(p);
  }
};

/** Device state, probed once. A machine with no driver must not cost a CUDA call
 *  on every query, and must not throw -- "no GPU" is a normal configuration. */
struct DeviceInfo
{
  bool        usable{ false };
  std::string description;
  std::string reason;
  std::size_t sharedBytesPerBlock{ 48 * 1024 };
};

const DeviceInfo &
Info()
{
  static const DeviceInfo info = [] {
    DeviceInfo d;
    int        count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    if (err != cudaSuccess)
    {
      d.reason = std::string("no usable CUDA device (") + cudaGetErrorString(err) + ")";
      return d;
    }
    if (count == 0)
    {
      d.reason = "no CUDA device found";
      return d;
    }
    cudaDeviceProp prop{};
    err = cudaGetDeviceProperties(&prop, 0);
    if (err != cudaSuccess)
    {
      d.reason = std::string("could not query CUDA device 0 (") + cudaGetErrorString(err) + ")";
      return d;
    }
    std::ostringstream text;
    text.setf(std::ios::fixed);
    text.precision(1);
    text << prop.name << " ("
         << static_cast<double>(prop.totalGlobalMem) / (1024.0 * 1024.0 * 1024.0)
         << " GiB, cc " << prop.major << "." << prop.minor << ")";
    d.usable = true;
    d.description = text.str();
    d.sharedBytesPerBlock = prop.sharedMemPerBlock;
    return d;
  }();
  return info;
}

__device__ inline long long
ClampIdx(long long v, long long hi)
{
  return v < 0 ? 0 : (v > hi ? hi : v);
}

/** One block per control node; threads stride over labels, each owning one label's
 *  two accumulators for the whole sample loop. That is what makes this bit-exact:
 *  the sum over samples stays sequential within a thread, and no value is ever
 *  summed across threads. */
__global__ void
OptimizationKernel(const float * __restrict__ fixed,
                   long long m,
                   long long n,
                   long long o,
                   const float * __restrict__ movingi,
                   long long mi,
                   long long ni,
                   long long oi,
                   int   gridx,
                   int   gridy,
                   float gridstep,
                   int   len,
                   int   lengthfull,
                   const int * __restrict__ yOff,
                   const int * __restrict__ xOff,
                   const int * __restrict__ zOff,
                   int maxsamp,
                   int randommode,
                   const int * __restrict__ sampleIdx, // batch-local when randommode
                   int subpixel,
                   int hw_local,
                   const unsigned char * __restrict__ activeNodes, // full array, node-indexed
                   float     epsilon,
                   long long nodeBegin,
                   float * __restrict__ energyOut) // batch-local, slot-major
{
  const long long slot = blockIdx.x;
  const long long i = nodeBegin + slot;
  float * myEnergy = energyOut + slot * static_cast<long long>(lengthfull);

  // --- node position on the control grid. Transcribed from optimization(); the
  // multiplications by `gridstep` go through float and truncate, which is where the
  // grid actually lands, so they must not become integer arithmetic here.
  long long z1 = i / (static_cast<long long>(gridx) * gridy);
  long long x1 = (i - z1 * gridx * gridy) / gridx;
  long long y1 = i - z1 * gridx * gridy - x1 * gridx;
  z1 *= gridstep;
  x1 *= gridstep;
  y1 *= gridstep;

  // A node with no evidence behind it gets a flat zero data cost.
  if (activeNodes != nullptr && activeNodes[i] == 0)
  {
    for (int l = threadIdx.x; l < lengthfull; l += blockDim.x) { myEnergy[l] = 0.0f; }
    return;
  }

  // --- boundary test, verbatim including the `|` and the x/ni, y/mi pairing.
  bool boundaries = true;
  if (subpixel)
  {
    if (x1 * 2 + (gridstep - 1) * 2 + hw_local >= ni | y1 * 2 + (gridstep - 1) * 2 + hw_local >= mi | z1 * 2 + (gridstep - 1) * 2 + hw_local >= oi)
      boundaries = false;
    if (x1 * 2 - hw_local < 0 | y1 * 2 - hw_local < 0 | z1 * 2 - hw_local < 0)
      boundaries = false;
  }
  else
  {
    if (x1 + (gridstep - 1) + hw_local >= ni | y1 + (gridstep - 1) + hw_local >= mi | z1 + (gridstep - 1) + hw_local >= oi)
      boundaries = false;
    if (x1 - hw_local < 0 | y1 - hw_local < 0 | z1 - hw_local < 0)
      boundaries = false;
  }

  // --- per-node sample table, built once and shared by every label in the block.
  // The CPU recomputes this inside its label loop; hoisting it changes no
  // arithmetic, only how often it is done. sbase comes first so its 8-byte
  // alignment holds whatever maxsamp is.
  extern __shared__ unsigned char smem[];
  long long * sbase = reinterpret_cast<long long *>(smem);
  float *     sfv = reinterpret_cast<float *>(sbase + maxsamp);
  int *       sbx = reinterpret_cast<int *>(sfv + maxsamp);
  int *       sby = sbx + maxsamp;
  int *       sbz = sby + maxsamp;

  for (int j1 = threadIdx.x; j1 < maxsamp; j1 += blockDim.x)
  {
    const int i1 = randommode ? sampleIdx[slot * static_cast<long long>(maxsamp) + j1] : j1;
    // Float divisions by the grid step, exactly as optimization() does them.
    long long zz = i1 / (gridstep * gridstep);
    long long xx = (i1 - zz * gridstep * gridstep) / gridstep;
    long long yy = i1 - zz * gridstep * gridstep - xx * gridstep;

    xx += x1;
    yy += y1;
    zz += z1;

    sfv[j1] = fixed[yy + xx * m + zz * m * n];

    const long long bx = subpixel ? xx * 2 : xx;
    const long long by = subpixel ? yy * 2 : yy;
    const long long bz = subpixel ? zz * 2 : zz;
    sbx[j1] = static_cast<int>(bx);
    sby[j1] = static_cast<int>(by);
    sbz[j1] = static_cast<int>(bz);
    sbase[j1] = by + bx * mi + bz * mi * ni;
  }
  __syncthreads();

  // --- normA. Every thread recomputes the same sequential sum rather than reducing
  // across threads: a tree reduction would reassociate it and change the last bits.
  // All threads read the same shared address each iteration, so it broadcasts.
  float normA = 0.0f;
  for (int j1 = 0; j1 < maxsamp; j1++)
  {
    const float fv = sfv[j1];
    normA += fv * fv;
  }

  // --- one label per thread, sample loop sequential inside it.
  for (int l = threadIdx.x; l < lengthfull; l += blockDim.x)
  {
    // Label l is i2 + j2*len + k2*len*len, and the three offset tables are indexed
    // one per component -- ys by i2, xs by j2, zs by k2.
    const int i2 = l % len;
    const int j2 = (l / len) % len;
    const int k2 = l / (len * len);

    float c = 0.0f;
    float nb = 0.0f;

    if (boundaries)
    {
      // Interior: the label offset is constant across samples, so it is one add
      // onto the precomputed per-sample base. Consecutive threads hold consecutive
      // i2, hence consecutive addresses, so the warp's reads coalesce.
      const long long lo = static_cast<long long>(zOff[k2]) * mi * ni +
                           static_cast<long long>(xOff[j2]) * mi + yOff[i2];
      for (int j1 = 0; j1 < maxsamp; j1++)
      {
        const float mv = movingi[sbase[j1] + lo];
        c += sfv[j1] * mv;
        nb += mv * mv;
      }
    }
    else
    {
      // Edge node: each axis clamps independently, and the clamp depends on both
      // the sample and the label, so it cannot be hoisted out of either loop.
      const int zo = zOff[k2];
      const int xo = xOff[j2];
      const int yo = yOff[i2];
      for (int j1 = 0; j1 < maxsamp; j1++)
      {
        const long long idx = ClampIdx(static_cast<long long>(sbz[j1]) + zo, oi - 1) * mi * ni +
                              ClampIdx(static_cast<long long>(sbx[j1]) + xo, ni - 1) * mi +
                              ClampIdx(static_cast<long long>(sby[j1]) + yo, mi - 1);
        const float mv = movingi[idx];
        c += sfv[j1] * mv;
        nb += mv * mv;
      }
    }

    // Float ratio then a double square, exactly as the CPU readout does it.
    const float q = c * c / normA / nb;
    myEnergy[l] = static_cast<float>(epsilon * (1.0 - static_cast<double>(q) * static_cast<double>(q)));
  }
}

} // namespace

bool
BuiltWithCuda()
{
  return true;
}

bool
DeviceAvailable()
{
  return Info().usable;
}

std::string
DeviceDescription()
{
  return Info().description;
}

std::string
UnavailableReason()
{
  return Info().usable ? std::string() : Info().reason;
}

bool
RunOptimization(const OptimizationInputs & in)
{
  const DeviceInfo & dev = Info();
  if (!dev.usable) { return false; }
  if (in.sz1 <= 0 || in.lengthfull <= 0 || in.maxsamp <= 0) { return false; }

  // The per-node sample table lives in shared memory. Beyond what fits, this kernel
  // has nothing to offer and the CPU path is the right answer -- only reachable with
  // a --samples far above the default 1024.
  const std::size_t sharedBytes = static_cast<std::size_t>(in.maxsamp) * 24;
  if (sharedBytes > dev.sharedBytesPerBlock) { return false; }

  const std::size_t fixedBytes = static_cast<std::size_t>(in.m * in.n * in.o) * sizeof(float);
  const std::size_t movingBytes = static_cast<std::size_t>(in.mi * in.ni * in.oi) * sizeof(float);

  // Both volumes stay resident for the whole call; only the cost volume is batched.
  // If they do not fit there is no point starting -- fall back rather than thrash.
  std::size_t freeBytes = 0;
  std::size_t totalBytes = 0;
  if (cudaMemGetInfo(&freeBytes, &totalBytes) != cudaSuccess) { return false; }
  const std::size_t perNodeEnergy = static_cast<std::size_t>(in.lengthfull) * sizeof(float);
  if (fixedBytes + movingBytes + perNodeEnergy + (64u << 20) > freeBytes) { return false; }

  // Batch the cost volume so it never has to be resident whole: at CCF 10 um the
  // finest level's volume is ~9 GB, which no consumer card holds.
  const std::size_t energyBudget = std::min<std::size_t>(
    256u << 20, freeBytes - fixedBytes - movingBytes - (64u << 20));
  long long batchNodes = static_cast<long long>(energyBudget / perNodeEnergy);
  batchNodes = std::max<long long>(1, std::min<long long>(batchNodes, in.sz1));

  DeviceInfo const & d = dev;
  (void)d;

  DeviceBuffer dFixed(fixedBytes);
  DeviceBuffer dMoving(movingBytes);
  DeviceBuffer dYOff(static_cast<std::size_t>(in.len) * sizeof(int));
  DeviceBuffer dXOff(static_cast<std::size_t>(in.len) * sizeof(int));
  DeviceBuffer dZOff(static_cast<std::size_t>(in.len) * sizeof(int));
  DeviceBuffer dActive(in.activeNodes ? static_cast<std::size_t>(in.sz1) : 0);
  DeviceBuffer dEnergy(static_cast<std::size_t>(batchNodes) * perNodeEnergy);
  DeviceBuffer dSamples(in.randommode ? static_cast<std::size_t>(batchNodes) *
                                          in.maxsamp * sizeof(int)
                                      : 0);

  Check(cudaMemcpy(dFixed.p, in.fixed, fixedBytes, cudaMemcpyHostToDevice), "copy fixed");
  Check(cudaMemcpy(dMoving.p, in.movingi, movingBytes, cudaMemcpyHostToDevice), "copy moving");
  const std::size_t offBytes = static_cast<std::size_t>(in.len) * sizeof(int);
  Check(cudaMemcpy(dYOff.p, in.yOff, offBytes, cudaMemcpyHostToDevice), "copy yOff");
  Check(cudaMemcpy(dXOff.p, in.xOff, offBytes, cudaMemcpyHostToDevice), "copy xOff");
  Check(cudaMemcpy(dZOff.p, in.zOff, offBytes, cudaMemcpyHostToDevice), "copy zOff");
  if (in.activeNodes)
  {
    Check(cudaMemcpy(dActive.p, in.activeNodes, static_cast<std::size_t>(in.sz1),
                     cudaMemcpyHostToDevice),
          "copy activeNodes");
  }

  // Enough threads to cover a warp's worth of labels, capped at 256: beyond that the
  // per-thread sample loop, not occupancy, is the limit.
  int threads = 256;
  while (threads > 32 && threads >= in.lengthfull * 2) { threads /= 2; }

  for (long long begin = 0; begin < in.sz1; begin += batchNodes)
  {
    const long long count = std::min<long long>(batchNodes, in.sz1 - begin);

    if (in.randommode)
    {
      Check(cudaMemcpy(dSamples.p, in.sampleIdx + begin * in.maxsamp,
                       static_cast<std::size_t>(count) * in.maxsamp * sizeof(int),
                       cudaMemcpyHostToDevice),
            "copy sampleIdx");
    }

    OptimizationKernel<<<static_cast<unsigned int>(count), threads, sharedBytes>>>(
      dFixed.as<const float>(), in.m, in.n, in.o,
      dMoving.as<const float>(), in.mi, in.ni, in.oi,
      in.gridx, in.gridy, in.gridstep,
      in.len, in.lengthfull,
      dYOff.as<const int>(), dXOff.as<const int>(), dZOff.as<const int>(),
      in.maxsamp, in.randommode ? 1 : 0, dSamples.as<const int>(),
      in.subpixel ? 1 : 0, in.hw_local,
      in.activeNodes ? dActive.as<const unsigned char>() : nullptr,
      in.epsilon, begin, dEnergy.as<float>());

    Check(cudaGetLastError(), "kernel launch");
    Check(cudaMemcpy(in.energy + begin * in.lengthfull, dEnergy.p,
                     static_cast<std::size_t>(count) * perNodeEnergy,
                     cudaMemcpyDeviceToHost),
          "copy energy back");
  }

  return true;
}

} // namespace cuda
} // namespace stpt
