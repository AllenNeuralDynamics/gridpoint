/*=========================================================================
 *
 *  StptDeformableMRF.h
 *
 *  Discrete symmetric deformable registration -- the distinctive part of the Allen
 *  STPT/CCF method. Displacements live on a B-spline-style control grid and are chosen
 *  by discrete optimization (minimum-spanning-tree message passing over a label space)
 *  rather than continuous gradient descent, and the registration is run forward and
 *  backward each level and composed, giving a symmetric, invertible field.
 *
 *  Derived from AllenInstitute/stpt_registration:
 *    lims2_modules/common/register/Registration/idpProjectionLocalAlignmentDispatcher.txx
 *    lims2_modules/common/register/Common/idpRegistrationUtilities.{h,cxx}  (Relation)
 *
 *  Copyright (c) Allen Institute for Brain Science. All rights reserved.
 *  Licensed under the Allen Institute Software License; see LICENSE.
 *
 *  Upstream these are methods on a class that inherits the LIMS-bound registration
 *  dispatcher. Here they are free functions over plain buffers, so this translation unit
 *  pulls in no ITK and no pipeline machinery at all.
 *
 *=========================================================================*/
#ifndef StptDeformableMRF_h
#define StptDeformableMRF_h

#include <vector>

namespace stpt
{

/** Tuning knobs for the discrete deformable stage. Defaults are the upstream values. */
struct DeformableParameters
{
  int          numLevels{ 4 };      // coarse-to-fine levels
  int          gridStart{ 30 };     // coarsest control-grid spacing, voxels; halved per level
  float        lambda{ 1.0e-5f };   // smoothness weight on the MRF pairwise term
  int          sampleNum{ 1024 };   // random samples per node, quartered each level (floor 64)
  bool         verbose{ false };
};

/** Result of the deformable stage. All buffers are x-fastest, size sizeX*sizeY*sizeZ. */
struct DeformableResult
{
  std::vector<float> warped;  // moving image resampled into fixed space
  // Forward displacement field in VOXELS, already mapped to x-fastest axis order.
  // Maps fixed positions to moving positions: use it to pull moving into fixed space.
  std::vector<float> dispX;
  std::vector<float> dispY;
  std::vector<float> dispZ;
  // Backward displacement field, same units and axis order. The symmetric optimization
  // produces this alongside the forward one -- it is the other half of the composed
  // path, not a numerical inversion -- and it maps moving positions back to fixed ones.
  // Both fields live on the fixed grid, because the deformable stage runs on the pair
  // after linear alignment has already resampled moving onto that grid.
  std::vector<float> invX;
  std::vector<float> invY;
  std::vector<float> invZ;
  float              metricBefore{ 0.0f };  // mean squared difference before
  float              metricAfter{ 0.0f };   // ... and after
  int                levelsUsed{ 0 };
  int                gridStartUsed{ 0 };
};

/** Rough peak resident bytes for RunDeformableRegistration on this shape.
 *
 *  Calibrated against measured peak working set at 64/96/128/160/192 cubed: about
 *  85 bytes per voxel of persistent and scratch buffers, plus the largest single
 *  level's discrete cost volumes (which depend on the grid schedule, so they are
 *  computed rather than folded into the per-voxel figure), plus a small baseline.
 *  Deliberately biased high: over 96..192 cubed it lands between 0.98x and 1.12x of
 *  measurement, and a pre-flight that under-predicts is worse than useless. The 64
 *  cubed end reads 1.24x, where the fixed baseline dominates. A planning number, not
 *  a guarantee -- it counts allocations, not allocator slack or page behaviour.
 */
long long
EstimateDeformablePeakBytes(int sizeX, int sizeY, int sizeZ, const DeformableParameters & params);

/** Run the discrete symmetric deformable registration.
 *
 *  `imgFix` and `imgMov` are x-fastest buffers of size sizeX*sizeY*sizeZ, already
 *  linearly aligned (run the rigid/affine stages first). The returned field maps fixed
 *  positions to moving positions, in voxels.
 *
 *  `fixedMask` and `movingMask` are optional buffers of the same size and layout as
 *  `imgFix`, non-zero meaning inside, both already on the fixed grid -- for the moving
 *  mask that means the caller has carried it through the linear transform, since this
 *  stage sees the moving image only after that transform resampled it.
 *
 *  They are kept separate rather than intersected by the caller because this stage
 *  re-derives the masked control grid at every pyramid level: the moving mask is
 *  displaced by the current field before being tested, and in the backward half of the
 *  symmetric optimization the two swap roles, the fixed mask being displaced by the
 *  inverse field. Intersecting them up front would freeze both at the linear
 *  alignment. A control-grid node whose cell contains no in-mask
 *  voxel is given a flat zero data cost, so the regularizer alone decides its
 *  displacement and it follows its neighbours smoothly rather than chasing whatever lies
 *  outside the mask. Such nodes stay in the graph, so the spanning tree cannot fragment,
 *  and the field stays defined everywhere.
 *
 *  Note what this does NOT mean: masked-out regions are not held still. Their
 *  displacement is a smooth extrapolation from the masked-in neighbourhood, which is
 *  usually what is wanted, but it is not zero. */
DeformableResult
RunDeformableRegistration(const float *                imgFix,
                          const float *                imgMov,
                          int                          sizeX,
                          int                          sizeY,
                          int                          sizeZ,
                          const DeformableParameters & params,
                          const unsigned char *        fixedMask = nullptr,
                          const unsigned char *        movingMask = nullptr);

} // namespace stpt

#endif
