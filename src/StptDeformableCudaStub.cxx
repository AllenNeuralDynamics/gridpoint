/*=========================================================================
 *
 *  StptDeformableCudaStub.cxx
 *
 *  The CUDA backend as it exists in a build without CUDA. Compiled INSTEAD OF
 *  StptDeformableCuda.cu when GRIDPOINT_ENABLE_CUDA is off.
 *
 *  Having a stub rather than #ifdefs at the call site is deliberate: the deformable
 *  stage and the CLI call exactly the same functions either way, so there is one
 *  code path to read and no configuration in which an untested branch compiles for
 *  the first time.
 *
 *  Copyright (c) Allen Institute for Brain Science. All rights reserved.
 *  Licensed under the Allen Institute Software License; see LICENSE.
 *
 *=========================================================================*/
#include "StptDeformableCuda.h"

namespace stpt
{
namespace cuda
{

bool
BuiltWithCuda()
{
  return false;
}

bool
DeviceAvailable()
{
  return false;
}

std::string
DeviceDescription()
{
  return std::string();
}

std::string
UnavailableReason()
{
  return "this build has no CUDA support (configure with -DGRIDPOINT_ENABLE_CUDA=ON)";
}

bool
RunOptimization(const OptimizationInputs &)
{
  return false;
}

} // namespace cuda
} // namespace stpt
