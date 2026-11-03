/*=========================================================================
 *
 *  StptLinearRegistration.h
 *
 *  Linear (rigid + affine) mutual-information registration, derived from
 *  AllenInstitute/stpt_registration:
 *
 *    lims2_modules/common/register/Registration/
 *      idpMutualInfoVersorRigidVolumeRegistration.{h,txx}
 *      idpMutualInfoCenteredAffineRegularGradientDescentRegistration.{h,txx}
 *
 *  Copyright (c) Allen Institute for Brain Science. All rights reserved.
 *  Licensed under the Allen Institute Software License; see LICENSE.
 *
 *  Reworked to depend on ITK alone (no Kakadu/JPEG2000, TinyXML, Boost or VTK)
 *  and to register two arbitrary 3D volumes rather than a LIMS image series.
 *
 *=========================================================================*/
#ifndef StptLinearRegistration_h
#define StptLinearRegistration_h

#include "itkAffineTransform.h"
#include "itkImage.h"
#include "itkImageMaskSpatialObject.h"
#include "itkVersorRigid3DTransform.h"

namespace stpt
{

constexpr unsigned int Dimension = 3;

using PixelType = float;
using ImageType = itk::Image<PixelType, Dimension>;
using MaskPixelType = unsigned char;
using MaskImageType = itk::Image<MaskPixelType, Dimension>;
using MaskSpatialObjectType = itk::ImageMaskSpatialObject<Dimension>;

using RigidTransformType = itk::VersorRigid3DTransform<double>;
using AffineTransformType = itk::AffineTransform<double, Dimension>;

/** Tuning knobs for the linear stages. Defaults reproduce the upstream Allen values
 *  (64 histogram bins, 250000 spatial samples, 100 iterations, step 60 -> 0.1,
 *  relaxation 0.8, translation scale 1e7, fixed RNG seed). */
struct LinearParameters
{
  unsigned int numberOfLevels{ 4 };
  unsigned int numberOfIterations{ 100 };
  unsigned int numberOfHistogramBins{ 64 };
  unsigned int numberOfSpatialSamples{ 250000 };
  double       maximumStepLength{ 60.0 };
  double       minimumStepLength{ 0.1 };
  double       maximumStepLengthRate{ 0.3 };
  double       minimumStepLengthRate{ 0.15 };
  double       relaxationFactor{ 0.8 };
  double       relaxationFactorRate{ 1.1 };
  double       angleScale{ 1.0 };
  double       matrixScale{ 1.0 };
  double       translationScale{ 1.0e7 };
  int          randomSeed{ 76926294 };
  bool         verbose{ false };
};

/** Optional metric masks; either member may be null. */
struct MaskPair
{
  MaskSpatialObjectType::Pointer fixedMask;
  MaskSpatialObjectType::Pointer movingMask;
};

/** Rigid (versor) mutual-information registration.
 *
 *  Initialized from image moments (centre of gravity), mirroring the upstream
 *  CenteredVersorTransformInitializer with MomentsOn(). */
RigidTransformType::Pointer
RunRigidRegistration(const ImageType *       fixed,
                     const ImageType *       moving,
                     const MaskPair &        masks,
                     const LinearParameters & params);

/** Affine mutual-information registration.
 *
 *  When `initial` is non-null the affine starts from that rigid result, so the
 *  returned transform is the collapsed rigid+affine mapping. Otherwise it is
 *  initialized from image moments. */
AffineTransformType::Pointer
RunAffineRegistration(const ImageType *          fixed,
                      const ImageType *          moving,
                      const MaskPair &           masks,
                      const LinearParameters &   params,
                      const RigidTransformType * initial);

} // namespace stpt

#endif
