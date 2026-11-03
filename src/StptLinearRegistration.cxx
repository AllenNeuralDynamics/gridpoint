/*=========================================================================
 *
 *  StptLinearRegistration.cxx
 *
 *  Derived from AllenInstitute/stpt_registration (see the header for the exact
 *  upstream files). Copyright (c) Allen Institute for Brain Science.
 *  Licensed under the Allen Institute Software License; see LICENSE.
 *
 *=========================================================================*/
#include "StptLinearRegistration.h"

#include "itkCenteredTransformInitializer.h"
#include "itkCenteredVersorTransformInitializer.h"
#include "itkCommand.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkMultiResolutionPyramidImageFilter.h"
#include "itkRegularStepGradientDescentBaseOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkVersorRigid3DTransformOptimizer.h"

#include <iostream>

namespace stpt
{
namespace
{

using MetricType = itk::MattesMutualInformationImageToImageMetric<ImageType, ImageType>;
using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, double>;
using PyramidType = itk::MultiResolutionPyramidImageFilter<ImageType, ImageType>;

/** Reduce the optimizer step sizes once past the coarsest level.
 *
 *  Mirrors the upstream BeforeIteration(): at every level > 0 the step lengths are
 *  scaled by their *Rate factors and the relaxation factor is loosened, all computed
 *  from the original settings rather than cumulatively. */
template <typename TRegistration>
class LevelObserver : public itk::Command
{
public:
  using Self = LevelObserver;
  using Superclass = itk::Command;
  using Pointer = itk::SmartPointer<Self>;
  itkNewMacro(Self);

  void
  SetOptimizer(itk::RegularStepGradientDescentBaseOptimizer * optimizer)
  {
    m_Optimizer = optimizer;
  }
  void
  SetParameters(const LinearParameters & params)
  {
    m_Params = params;
  }

  void
  Execute(itk::Object * caller, const itk::EventObject & event) override
  {
    if (!itk::IterationEvent().CheckEvent(&event) || m_Optimizer == nullptr)
    {
      return;
    }
    const auto * registration = static_cast<const TRegistration *>(caller);
    const unsigned int level = registration->GetCurrentLevel();
    if (level == 0)
    {
      return;
    }
    m_Optimizer->SetMaximumStepLength(m_Params.maximumStepLength * m_Params.maximumStepLengthRate);
    m_Optimizer->SetMinimumStepLength(m_Params.minimumStepLength * m_Params.minimumStepLengthRate);
    m_Optimizer->SetRelaxationFactor(m_Params.relaxationFactor * m_Params.relaxationFactorRate);
    if (m_Params.verbose)
    {
      std::cout << "  level " << level << ": step "
                << m_Optimizer->GetMaximumStepLength() << " -> "
                << m_Optimizer->GetMinimumStepLength() << std::endl;
    }
  }

  void
  Execute(const itk::Object *, const itk::EventObject &) override
  {}

protected:
  LevelObserver() = default;

private:
  itk::RegularStepGradientDescentBaseOptimizer * m_Optimizer{ nullptr };
  LinearParameters                               m_Params;
};

/** Metric configured the upstream way, with optional masks.
 *
 *  The upstream sample count (250000) assumes a full-size brain volume. Asking for more
 *  samples than the fixed image has voxels makes the Mattes metric fail outright, so fall
 *  back to using every pixel in that case -- which is what the request really meant. */
MetricType::Pointer
MakeMetric(const ImageType * fixed, const MaskPair & masks, const LinearParameters & params)
{
  auto metric = MetricType::New();
  metric->SetNumberOfHistogramBins(params.numberOfHistogramBins);

  const auto availablePixels = fixed->GetBufferedRegion().GetNumberOfPixels();
  if (static_cast<itk::SizeValueType>(params.numberOfSpatialSamples) >= availablePixels)
  {
    metric->SetUseAllPixels(true);
    if (params.verbose)
    {
      std::cout << "  metric: using all " << availablePixels << " pixels (requested "
                << params.numberOfSpatialSamples << " samples)" << std::endl;
    }
  }
  else
  {
    metric->SetNumberOfSpatialSamples(params.numberOfSpatialSamples);
  }
  metric->ReinitializeSeed(params.randomSeed);
  if (masks.fixedMask)
  {
    metric->SetFixedImageMask(masks.fixedMask);
  }
  if (masks.movingMask)
  {
    metric->SetMovingImageMask(masks.movingMask);
  }
  return metric;
}

} // namespace

RigidTransformType::Pointer
RunRigidRegistration(const ImageType *        fixed,
                     const ImageType *        moving,
                     const MaskPair &         masks,
                     const LinearParameters & params)
{
  using RegistrationType = itk::MultiResolutionImageRegistrationMethod<ImageType, ImageType>;
  using OptimizerType = itk::VersorRigid3DTransformOptimizer;

  auto transform = RigidTransformType::New();

  // Centre-of-gravity initialization (upstream MomentsOn()).
  using InitializerType = itk::CenteredVersorTransformInitializer<ImageType, ImageType>;
  auto initializer = InitializerType::New();
  initializer->SetTransform(transform);
  initializer->SetFixedImage(fixed);
  initializer->SetMovingImage(moving);
  initializer->MomentsOn();
  initializer->InitializeTransform();
  if (params.verbose)
  {
    std::cout << "  rigid init: " << transform->GetParameters() << std::endl;
  }

  auto optimizer = OptimizerType::New();
  OptimizerType::ScalesType scales(transform->GetNumberOfParameters());
  scales[0] = 1.0 / params.angleScale;
  scales[1] = 1.0 / params.angleScale;
  scales[2] = 1.0 / params.angleScale;
  scales[3] = 1.0 / params.translationScale;
  scales[4] = 1.0 / params.translationScale;
  scales[5] = 1.0 / params.translationScale;
  optimizer->SetScales(scales);
  optimizer->SetNumberOfIterations(params.numberOfIterations);
  optimizer->SetMaximumStepLength(params.maximumStepLength);
  optimizer->SetMinimumStepLength(params.minimumStepLength);
  optimizer->SetRelaxationFactor(params.relaxationFactor);
  optimizer->MinimizeOn();

  auto registration = RegistrationType::New();
  registration->SetTransform(transform);
  registration->SetMetric(MakeMetric(fixed, masks, params));
  registration->SetOptimizer(optimizer);
  registration->SetInterpolator(InterpolatorType::New());
  registration->SetFixedImagePyramid(PyramidType::New());
  registration->SetMovingImagePyramid(PyramidType::New());
  registration->SetFixedImage(fixed);
  registration->SetMovingImage(moving);
  registration->SetFixedImageRegion(fixed->GetBufferedRegion());
  registration->SetNumberOfLevels(params.numberOfLevels);
  registration->SetInitialTransformParameters(transform->GetParameters());

  auto observer = LevelObserver<RegistrationType>::New();
  observer->SetOptimizer(optimizer);
  observer->SetParameters(params);
  registration->AddObserver(itk::IterationEvent(), observer);

  registration->Update();

  transform->SetParameters(registration->GetLastTransformParameters());
  if (params.verbose)
  {
    std::cout << "rigid: " << optimizer->GetStopConditionDescription() << std::endl;
  }
  return transform;
}

AffineTransformType::Pointer
RunAffineRegistration(const ImageType *          fixed,
                      const ImageType *          moving,
                      const MaskPair &           masks,
                      const LinearParameters &   params,
                      const RigidTransformType * initial)
{
  using RegistrationType = itk::MultiResolutionImageRegistrationMethod<ImageType, ImageType>;
  using OptimizerType = itk::RegularStepGradientDescentOptimizer;

  auto transform = AffineTransformType::New();
  if (initial != nullptr)
  {
    // Warm-start from the rigid result, so the affine output is the collapsed
    // rigid+affine transform.
    transform->SetCenter(initial->GetCenter());
    transform->SetMatrix(initial->GetMatrix());
    transform->SetTranslation(initial->GetTranslation());
  }
  else
  {
    using InitializerType = itk::CenteredTransformInitializer<AffineTransformType, ImageType, ImageType>;
    auto initializer = InitializerType::New();
    initializer->SetTransform(transform);
    initializer->SetFixedImage(fixed);
    initializer->SetMovingImage(moving);
    initializer->MomentsOn();
    initializer->InitializeTransform();
  }
  if (params.verbose)
  {
    std::cout << "  affine init: " << transform->GetParameters() << std::endl;
  }

  auto optimizer = OptimizerType::New();
  OptimizerType::ScalesType scales(transform->GetNumberOfParameters());
  // 9 matrix entries then 3 translation entries (upstream indexes 0..11).
  for (unsigned int i = 0; i < 9; ++i)
  {
    scales[i] = 1.0 / params.matrixScale;
  }
  for (unsigned int i = 9; i < transform->GetNumberOfParameters(); ++i)
  {
    scales[i] = 1.0 / params.translationScale;
  }
  optimizer->SetScales(scales);
  optimizer->SetNumberOfIterations(params.numberOfIterations);
  optimizer->SetMaximumStepLength(params.maximumStepLength);
  optimizer->SetMinimumStepLength(params.minimumStepLength);
  optimizer->SetRelaxationFactor(params.relaxationFactor);
  optimizer->MinimizeOn();

  auto registration = RegistrationType::New();
  registration->SetTransform(transform);
  registration->SetMetric(MakeMetric(fixed, masks, params));
  registration->SetOptimizer(optimizer);
  registration->SetInterpolator(InterpolatorType::New());
  registration->SetFixedImagePyramid(PyramidType::New());
  registration->SetMovingImagePyramid(PyramidType::New());
  registration->SetFixedImage(fixed);
  registration->SetMovingImage(moving);
  registration->SetFixedImageRegion(fixed->GetBufferedRegion());
  registration->SetNumberOfLevels(params.numberOfLevels);
  registration->SetInitialTransformParameters(transform->GetParameters());

  auto observer = LevelObserver<RegistrationType>::New();
  observer->SetOptimizer(optimizer);
  observer->SetParameters(params);
  registration->AddObserver(itk::IterationEvent(), observer);

  registration->Update();

  transform->SetParameters(registration->GetLastTransformParameters());
  if (params.verbose)
  {
    std::cout << "affine: " << optimizer->GetStopConditionDescription() << std::endl;
  }
  return transform;
}

} // namespace stpt
