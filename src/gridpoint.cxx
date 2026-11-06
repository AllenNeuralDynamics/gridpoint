/*=========================================================================
 *
 *  gridpoint.cxx -- register two 3D volumes with the Allen STPT/CCF algorithm.
 *
 *  Derived from AllenInstitute/stpt_registration.
 *  Copyright (c) Allen Institute for Brain Science. All rights reserved.
 *  Licensed under the Allen Institute Software License; see LICENSE.
 *
 *  Depends on ITK alone. Unlike the upstream idpProjection* modules, this takes two
 *  ordinary 3D volumes rather than a LIMS image series + model directory.
 *
 *=========================================================================*/
#include "StptDeformableMRF.h"
#include "StptLinearRegistration.h"

#include "itkIdentityTransform.h"
#include "itkImage.h"
#include "itkMultiThreaderBase.h"
#include "itkVector.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkTransformFileWriter.h"

#include <itksys/SystemInformation.hxx>
#include <itksys/SystemTools.hxx>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{

void
PrintUsage(const char * argv0)
{
  std::cerr
    << "Usage: " << argv0 << " -f <fixed> -m <moving> -o <out_dir> [options]\n\n"
    << "Registers two 3D volumes with the Allen STPT/CCF algorithm.\n\n"
    << "Required:\n"
    << "  -f <path>                    fixed (reference) image\n"
    << "  -m <path>                    moving image\n"
    << "  -o <dir>                     output directory\n\n"
    << "Options:\n"
    << "  --stages a,b,c               rigid,affine,deformable (default: rigid,affine,deformable)\n"
    << "  --fixed-mask <path>          metric mask for the fixed image\n"
    << "  --moving-mask <path>         metric mask for the moving image\n"
    << "  --threads <n>                maximum number of threads\n"
    << "  --bins <n>                   MI histogram bins (default 64)\n"
    << "  --linear-samples <n>         MI spatial samples (default 250000)\n"
    << "  --linear-iterations <n>      optimizer iterations per level (default 100)\n"
    << "  --levels <n>                 multi-resolution levels (default 4 => shrink 8,4,2,1)\n"
    << "  --translation-scale <x>      optimizer translation scale (default 1e7, upstream).\n"
    << "                               NOTE: tuned for micron-spaced CCF volumes; on images\n"
    << "                               with ~unit spacing use a value near the physical\n"
    << "                               extent (e.g. 1000) or the optimizer will diverge.\n"
    << "  --angle-scale <x>            optimizer rotation scale (default 1)\n"
    << "  --matrix-scale <x>           optimizer affine-matrix scale (default 1)\n"
    << "  --max-step <x>               optimizer maximum step length (default 60)\n"
    << "  --min-step <x>               optimizer minimum step length (default 0.1)\n"
    << "  --grid a,b,c,d               deformable control-grid schedule, voxels\n"
    << "  --lambda <x>                 deformable smoothness weight\n"
    << "  --samples <n>                deformable random samples\n"
    << "  --warped-ext <.ext>          warped image extension (default .nii.gz)\n"
    << "  --write-deformation-field    also write the dense deformation fields\n"
    << "  -v                           verbose\n\n"
    << "Outputs, written into <out_dir>:\n"
    << "  stpt_Affine.txt              collapsed rigid+affine transform (ITK, maps fixed->moving)\n"
    << "  stpt_Warped<ext>             moving image resampled onto the fixed grid\n"
    << "  stpt_Warp.nii.gz             forward deformation field, fixed->moving\n"
    << "  stpt_InverseWarp.nii.gz      backward deformation field, moving->fixed\n"
    << "                               (both need --write-deformation-field)\n"
    << "  stpt.log                     settings and outcome\n";
}

std::string
HumanBytes(long long bytes)
{
  const char * units[] = { "B", "KiB", "MiB", "GiB", "TiB" };
  double       value = static_cast<double>(bytes);
  int          unit = 0;
  while (value >= 1024.0 && unit < 4)
  {
    value /= 1024.0;
    ++unit;
  }
  std::ostringstream out;
  out.setf(std::ios::fixed);
  out.precision(value < 10.0 ? 2 : 1);
  out << value << " " << units[unit];
  return out.str();
}

std::vector<std::string>
Split(const std::string & text, char delimiter)
{
  std::vector<std::string> parts;
  std::stringstream        stream(text);
  std::string              item;
  while (std::getline(stream, item, delimiter))
  {
    if (!item.empty())
    {
      parts.push_back(item);
    }
  }
  return parts;
}

/** Read an optional mask and wrap it as a metric mask. Returns null when path is empty. */
stpt::MaskSpatialObjectType::Pointer
ReadMask(const std::string & path)
{
  if (path.empty())
  {
    return nullptr;
  }
  using ReaderType = itk::ImageFileReader<stpt::MaskImageType>;
  auto reader = ReaderType::New();
  reader->SetFileName(path);
  reader->Update();

  auto mask = stpt::MaskSpatialObjectType::New();
  mask->SetImage(reader->GetOutput());
  mask->Update();
  return mask;
}

} // namespace

int
main(int argc, char * argv[])
{
  std::string fixedPath, movingPath, outDir;
  std::string fixedMaskPath, movingMaskPath;
  std::string stagesText = "rigid,affine,deformable";
  std::string warpedExt = ".nii.gz";
  std::string gridText;
  bool        writeDeformationField = false;
  int         threads = 0;

  stpt::LinearParameters     params;
  stpt::DeformableParameters defParams;

  // The upstream optimizer constants are calibrated to the Allen data's physical units and
  // diverge on images of other sizes, so by default we derive them from the image geometry.
  // --upstream-optimizer restores the original values; explicit flags always win.
  bool upstreamOptimizer = false;
  bool userSetTranslationScale = false;
  bool userSetMaxStep = false;
  bool userSetMinStep = false;

  for (int i = 1; i < argc; ++i)
  {
    const std::string arg = argv[i];
    auto              next = [&](const char * flag) -> std::string {
      if (i + 1 >= argc)
      {
        throw std::runtime_error(std::string("missing value for ") + flag);
      }
      return argv[++i];
    };

    try
    {
      if (arg == "-f") { fixedPath = next("-f"); }
      else if (arg == "-m") { movingPath = next("-m"); }
      else if (arg == "-o") { outDir = next("-o"); }
      else if (arg == "--stages") { stagesText = next("--stages"); }
      else if (arg == "--fixed-mask") { fixedMaskPath = next("--fixed-mask"); }
      else if (arg == "--moving-mask") { movingMaskPath = next("--moving-mask"); }
      else if (arg == "--threads") { threads = std::stoi(next("--threads")); }
      else if (arg == "--bins") { params.numberOfHistogramBins = std::stoul(next("--bins")); }
      else if (arg == "--linear-samples") { params.numberOfSpatialSamples = std::stoul(next("--linear-samples")); }
      else if (arg == "--linear-iterations") { params.numberOfIterations = std::stoul(next("--linear-iterations")); }
      else if (arg == "--levels") { params.numberOfLevels = std::stoul(next("--levels")); defParams.numLevels = static_cast<int>(params.numberOfLevels); }
      else if (arg == "--angle-scale") { params.angleScale = std::stod(next("--angle-scale")); }
      else if (arg == "--matrix-scale") { params.matrixScale = std::stod(next("--matrix-scale")); }
      else if (arg == "--translation-scale") { params.translationScale = std::stod(next("--translation-scale")); userSetTranslationScale = true; }
      else if (arg == "--max-step") { params.maximumStepLength = std::stod(next("--max-step")); userSetMaxStep = true; }
      else if (arg == "--min-step") { params.minimumStepLength = std::stod(next("--min-step")); userSetMinStep = true; }
      else if (arg == "--upstream-optimizer") { upstreamOptimizer = true; }
      else if (arg == "--grid") { gridText = next("--grid"); }
      else if (arg == "--lambda") { defParams.lambda = std::stof(next("--lambda")); }
      else if (arg == "--samples") { defParams.sampleNum = std::stoi(next("--samples")); }
      else if (arg == "--warped-ext") { warpedExt = next("--warped-ext"); }
      else if (arg == "--write-deformation-field") { writeDeformationField = true; }
      else if (arg == "-v") { params.verbose = true; }
      else if (arg == "-h" || arg == "--help") { PrintUsage(argv[0]); return EXIT_SUCCESS; }
      else
      {
        std::cerr << "error: unknown option " << arg << "\n";
        PrintUsage(argv[0]);
        return EXIT_FAILURE;
      }
    }
    catch (const std::exception & e)
    {
      std::cerr << "error: " << e.what() << "\n";
      return EXIT_FAILURE;
    }
  }

  if (fixedPath.empty() || movingPath.empty() || outDir.empty())
  {
    std::cerr << "error: -f, -m and -o are required\n";
    PrintUsage(argv[0]);
    return EXIT_FAILURE;
  }

  const std::vector<std::string> stages = Split(stagesText, ',');
  bool                           wantRigid = false, wantAffine = false, wantDeformable = false;
  for (const auto & stage : stages)
  {
    if (stage == "rigid") { wantRigid = true; }
    else if (stage == "affine") { wantAffine = true; }
    else if (stage == "deformable") { wantDeformable = true; }
    else
    {
      std::cerr << "error: unknown stage '" << stage << "' (use rigid, affine, deformable)\n";
      return EXIT_FAILURE;
    }
  }
  if (!wantRigid && !wantAffine && !wantDeformable)
  {
    std::cerr << "error: at least one stage is required\n";
    return EXIT_FAILURE;
  }
  if (!gridText.empty())
  {
    // A full schedule is accepted, but the algorithm derives it by halving: take the
    // coarsest spacing and the number of entries as the level count.
    const std::vector<std::string> grid = Split(gridText, ',');
    if (grid.empty())
    {
      std::cerr << "error: --grid needs at least one spacing" << std::endl;
      return EXIT_FAILURE;
    }
    defParams.gridStart = std::stoi(grid[0]);
    defParams.numLevels = static_cast<int>(grid.size());
  }
  if (writeDeformationField && !wantDeformable)
  {
    std::cerr << "error: --write-deformation-field needs the deformable stage" << std::endl;
    return EXIT_FAILURE;
  }

  try
  {
    if (threads > 0)
    {
      itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(threads);
      itk::MultiThreaderBase::SetGlobalDefaultNumberOfThreads(threads);
    }

    itksys::SystemTools::MakeDirectory(outDir);

    using ReaderType = itk::ImageFileReader<stpt::ImageType>;
    auto fixedReader = ReaderType::New();
    fixedReader->SetFileName(fixedPath);
    fixedReader->Update();
    auto movingReader = ReaderType::New();
    movingReader->SetFileName(movingPath);
    movingReader->Update();

    stpt::ImageType::Pointer fixed = fixedReader->GetOutput();
    stpt::ImageType::Pointer moving = movingReader->GetOutput();

    stpt::MaskPair masks;
    masks.fixedMask = ReadMask(fixedMaskPath);
    masks.movingMask = ReadMask(movingMaskPath);

    // Pre-flight the deformable stage's memory before doing any work. The buffers scale
    // strictly with voxel count, so a job that cannot fit is worth saying so about now
    // rather than after the linear stages have run.
    if (wantDeformable)
    {
      const auto      fixedSize = fixed->GetBufferedRegion().GetSize();
      const long long estimate = stpt::EstimateDeformablePeakBytes(static_cast<int>(fixedSize[0]),
                                                                   static_cast<int>(fixedSize[1]),
                                                                   static_cast<int>(fixedSize[2]),
                                                                   defParams);
      const long long voxels = static_cast<long long>(fixedSize[0]) * fixedSize[1] * fixedSize[2];
      // Each mask is rasterized onto the fixed grid as one byte per voxel and held for
      // the whole stage. The estimator itself only knows about image size, so the
      // driver adds what it alone can see.
      long long maskBytes = 0;
      if (!fixedMaskPath.empty()) { maskBytes += voxels; }
      if (!movingMaskPath.empty()) { maskBytes += voxels; }
      std::cout << "deformable stage: " << voxels << " voxels (" << fixedSize[0] << "x"
                << fixedSize[1] << "x" << fixedSize[2] << "), estimated peak "
                << HumanBytes(estimate + maskBytes) << std::endl;

      itksys::SystemInformation sysInfo;
      sysInfo.RunMemoryCheck();
      const long long physical = static_cast<long long>(sysInfo.GetTotalPhysicalMemory()) * 1024 * 1024;
      if (physical > 0 && estimate + maskBytes > physical)
      {
        std::cerr << std::endl;
        std::cerr << "WARNING: the estimate above exceeds this machine's physical memory ("
                  << HumanBytes(physical) << ")." << std::endl;
        std::cerr << "         Expect this run to thrash or to be killed. Memory scales with"
                  << std::endl;
        std::cerr << "         voxel count, so halving the resolution costs about 8x less."
                  << std::endl;
        std::cerr << "         Register coarse, then apply at full resolution:" << std::endl;
        std::cerr << "           gridpoint -f fixed_coarse -m moving_coarse -o out "
                  << "--write-deformation-field" << std::endl;
        std::cerr << "           antsApplyTransforms -d 3 -i moving_full -r fixed_full" << std::endl;
        std::cerr << "               -o warped.nii.gz -t out/stpt_Warp.nii.gz -t out/stpt_Affine.txt"
                  << std::endl;
        std::cerr << "         Note --grid is in VOXELS, so scale it by the same factor."
                  << std::endl;
        std::cerr << "         Proceeding anyway." << std::endl;
        std::cerr << std::endl;
      }
    }

    // Scale the optimizer to the image unless the caller asked for the upstream constants.
    // scales[translation] = 1/extent makes a unit step in the scaled parameter space
    // correspond to translating by the whole image, which behaves the same whether the data
    // is in microns or millimetres. Verified to converge on 48^3 / 96^3 / 128^3 synthetic
    // pairs, where the raw upstream values diverge on two of the three.
    if (!upstreamOptimizer)
    {
      double extent = 0.0;
      for (unsigned int i = 0; i < stpt::Dimension; ++i)
      {
        extent = std::max(extent,
                          fixed->GetBufferedRegion().GetSize()[i] * fixed->GetSpacing()[i]);
      }
      if (!userSetTranslationScale && extent > 0.0) { params.translationScale = extent; }
      if (!userSetMaxStep) { params.maximumStepLength = 0.2; }
      if (!userSetMinStep) { params.minimumStepLength = 0.0001; }
      if (params.verbose)
      {
        std::cout << "  optimizer: translation scale " << params.translationScale
                  << ", step " << params.maximumStepLength << " -> "
                  << params.minimumStepLength << " (auto from image extent)" << std::endl;
      }
    }

    // Keep the coarsest pyramid level from collapsing. The upstream schedule (4 levels,
    // coarsest shrink 8) assumes a full-size brain volume; on a small image it can leave a
    // handful of voxels per axis, which makes the Mattes metric degenerate and abort with
    // "too many samples map outside moving image buffer". Cap the levels so the coarsest
    // level keeps at least kMinCoarsestExtent voxels on its shortest axis.
    constexpr itk::SizeValueType kMinCoarsestExtent = 16;
    const auto                   fixedSize = fixed->GetBufferedRegion().GetSize();
    itk::SizeValueType           shortestAxis = fixedSize[0];
    for (unsigned int i = 1; i < stpt::Dimension; ++i)
    {
      shortestAxis = std::min(shortestAxis, fixedSize[i]);
    }
    unsigned int allowedLevels = 1;
    while (allowedLevels < params.numberOfLevels &&
           (shortestAxis >> allowedLevels) >= kMinCoarsestExtent)
    {
      ++allowedLevels;
    }
    if (allowedLevels < params.numberOfLevels)
    {
      std::cerr << "warning: reducing multi-resolution levels " << params.numberOfLevels
                << " -> " << allowedLevels << " so the coarsest level keeps at least "
                << kMinCoarsestExtent << " voxels per axis (shortest axis is "
                << shortestAxis << ")." << std::endl;
      params.numberOfLevels = allowedLevels;
    }

    stpt::RigidTransformType::Pointer  rigid;
    stpt::AffineTransformType::Pointer affine;

    if (wantRigid)
    {
      if (params.verbose) { std::cout << "stage: rigid" << std::endl; }
      rigid = stpt::RunRigidRegistration(fixed, moving, masks, params);
    }
    if (wantAffine)
    {
      if (params.verbose) { std::cout << "stage: affine" << std::endl; }
      affine = stpt::RunAffineRegistration(fixed, moving, masks, params, rigid.GetPointer());
    }

    // The affine already contains the rigid when both ran, so the final linear
    // transform is whichever stage ran last.
    // Deformable-only runs have no linear transform: resample through identity and
    // do not write a transform file that would claim otherwise.
    const bool haveLinear = (rigid.IsNotNull() || affine.IsNotNull());
    itk::Transform<double, stpt::Dimension, stpt::Dimension>::Pointer finalTransform;
    if (affine) { finalTransform = affine.GetPointer(); }
    else if (rigid) { finalTransform = rigid.GetPointer(); }
    else { finalTransform = itk::IdentityTransform<double, stpt::Dimension>::New().GetPointer(); }

    const std::string transformPath = outDir + "/stpt_Affine.txt";
    if (haveLinear)
    {
      auto transformWriter = itk::TransformFileWriterTemplate<double>::New();
      transformWriter->SetFileName(transformPath);
      transformWriter->SetInput(finalTransform);
      transformWriter->Update();
    }

    using ResampleType = itk::ResampleImageFilter<stpt::ImageType, stpt::ImageType>;
    auto resample = ResampleType::New();
    resample->SetInput(moving);
    resample->SetTransform(finalTransform);
    resample->SetInterpolator(itk::LinearInterpolateImageFunction<stpt::ImageType, double>::New());
    resample->SetOutputParametersFromImage(fixed);
    resample->SetDefaultPixelValue(0);
    resample->Update();

    // --- deformable: discrete symmetric MRF on the linearly-aligned pair ----------
    stpt::ImageType::Pointer linearWarped = resample->GetOutput();
    stpt::ImageType::Pointer finalWarped = linearWarped;
    stpt::DeformableResult   deformable;
    bool                     ranDeformable = false;

    if (wantDeformable)
    {
      defParams.verbose = params.verbose;
      if (params.verbose) { std::cout << "stage: deformable (discrete MRF)" << std::endl; }

      // Mask the control grid.
      //
      // A fixed-image mask lands on the grid directly. A moving-image mask has to be
      // carried through the linear transform first, because the deformable stage sees
      // the moving image only AFTER that transform resampled it onto the fixed grid --
      // the mask has to make the same journey as the image it describes.
      //
      // The two are intersected into one buffer here, matching the ANTs rule that a
      // sample counts only if it is valid in both images. Doing it in the driver means
      // the MRF has to know about a single mask and its hot loop is untouched.
      //
      // Nearest neighbour throughout: a mask is a label image, not an intensity.
      std::vector<unsigned char> gridMaskFixed;
      std::vector<unsigned char> gridMaskMoving;
      if (masks.fixedMask || masks.movingMask)
      {
        const itk::SizeValueType voxels = fixed->GetBufferedRegion().GetNumberOfPixels();

        auto rasterize =
          [&](const stpt::MaskSpatialObjectType *                             mask,
              const itk::Transform<double, stpt::Dimension, stpt::Dimension> * toMaskSpace,
              const char *                                                    label,
              std::vector<unsigned char> &                                    out) {
            using MaskResampleType =
              itk::ResampleImageFilter<stpt::MaskImageType, stpt::MaskImageType>;
            auto resampleMask = MaskResampleType::New();
            resampleMask->SetInput(mask->GetImage());
            resampleMask->SetInterpolator(
              itk::NearestNeighborInterpolateImageFunction<stpt::MaskImageType, double>::New());
            resampleMask->SetTransform(toMaskSpace);
            resampleMask->SetOutputParametersFromImage(fixed);
            resampleMask->SetDefaultPixelValue(0);
            resampleMask->Update();
            const unsigned char * rasterized = resampleMask->GetOutput()->GetBufferPointer();
            out.assign(rasterized, rasterized + voxels);
            if (params.verbose)
            {
              itk::SizeValueType inside = 0;
              for (itk::SizeValueType v = 0; v < voxels; ++v) { inside += (out[v] != 0); }
              std::cout << "  deformable: " << label << " mask covers "
                        << (100.0 * static_cast<double>(inside) / static_cast<double>(voxels))
                        << "% of the volume" << std::endl;
            }
          };

        if (masks.fixedMask)
        {
          rasterize(masks.fixedMask, itk::IdentityTransform<double, stpt::Dimension>::New(),
                    "fixed", gridMaskFixed);
        }
        if (masks.movingMask) { rasterize(masks.movingMask, finalTransform, "moving", gridMaskMoving); }
      }

      const auto size = fixed->GetBufferedRegion().GetSize();
      deformable = stpt::RunDeformableRegistration(fixed->GetBufferPointer(),
                                                   linearWarped->GetBufferPointer(),
                                                   static_cast<int>(size[0]),
                                                   static_cast<int>(size[1]),
                                                   static_cast<int>(size[2]),
                                                   defParams,
                                                   gridMaskFixed.empty()
                                                     ? nullptr
                                                     : gridMaskFixed.data(),
                                                   gridMaskMoving.empty()
                                                     ? nullptr
                                                     : gridMaskMoving.data());
      ranDeformable = true;

      // The algorithm already resampled the moving image; copy that out rather than
      // re-warping, so the written image is exactly what the metric was computed on.
      auto warpedImage = stpt::ImageType::New();
      warpedImage->CopyInformation(fixed);
      warpedImage->SetRegions(fixed->GetBufferedRegion());
      warpedImage->Allocate();
      std::copy(deformable.warped.begin(), deformable.warped.end(), warpedImage->GetBufferPointer());
      finalWarped = warpedImage;

      if (writeDeformationField)
      {
        using VectorType = itk::Vector<float, stpt::Dimension>;
        using FieldType = itk::Image<VectorType, stpt::Dimension>;

        // Both fields go through this one helper so their conventions cannot drift
        // apart: voxel displacements -> physical, honouring spacing and direction
        // cosines, on the fixed grid.
        auto writeField = [&](const std::vector<float> & compX,
                              const std::vector<float> & compY,
                              const std::vector<float> & compZ,
                              const std::string &        path) {
          auto field = FieldType::New();
          field->CopyInformation(fixed);
          field->SetRegions(fixed->GetBufferedRegion());
          field->Allocate();

          const auto   spacing = fixed->GetSpacing();
          const auto   direction = fixed->GetDirection();
          VectorType * out = field->GetBufferPointer();
          for (size_t i = 0; i < compX.size(); ++i)
          {
            const double local[stpt::Dimension] = { compX[i] * spacing[0],
                                                    compY[i] * spacing[1],
                                                    compZ[i] * spacing[2] };
            VectorType   v;
            for (unsigned int r = 0; r < stpt::Dimension; ++r)
            {
              double acc = 0.0;
              for (unsigned int c = 0; c < stpt::Dimension; ++c)
              {
                acc += direction[r][c] * local[c];
              }
              v[r] = static_cast<float>(acc);
            }
            out[i] = v;
          }

          auto fieldWriter = itk::ImageFileWriter<FieldType>::New();
          fieldWriter->SetFileName(path);
          fieldWriter->SetInput(field);
          fieldWriter->Update();
          if (params.verbose) { std::cout << "wrote " << path << std::endl; }
        };

        writeField(deformable.dispX, deformable.dispY, deformable.dispZ,
                   outDir + "/stpt_Warp.nii.gz");
        // The backward half of the symmetric optimization, for mapping fixed-space data
        // (atlas labels, point sets) back onto the moving image.
        writeField(deformable.invX, deformable.invY, deformable.invZ,
                   outDir + "/stpt_InverseWarp.nii.gz");
      }
    }

    const std::string warpedPath = outDir + "/stpt_Warped" + warpedExt;
    auto              writer = itk::ImageFileWriter<stpt::ImageType>::New();
    writer->SetFileName(warpedPath);
    writer->SetInput(finalWarped);
    writer->Update();

    std::ofstream log(outDir + "/stpt.log");
    log << "gridpoint\n"
        << "fixed: " << fixedPath << "\n"
        << "moving: " << movingPath << "\n"
        << "stages: " << stagesText << "\n"
        << "levels: " << params.numberOfLevels << "\n"
        << "bins: " << params.numberOfHistogramBins << "\n"
        << "spatial samples: " << params.numberOfSpatialSamples << "\n"
        << "iterations: " << params.numberOfIterations << "\n"
        << "transform: " << transformPath << "\n"
        << "warped: " << warpedPath << "\n"
        << "final parameters: " << finalTransform->GetParameters() << "\n";
    if (ranDeformable)
    {
      log << "deformable levels: " << deformable.levelsUsed << "\n"
          << "deformable grid start: " << deformable.gridStartUsed << "\n"
          << "deformable lambda: " << defParams.lambda << "\n"
          << "mean squared difference: " << deformable.metricBefore << " -> "
          << deformable.metricAfter << "\n";
    }
    log.close();

    if (params.verbose)
    {
      std::cout << "wrote " << transformPath << "\n"
                << "wrote " << warpedPath << std::endl;
    }
  }
  catch (const itk::ExceptionObject & e)
  {
    std::cerr << "error: " << e << std::endl;
    return EXIT_FAILURE;
  }
  catch (const std::exception & e)
  {
    std::cerr << "error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
