cmake_minimum_required(VERSION 3.16)
project(gridpoint VERSION 0.4.0 LANGUAGES C CXX)

# Bake the version and the exact commit into the binary. Without this there is no way
# to tell which build you are running: the Code Ocean postInstall clones, builds and
# then deletes the source tree, so nothing on that machine records what went in.
find_package(Git QUIET)
set(GRIDPOINT_GIT_SHA "unknown")
if(GIT_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse --short=12 HEAD
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE GRIDPOINT_GIT_SHA
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
  if(NOT GRIDPOINT_GIT_SHA)
    set(GRIDPOINT_GIT_SHA "unknown")
  endif()
endif()
message(STATUS "gridpoint ${PROJECT_VERSION} (git ${GRIDPOINT_GIT_SHA})")
# C is enabled only because ITK's HDF5 module reaches FindHDF5, which probes the
# toolchain by compiling a .c file. A CXX-only project fails that probe outright.

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ITK is the ONLY dependency: no Kakadu/JPEG2000, TinyXML, Boost or VTK.
#
# Request just the modules we use rather than every ITK module: it keeps the link
# small and avoids dragging in optional image-format backends (JPEG2000 and friends)
# that a given ITK build may not have compiled.
find_package(ITK REQUIRED COMPONENTS
  ITKCommon
  ITKRegistrationCommon
  ITKOptimizers
  ITKTransform
  ITKSpatialObjects
  ITKImageGrid
  ITKImageFunction
  ITKImageIntensity
  ITKSmoothing
  ITKIOImageBase
  ITKIONIFTI
  ITKIOTransformBase
  ITKIOTransformInsightLegacy
)
include(${ITK_USE_FILE})

# ---------------------------------------------------------------------------
# Optional CUDA backend for the deformable cost volume.
#
# Off by default: a machine with no CUDA toolkit must still build and run, and with
# this off the binary is byte-for-byte the same program it was before the backend
# existed -- the stub below supplies the same symbols and always declines.
#
# The GPU kernel is bit-identical to the CPU path, so switching backends cannot
# change a result. See the contract at the top of src/StptDeformableCuda.h; the
# --fmad=false below is load-bearing, not decorative.
# ---------------------------------------------------------------------------
option(GRIDPOINT_ENABLE_CUDA "Build the CUDA backend for the deformable cost volume" OFF)

if(GRIDPOINT_ENABLE_CUDA)
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    # Volta through Ampere as real code, plus PTX so a newer card still runs. Set
    # explicitly rather than "native" for the same reason -march=native is avoided
    # above: the Code Ocean capsule builds on one host and runs on another.
    set(CMAKE_CUDA_ARCHITECTURES 70 75 80 86 86-virtual)
  endif()
  set(GRIDPOINT_DEVICE_SOURCE src/StptDeformableCuda.cu)
  message(STATUS "gridpoint: CUDA backend ON (architectures ${CMAKE_CUDA_ARCHITECTURES})")
else()
  set(GRIDPOINT_DEVICE_SOURCE src/StptDeformableCudaStub.cxx)
  message(STATUS "gridpoint: CUDA backend OFF (-DGRIDPOINT_ENABLE_CUDA=ON to enable)")
endif()

add_executable(gridpoint
  src/gridpoint.cxx
  src/StptLinearRegistration.cxx
  src/StptDeformableMRF.cxx
  ${GRIDPOINT_DEVICE_SOURCE}
)
target_include_directories(gridpoint PRIVATE src)
target_compile_definitions(gridpoint PRIVATE
  GRIDPOINT_VERSION="${PROJECT_VERSION}"
  GRIDPOINT_GIT_SHA="${GRIDPOINT_GIT_SHA}")
target_link_libraries(gridpoint PRIVATE ${ITK_LIBRARIES})

# ---------------------------------------------------------------------------
# Vectorization.
#
# With no architecture flag, Release targets baseline x86-64 -- SSE2, from 2003 --
# so the deformable stage's dominant loop ran 4-wide at best on hardware that does
# 8. x86-64-v3 (AVX2+FMA, Haswell 2013 / Zen 2017) is the floor assumed here; it is
# set explicitly rather than via -march=native, because the Code Ocean capsule
# builds in postInstall on one host and runs on another, and native would bake in
# the builder's ISA and fault with SIGILL when the fleet rolls over.
#
# Contraction is OFF deliberately. The label loop in optimization() accumulates
# `cost1[l] += fv * mv`, and letting the compiler fuse that into an FMA changes the
# rounding -- the result would stay correct but stop being bit-identical to every
# template built before this flag existed. Vectorizing ACROSS labels is bit-exact
# on its own (each label has a private accumulator, so no sum is reassociated), so
# turning contraction off is the whole of what bit-compatibility costs: 8-wide mul
# plus 8-wide add instead of 8-wide FMA.
#
# This is not a precaution, it is a measured requirement. With contraction left on,
# the bit-exactness matrix (5 configurations, 23 outputs) drops from 23/23 to 14/23.
#
# aarch64 is left alone: NEON is baseline there, so no flag is needed.
# ---------------------------------------------------------------------------
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|x86|i686")
  if(MSVC)
    # /fp:strict, not /fp:precise. MSVC contracts a*b+c into vfmadd under
    # /fp:precise once /arch:AVX2 is on -- not only in the label loop but in the
    # trilinear interpolation inside interpolate3d/updateImage, which changes the
    # warped image and therefore every field downstream. MSVC 16.11 has no source-
    # level contraction switch (#pragma STDC FP_CONTRACT is not recognized), so the
    # only way to keep bit-compatibility is to take it off the whole target.
    #
    # /fp:strict does not cost the vectorization: MSVC still reports C5001 on the
    # hot label loop, and strict measures about 6% slower than precise on the whole
    # deformable stage. That 6% is what the bit-exactness guarantee costs.
    # Scoped to CXX: nvcc does not take MSVC switches directly, and the CUDA
    # backend has its own float-contraction control (--fmad=false).
    target_compile_options(gridpoint PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/arch:AVX2>
                                             $<$<COMPILE_LANGUAGE:CXX>:/fp:strict>)
  else()
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-march=x86-64-v3" GRIDPOINT_HAS_X86_64_V3)
    if(GRIDPOINT_HAS_X86_64_V3)
      target_compile_options(gridpoint PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=x86-64-v3>)
    else()
      # Older toolchains predate the -v3 microarchitecture levels.
      check_cxx_compiler_flag("-mavx2" GRIDPOINT_HAS_AVX2)
      if(GRIDPOINT_HAS_AVX2)
        target_compile_options(gridpoint PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-mavx2>
                                                 $<$<COMPILE_LANGUAGE:CXX>:-mfma>)
      endif()
    endif()
  endif()
endif()
if(NOT MSVC)
  target_compile_options(gridpoint PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-ffp-contract=off>)
endif()

if(GRIDPOINT_ENABLE_CUDA)
  # The exact counterpart of -ffp-contract=off / MSVC's /fp:strict, and required for
  # the same reason: nvcc fuses a*b+c into an FMA by default, which rounds once
  # instead of twice and would stop the GPU matching the CPU. Never add
  # -use_fast_math here -- it re-enables contraction AND switches division to the
  # approximate form, and the node index arithmetic divides by a float grid step.
  target_compile_options(gridpoint PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--fmad=false>)

  # On MSVC the CUDA runtime has to match the C runtime everything else uses. ITK is
  # built against the dynamic CRT, and nvcc defaults to the STATIC cudart, which drags
  # in LIBCMT alongside MSVCRT -- two C runtimes, two heaps, two sets of global state,
  # announced only as LNK4098. Elsewhere the static runtime is preferred, since it is
  # one fewer shared object to ship.
  if(MSVC)
    set_target_properties(gridpoint PROPERTIES CUDA_RUNTIME_LIBRARY Shared)
  else()
    set_target_properties(gridpoint PROPERTIES CUDA_RUNTIME_LIBRARY Static)
  endif()
endif()


install(TARGETS gridpoint RUNTIME DESTINATION bin)
