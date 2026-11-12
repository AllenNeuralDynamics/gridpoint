cmake_minimum_required(VERSION 3.16)
project(gridpoint VERSION 0.3.0 LANGUAGES C CXX)

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

add_executable(gridpoint
  src/gridpoint.cxx
  src/StptLinearRegistration.cxx
  src/StptDeformableMRF.cxx
)
target_include_directories(gridpoint PRIVATE src)
target_compile_definitions(gridpoint PRIVATE
  GRIDPOINT_VERSION="${PROJECT_VERSION}"
  GRIDPOINT_GIT_SHA="${GRIDPOINT_GIT_SHA}")
target_link_libraries(gridpoint PRIVATE ${ITK_LIBRARIES})

install(TARGETS gridpoint RUNTIME DESTINATION bin)
