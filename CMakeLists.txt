cmake_minimum_required(VERSION 3.16)
project(gridpoint CXX)

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
target_link_libraries(gridpoint PRIVATE ${ITK_LIBRARIES})

install(TARGETS gridpoint RUNTIME DESTINATION bin)
