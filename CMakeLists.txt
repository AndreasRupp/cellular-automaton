cmake_minimum_required(VERSION 3.16)

project(CAM LANGUAGES CXX)

## Check whether appropriate compiler version is used
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10)
    message(
      FATAL_ERROR
      "Insufficient gcc/g++ version! Minimum version required is 10.
       Consider setting variable CXX in Makefile or CMAKE_CXX_COMPILER in CMake, respectively."
    )
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 10)
    message(
      FATAL_ERROR
      "Insufficient clang/clang++ version! Minimum version required is 10.
       Consider setting variable CXX in Makefile or CMAKE_CXX_COMPILER in CMake, respectively."
    )
  endif()
endif()


## Create a target that transports all the necessary information to downstream projects
add_library(CAM INTERFACE)

## Add the correct include directories
target_include_directories(
  CAM
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

## Add the necessary C++ standard requirement
target_compile_features(
  CAM
  INTERFACE
  cxx_std_20
)

## Add the examples subdirectory
add_subdirectory(examples_cpp)

## Add a namespaced alias target for the export
add_library(CAM::CAM ALIAS CAM)

## Install targets, headers and configuration
include(GNUInstallDirs)
install(
  TARGETS CAM
  EXPORT CAM-config
)

install(
  EXPORT CAM-config
  NAMESPACE CAM::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/CAM
)

install(
  DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

## Write a configuration file for Cython
configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/mex/CAM/run_cam.m.in
  ${CMAKE_CURRENT_BINARY_DIR}/run_cam.m
)

## Write a configuration file for Cython
configure_file(
  ${CMAKE_CURRENT_LIST_DIR}/mex/CAM/run_cam.cxx.in
  ${CMAKE_CURRENT_BINARY_DIR}/run_cam.cxx.in
)

## Print a summary
include(FeatureSummary)
feature_summary(WHAT ALL)
