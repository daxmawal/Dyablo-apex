# set minimal version the one requested by kokkos
cmake_minimum_required(VERSION 3.16)

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.19.0")
  message(STATUS "Setting policy CMP0114 to use new behavior")
  cmake_policy(SET CMP0114 OLD)
endif()

###################################################################
# Use the CMake Superbuild pattern to help build dependencies 
# before building Dyablo.
# Uses ExternalProject_add to build dependencies 
# (Kokkos is compiled as a subdyrectory in Dyablo)
# external dependencies are also verified
# - Dyablo : MPI, HDF5
# Dependencies may be verified twice (here and when compilig Dyablo externalproject)
# The intent is to find as much external dependency errors at configure time as possible
###################################################################
project (dyablo_superbuild LANGUAGES C CXX)

include(GNUInstallDirs)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/core/cmake)
include(preventBuildInSource)

###################################################################
# Global options for Dyablo supported by Superbuild
# Additional options can be accessed from actual Dyablo project core/CMakeLists.txt
###################################################################
option (DYABLO_ENABLE_UNIT_TESTING "Enable unit testing" OFF)
option (DYABLO_USE_MPI "Activate / want MPI build" ON)
option (DYABLO_USE_HDF5 "build HDF5 input/output support" ON)
option (DYABLO_USE_MPI_CUDA_AWARE_ENFORCED "Some MPI cuda-aware implementation are not well detected; use this to enforce" OFF)
set( DYABLO_CXX_COMPILER  ${CMAKE_CXX_COMPILER} CACHE STRING "compiler used to build dyablo")
option (DYABLO_USE_INTERNAL_KOKKOS "Use internal Kokkos from external/kokkos" ON)
option (Kokkos_ENABLE_CUDA "Enable Kokkos::CUDA backend" OFF)
option (Kokkos_ENABLE_SYCL "Enable Kokkos::SYCL backend" OFF)
option (Kokkos_ENABLE_HIP "Enable Kokkos::HIP backend" OFF)
option (Kokkos_ENABLE_TUNING "Create bindings for Kokkos tuning tools" OFF)
option (DYABLO_ENABLE_APEX "Build the APEX Kokkos tool" OFF)
set(DYABLO_APEX_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/apex" CACHE PATH
    "Path to a local APEX checkout")
set(DYABLO_APEX_GIT_REPOSITORY "" CACHE STRING
    "Optional APEX Git repository to clone instead of using DYABLO_APEX_SOURCE_DIR")
set(DYABLO_APEX_GIT_TAG "" CACHE STRING
    "APEX Git tag, branch, or commit to use with DYABLO_APEX_GIT_REPOSITORY")
set(DYABLO_CUDATOOLKIT_ROOT "" CACHE PATH
    "Optional CUDA toolkit root forwarded to the Dyablo subproject")
set(DYABLO_CUDA_COMPILER "" CACHE FILEPATH
    "Optional CUDA compiler forwarded to the Dyablo subproject")
set( Kokkos_ARCH  "" CACHE STRING "List of architectures to set for Kokkos : -DKokkos_ARCH_<arch>=ON")
set( DYABLO_CMAKE_ARGS "" CACHE STRING "Additional cmake command line arguments to forward to Dyablo subproject"  )

include (ExternalProject)
set (DYABLO_EXTERNALPROJECT_DEPENDENCIES)
set (DYABLO_CMAKE_ARGS_ ${DYABLO_CMAKE_ARGS})

if( DYABLO_USE_INTERNAL_KOKKOS )
  # Read list from ${Kokkos_ARCH} and set appropriate -DKokkos_ARCH_*=ON
  message( STATUS "Forced architectures for Kokkos:" )
  foreach(ARCH ${Kokkos_ARCH})
    list (APPEND DYABLO_CMAKE_ARGS_
      -DKokkos_ARCH_${ARCH}=ON
    )
    message( STATUS "- Kokkos_ARCH_${ARCH}=ON" )
  endforeach()
  set(DYABLO_KOKKOS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/kokkos")
else( DYABLO_USE_INTERNAL_KOKKOS )
  # Find installed Kokkos 
  find_package(Kokkos 3.3 REQUIRED)

  set(DYABLO_KOKKOS_DIR ${Kokkos_DIR})
endif()

if(NOT DYABLO_CUDATOOLKIT_ROOT AND DEFINED CUDAToolkit_ROOT)
  set(DYABLO_CUDATOOLKIT_ROOT "${CUDAToolkit_ROOT}" CACHE PATH
      "Optional CUDA toolkit root forwarded to the Dyablo subproject" FORCE)
elseif(NOT DYABLO_CUDATOOLKIT_ROOT AND DEFINED ENV{CUDAToolkit_ROOT})
  set(DYABLO_CUDATOOLKIT_ROOT "$ENV{CUDAToolkit_ROOT}" CACHE PATH
      "Optional CUDA toolkit root forwarded to the Dyablo subproject" FORCE)
endif()

if(NOT DYABLO_CUDA_COMPILER AND DEFINED CMAKE_CUDA_COMPILER)
  set(DYABLO_CUDA_COMPILER "${CMAKE_CUDA_COMPILER}" CACHE FILEPATH
      "Optional CUDA compiler forwarded to the Dyablo subproject" FORCE)
elseif(NOT DYABLO_CUDA_COMPILER AND DEFINED ENV{CUDACXX})
  set(DYABLO_CUDA_COMPILER "$ENV{CUDACXX}" CACHE FILEPATH
      "Optional CUDA compiler forwarded to the Dyablo subproject" FORCE)
endif()

if(DYABLO_ENABLE_APEX)
  if(DYABLO_APEX_GIT_REPOSITORY AND NOT DYABLO_APEX_GIT_TAG)
    message(FATAL_ERROR
      "DYABLO_APEX_GIT_TAG must be set when DYABLO_APEX_GIT_REPOSITORY is used.")
  endif()
  if(NOT DYABLO_APEX_GIT_REPOSITORY AND NOT EXISTS "${DYABLO_APEX_SOURCE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR
      "DYABLO_ENABLE_APEX=ON requires either a local APEX checkout at "
      "'${DYABLO_APEX_SOURCE_DIR}' or DYABLO_APEX_GIT_REPOSITORY/DYABLO_APEX_GIT_TAG.")
  endif()
  if(NOT Kokkos_ENABLE_TUNING)
    message(STATUS "Forcing Kokkos_ENABLE_TUNING=ON because DYABLO_ENABLE_APEX is ON")
    set(Kokkos_ENABLE_TUNING ON CACHE BOOL "Create bindings for Kokkos tuning tools" FORCE)
  endif()
endif()

#####################################
# Now we build dyablo
#####################################

message( STATUS "== Setup Dyablo ExternalProject ==" )

# Check external dependencies
if( DYABLO_USE_MPI )
  find_package(MPI REQUIRED)
  list (APPEND DYABLO_CMAKE_ARGS_
    -DMPI_C_COMPILER=${MPI_C_COMPILER}
    -DMPI_CXX_COMPILER=${MPI_CXX_COMPILER}
  )
endif()
# hdf5
if (DYABLO_USE_HDF5)
  set( HDF5_PREFER_PARALLEL ON )
  set( HDF5_ENABLE_PARALLEL ON ) # This is a workaround for a bug in CMake + HDF5 2.0
  find_package( HDF5 REQUIRED COMPONENTS C)
  if( DYABLO_USE_MPI AND NOT HDF5_IS_PARALLEL )
    message( FATAL_ERROR "HDF5 library must be parallel when compiling Dyablo with MPI" )
  endif()
  list (APPEND DYABLO_CMAKE_ARGS_
    -DHDF5_C_COMPILER_EXECUTABLE=${HDF5_C_COMPILER_EXECUTABLE}
  )
endif()

# default arguments - common to all platform
list (APPEND DYABLO_CMAKE_ARGS_
  -DCMAKE_CXX_COMPILER=${DYABLO_CXX_COMPILER}
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  -DDYABLO_ENABLE_UNIT_TESTING=${DYABLO_ENABLE_UNIT_TESTING}
  -DDYABLO_USE_MPI=${DYABLO_USE_MPI}
  -DDYABLO_USE_HDF5=${DYABLO_USE_HDF5}
  -DDYABLO_USE_MPI_CUDA_AWARE_ENFORCED=${DYABLO_USE_MPI_CUDA_AWARE_ENFORCED}
  -DDYABLO_USE_INTERNAL_KOKKOS=${DYABLO_USE_INTERNAL_KOKKOS}
  -DDYABLO_ENABLE_APEX=${DYABLO_ENABLE_APEX}
  -DDYABLO_APEX_SOURCE_DIR=${DYABLO_APEX_SOURCE_DIR}
  -DDYABLO_APEX_GIT_REPOSITORY=${DYABLO_APEX_GIT_REPOSITORY}
  -DDYABLO_APEX_GIT_TAG=${DYABLO_APEX_GIT_TAG}
  -DKokkos_ENABLE_CUDA=${Kokkos_ENABLE_CUDA}
  -DKokkos_ENABLE_SYCL=${Kokkos_ENABLE_SYCL}
  -DKokkos_ENABLE_HIP=${Kokkos_ENABLE_HIP}
  -DKokkos_ENABLE_TUNING=${Kokkos_ENABLE_TUNING}
  -DKokkos_DIR=${DYABLO_KOKKOS_DIR}
)

if(DYABLO_CUDATOOLKIT_ROOT)
  list(APPEND DYABLO_CMAKE_ARGS_
    -DCUDAToolkit_ROOT=${DYABLO_CUDATOOLKIT_ROOT}
  )
endif()

if(DYABLO_CUDA_COMPILER)
  list(APPEND DYABLO_CMAKE_ARGS_
    -DCMAKE_CUDA_COMPILER=${DYABLO_CUDA_COMPILER}
  )
endif()

set(DYABLO_SUBPROJECT_ENV)
if(DYABLO_CUDATOOLKIT_ROOT)
  list(APPEND DYABLO_SUBPROJECT_ENV
    CUDA_ROOT=${DYABLO_CUDATOOLKIT_ROOT}
  )
endif()
if(DYABLO_CUDA_COMPILER)
  get_filename_component(DYABLO_CUDA_COMPILER_BIN_DIR "${DYABLO_CUDA_COMPILER}" DIRECTORY)
  list(APPEND DYABLO_SUBPROJECT_ENV
    CUDACXX=${DYABLO_CUDA_COMPILER}
    PATH=${DYABLO_CUDA_COMPILER_BIN_DIR}:$ENV{PATH}
  )
endif()

# Build dyablo
ExternalProject_Add (dyablo
  DEPENDS ${DYABLO_EXTERNALPROJECT_DEPENDENCIES}
  BUILD_ALWAYS TRUE
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/core
  CMAKE_ARGS ${DYABLO_CMAKE_ARGS_}
  BUILD_COMMAND ${CMAKE_COMMAND} -E env ${DYABLO_SUBPROJECT_ENV} ${CMAKE_COMMAND} --build .
  BUILD_JOB_SERVER_AWARE TRUE
  INSTALL_COMMAND ""
  BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/dyablo)

# add target dyablo-test
if (DYABLO_ENABLE_UNIT_TESTING)
  ExternalProject_Add_Step (dyablo test
    COMMAND           make test
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dyablo
    COMMENT           "running dyablo unit tests"
    ALWAYS            TRUE
    EXCLUDE_FROM_MAIN TRUE
    )
  ExternalProject_Add_StepTargets(dyablo test)
endif()

message(  ${DYABLO_CMAKE_ARGS_} )

message( WARNING "This is not an error, just a reminder: 
CMake superbuild is meant to help you setup dyablo by compiling dependencies before configuring the actual Dyablo cmake project.
CMake variables are not forwarded to the actual Dyablo CMake project : if you wish to setup dyablo with custom variables (such as CMAKE_CXX_FLAGS), you need to modify the CMakeCache (with cmake-gui or ccmake for instance) from subdirectory ${CMAKE_CURRENT_BINARY_DIR}/dyablo" )

if( CMAKE_CXX_FLAGS )
  message( FATAL_ERROR "CMAKE_CXX_FLAGS is not forwarded to Dyablo subproject : see warning above" )
endif()
