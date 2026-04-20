# set minimal version the one requested by kokkos
cmake_minimum_required(VERSION 3.16)

project(dyablo C CXX)

# default local cmake macro repository
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
# Prevent from build in source tree
include(preventBuildInSource)
# Write a small header with build data, git version, etc...
include(write_version)
# Init build type: Release, Debug, ...
include(initBuildType)
# always export compile commands database
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(FetchContent)
include(ExternalProject)
include(GNUInstallDirs)

# C++14 is for Kokkos
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

# dyablo options from Superbuild
option (DYABLO_ENABLE_UNIT_TESTING "Enable unit testing" OFF)
option (DYABLO_USE_MPI "Activate / want MPI build" ON)
option (DYABLO_USE_HDF5 "build HDF5 input/output support" ON)
option (DYABLO_USE_MPI_CUDA_AWARE_ENFORCED "Some MPI cuda-aware implementation are not well detected; use this to enforce" OFF)
option (DYABLO_USE_INTERNAL_KOKKOS "Use internal Kokkos from external/kokkos" ON)
option (DYABLO_ENABLE_APEX "Build the APEX Kokkos tool" OFF)
set(DYABLO_APEX_SOURCE_DIR "${PROJECT_SOURCE_DIR}/../external/apex" CACHE PATH
    "Path to a local APEX checkout")
set(DYABLO_APEX_GIT_REPOSITORY "" CACHE STRING
    "Optional APEX Git repository to clone instead of using DYABLO_APEX_SOURCE_DIR")
set(DYABLO_APEX_GIT_TAG "" CACHE STRING
    "APEX Git tag, branch, or commit to use with DYABLO_APEX_GIT_REPOSITORY")
option (Kokkos_ENABLE_CUDA   "enable Kokkos::Cuda backend" OFF)
option (Kokkos_ENABLE_HIP   "enable Kokkos::HIP backend" OFF)
option (Kokkos_ENABLE_SYCL   "enable Kokkos::SyCL backend" OFF)
option (Kokkos_ENABLE_TUNING "Create bindings for Kokkos tuning tools" OFF)
option (Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC "Use CudaMallocAsync for Kokkos allocations" OFF)
set    (Kokkos_DIR "Kokkos installation directory or internal source" CACHE PATH ""  )


# Options not forwarded from superbuild
option (DYABLO_BUILD_DOC  "Enable / disable documentation build" OFF)
set    (DYABLO_DOC "doxygen" CACHE STRING "documentation type (doxygen or mkdocs)")
option (DYABLO_USE_DOUBLE "build with double precision" ON)
option (DYABLO_USE_FPE_DEBUG "build with floating point Nan tracing (signal handler)" OFF)
option (Kokkos_ENABLE_OPENMP "enable Kokkos::OpenMP backend" ON)
set (Kokkos_ENABLE_CUDA_LAMBDA ${Kokkos_ENABLE_CUDA} CACHE BOOL "Enable use of cuda lambdas" FORCE)
set (Kokkos_ENABLE_CUDA_CONSTEXPR ${Kokkos_ENABLE_CUDA} CACHE BOOL "Enable use of cuda constexpr" FORCE)

if(DYABLO_ENABLE_APEX AND NOT Kokkos_ENABLE_TUNING)
  message(STATUS "Forcing Kokkos_ENABLE_TUNING=ON because DYABLO_ENABLE_APEX is ON")
  set(Kokkos_ENABLE_TUNING ON CACHE BOOL "Create bindings for Kokkos tuning tools" FORCE)
endif()

# Use CMAKE_DL_LIBS to find the dl library
set( LIBDL_LIBRARY CACHE STRING "LIBDL library" "${CMAKE_DL_LIBS}" )

#####################################################################
# Kokkos
#####################################################################
if( DYABLO_USE_INTERNAL_KOKKOS )
  FetchContent_Declare( kokkos_external
    SOURCE_DIR ${Kokkos_DIR}
  )
  # Import kokkos targets
  FetchContent_MakeAvailable(kokkos_external)
else( DYABLO_USE_INTERNAL_KOKKOS )
  find_package( Kokkos 3.5 REQUIRED )
  kokkos_check( DEVICES "OpenMP" )
  if(Kokkos_ENABLE_CUDA)
    # kokkos_check is defined in KokkosConfigCommon.cmake
    kokkos_check( DEVICES "Cuda" )
    kokkos_check( OPTIONS CUDA_LAMBDA)
  endif()
endif()

set(DYABLO_APEX_CUDA_CMAKE_ARGS)
if(Kokkos_ENABLE_CUDA)
  find_package(CUDAToolkit REQUIRED)
  get_filename_component(DYABLO_CUDATOOLKIT_ROOT "${CUDAToolkit_BIN_DIR}/.." ABSOLUTE)
  set(DYABLO_APEX_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  list(APPEND DYABLO_APEX_CUDA_CMAKE_ARGS
    -DCUDAToolkit_ROOT=${DYABLO_CUDATOOLKIT_ROOT}
  )
  if(CMAKE_CUDA_COMPILER)
    get_filename_component(DYABLO_CUDA_COMPILER_BIN_DIR "${CMAKE_CUDA_COMPILER}" DIRECTORY)
    list(APPEND DYABLO_APEX_CUDA_CMAKE_ARGS
      -DCMAKE_CUDA_COMPILER=${CMAKE_CUDA_COMPILER}
    )
  endif()
  find_path(DYABLO_NVTX_INCLUDE_DIR
    NAMES nvToolsExt.h
    HINTS
      "${DYABLO_CUDATOOLKIT_ROOT}/include"
      "${DYABLO_CUDATOOLKIT_ROOT}/include/nvtx3"
      "${DYABLO_CUDA_COMPILER_BIN_DIR}/../include-stdpar/nvtx3"
  )
  if(DYABLO_NVTX_INCLUDE_DIR)
    string(APPEND DYABLO_APEX_CXX_FLAGS " -I${DYABLO_NVTX_INCLUDE_DIR}")
    list(APPEND DYABLO_APEX_CUDA_CMAKE_ARGS
      -DCMAKE_CXX_FLAGS=${DYABLO_APEX_CXX_FLAGS}
    )
    message(STATUS "Using NVTX headers from ${DYABLO_NVTX_INCLUDE_DIR} for APEX")
  endif()
endif()

# Check compiler version
message( "Detected kokkos compiler : ${KOKKOS_CXX_COMPILER_ID} ${KOKKOS_CXX_COMPILER_VERSION}" )
set( NVCC_FIRST_VERSION 11.2 )
if(   KOKKOS_CXX_COMPILER_ID STREQUAL "NVIDIA" AND KOKKOS_CXX_COMPILER_VERSION LESS NVCC_FIRST_VERSION )
  set( DYABLO_ENABLE_DEPRECATED_NVCC "Use deprecated NVCC" CACHE STRING "I want right results"  )
  if( DYABLO_ENABLE_DEPRECATED_NVCC STREQUAL "I want wrong results" )
    message( WARNING "nvcc version ${KOKKOS_CXX_COMPILER_VERSION} (< ${NVCC_FIRST_VERSION}) is deprecated but you chose to use it anyway. That may cause wrong numerical results." )
  else()
    message( FATAL_ERROR "nvcc version ${KOKKOS_CXX_COMPILER_VERSION} (< ${NVCC_FIRST_VERSION}) is deprecated because it may cause wrong numerical results. If you want to use it anyway, set DYABLO_ENABLE_DEPRECATED_NVCC to \"I want wrong results\" " )
  endif()
endif()

#####################################################################
# MPI
#####################################################################
if(DYABLO_USE_MPI)
  find_package(MPI REQUIRED)

  message(STATUS "MPI support found")
  message(STATUS "MPI compile flags: " ${MPI_CXX_COMPILE_FLAGS})
  message(STATUS "MPI include path: " ${MPI_CXX_INCLUDE_PATH})
  message(STATUS "MPI LINK flags path: " ${MPI_CXX_LINK_FLAGS})
  message(STATUS "MPI libraries: " ${MPI_CXX_LIBRARIES})

  if( DYABLO_USE_MPI_CUDA_AWARE_ENFORCED )
    message( STATUS "  MPI is CUDA Aware because DYABLO_USE_MPI_CUDA_AWARE_ENFORCED is ON" )
    set( MPI_IS_CUDA_AWARE ON )
  else( DYABLO_USE_MPI_CUDA_AWARE_ENFORCED )
    if(Kokkos_ENABLE_CUDA)
      find_program(OMPI_INFO
        NAMES ompi_info
        PATHS ${MPI_CXX_LIBRARIES}/../bin)

      message( STATUS "MPI checking CUDA support..." )
      # Full command line to probe if cuda support in MPI implementation is enabled
      # ompi_info --parsable --all | grep mpi_built_with_cuda_support:value
      if (OMPI_INFO)
        message( STATUS "  Found ompi_info" )
        execute_process(COMMAND ${OMPI_INFO}
          OUTPUT_VARIABLE _output)
        if ((_output MATCHES "smcuda"))
          message( STATUS "  OpenMPI is CUDA Aware" )
          set( MPI_IS_CUDA_AWARE ON )
        else((_output MATCHES "smcuda"))
          message( STATUS "  OpenMPI is NOT CUDA Aware" )
          set( MPI_IS_CUDA_AWARE OFF )
        endif((_output MATCHES "smcuda"))
      else(OMPI_INFO)
        message( STATUS "  MPI is supposed NOT CUDA Aware : could not detect (To override, set DYABLO_USE_MPI_CUDA_AWARE_ENFORCED) " )
        set( MPI_IS_CUDA_AWARE OFF )
      endif(OMPI_INFO)
    else(Kokkos_ENABLE_CUDA)
      message( STATUS "  MPI is supposed NOT CUDA Aware : CUDA not used (To override, set DYABLO_USE_MPI_CUDA_AWARE_ENFORCED) " )
      set( MPI_IS_CUDA_AWARE OFF )
    endif(Kokkos_ENABLE_CUDA)
  endif( DYABLO_USE_MPI_CUDA_AWARE_ENFORCED )

  message( STATUS "  MPI CUDA Aware = ${MPI_IS_CUDA_AWARE}" )

  if( MPI_IS_CUDA_AWARE )
    add_compile_options(-DMPI_IS_CUDA_AWARE)
  endif ( MPI_IS_CUDA_AWARE )
endif(DYABLO_USE_MPI)

#####################################################################
# HDF5
#####################################################################
# prefer using parallel HDF5 when build with mpi
if (DYABLO_USE_MPI)
  set(HDF5_PREFER_PARALLEL TRUE)
  set(HDF5_ENABLE_PARALLEL ON) # This is a workaround for a bug in CMake + HDF5 2.0
endif(DYABLO_USE_MPI)

if (DYABLO_USE_HDF5)
  find_package(HDF5)
  if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
    set(my_hdf5_libs hdf5 hdf5_cpp)
    add_compile_options(-DDYABLO_USE_HDF5)
    if (HDF5_IS_PARALLEL)
      add_compile_options(-DDYABLO_USE_HDF5_PARALLEL)
    endif()
  endif(HDF5_FOUND)
endif(DYABLO_USE_HDF5)

#####################################################################
# Backward-cpp
#####################################################################
FetchContent_Declare( backward-cpp_external
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/../external/backward-cpp
)
# Import kokkos targets
FetchContent_MakeAvailable(backward-cpp_external)

set(DYABLO_APEX_LIBRARY_PATH "")
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

  set(DYABLO_APEX_INSTALL_DIR "${CMAKE_BINARY_DIR}/apex/install")
  set(DYABLO_APEX_LIBRARY_PATH
      "${DYABLO_APEX_INSTALL_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}apex${CMAKE_SHARED_LIBRARY_SUFFIX}")
  set(DYABLO_APEX_KOKKOS_DIR "${CMAKE_BINARY_DIR}/cmake_packages/Kokkos")
  set(DYABLO_APEX_EXTERNALPROJECT_SOURCE_ARGS)
  if(DYABLO_APEX_GIT_REPOSITORY)
    message(STATUS
      "Fetching APEX from ${DYABLO_APEX_GIT_REPOSITORY} (${DYABLO_APEX_GIT_TAG})")
    list(APPEND DYABLO_APEX_EXTERNALPROJECT_SOURCE_ARGS
      GIT_REPOSITORY ${DYABLO_APEX_GIT_REPOSITORY}
      GIT_TAG ${DYABLO_APEX_GIT_TAG}
      GIT_SHALLOW TRUE
    )
  else()
    message(STATUS "Using local APEX checkout from ${DYABLO_APEX_SOURCE_DIR}")
    list(APPEND DYABLO_APEX_EXTERNALPROJECT_SOURCE_ARGS
      SOURCE_DIR ${DYABLO_APEX_SOURCE_DIR}
    )
  endif()

  ExternalProject_Add(apex
    ${DYABLO_APEX_EXTERNALPROJECT_SOURCE_ARGS}
    BINARY_DIR ${CMAKE_BINARY_DIR}/apex/build
    CMAKE_ARGS
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_INSTALL_PREFIX=${DYABLO_APEX_INSTALL_DIR}
      -DCMAKE_INSTALL_LIBDIR=lib
      -DAPEX_BUILD_TESTS=OFF
      -DAPEX_BUILD_EXAMPLES=OFF
      -DAPEX_WITH_KOKKOS=ON
      -DAPEX_WITH_MPI=OFF
      -DAPEX_WITH_OMPT=OFF
      -DAPEX_WITH_CUDA=${Kokkos_ENABLE_CUDA}
      -DAPEX_WITH_HIP=${Kokkos_ENABLE_HIP}
      -DAPEX_WITH_LEVEL0=OFF
      -DAPEX_WITH_PLUGINS=OFF
      -DKokkos_DIR=${DYABLO_APEX_KOKKOS_DIR}
      ${DYABLO_APEX_CUDA_CMAKE_ARGS}
    BUILD_BYPRODUCTS ${DYABLO_APEX_LIBRARY_PATH}
    INSTALL_DIR ${DYABLO_APEX_INSTALL_DIR}
    INSTALL_COMMAND ${CMAKE_COMMAND} --build . --target install
  )
endif()

# common flags
if (Kokkos_ENABLE_CUDA)
  add_compile_options(-DCUDA)
endif()
if (DYABLO_USE_DOUBLE)
  add_compile_options(-DDYABLO_USE_DOUBLE)
endif()
if (DYABLO_USE_MPI)
  add_compile_options(-DDYABLO_USE_MPI)    
endif()
if (DYABLO_USE_FPE_DEBUG)
  add_compile_options(-DDYABLO_USE_FPE_DEBUG)
endif()

##
## Using flags -Wextra, it's to strong for Kokkos, too many warnings
## But -Wall is really a minimum
##
#add_definitions(-Wall -Wextra -Wconversion -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wcast-qual -Wunused -Woverloaded-virtual  -pedantic )
#add_definitions( -Wall -Wextra )
add_definitions( -Wall )

add_subdirectory(bin)
add_subdirectory(src)

# Documentation type
if(DYABLO_BUILD_DOC)
  add_subdirectory(doc)
endif(DYABLO_BUILD_DOC)

# unit testing
if (DYABLO_ENABLE_UNIT_TESTING)
  enable_testing()
  add_subdirectory(unit_test)
endif()

##################### PRINT CONFIGURE STATUS ######################
message("//===================================================")
message("  ${PROJECT_NAME} build configuration:")
message("//===================================================")
message("  C++ Compiler : ${CMAKE_CXX_COMPILER_ID} " 
  "${CMAKE_CXX_COMPILER_VERSION} " 
  "${CMAKE_CXX_COMPILER_WRAPPER}")
message("    ${CMAKE_CXX_COMPILER}")
if (DYABLO_USE_MPI)
  message("  MPI enabled")
else(DYABLO_USE_MPI)
  message("  MPI not enabled")
endif(DYABLO_USE_MPI)
if (Kokkos_FOUND)
  message("  Kokkos_Found true")
  message("  Kokkos_ENABLE_HWLOC=${Kokkos_ENABLE_HWLOC}")
  message("  Kokkos_ENABLE_OPENMP=${Kokkos_ENABLE_OPENMP}")
  message("  Kokkos_ENABLE_DEBUG=${Kokkos_ENABLE_DEBUG}")
  message("  Kokkos_ENABLE_DEBUG_BOUNDS_CHECK=${Kokkos_ENABLE_DEBUG_BOUNDS_CHECK}")
  message("  Kokkos_ENABLE_CUDA=${Kokkos_ENABLE_CUDA}")
  message("  Kokkos_ENABLE_CUDA_LAMBDA=${Kokkos_ENABLE_CUDA_LAMBDA}")
  message("  Kokkos_ENABLE_CUDA_CONSTEXPR=${Kokkos_ENABLE_CUDA_CONSTEXPR}")
  message("  Kokkos_ENABLE_CUDA_UVM=${Kokkos_ENABLE_CUDA_UVM}")
  message("  Kokkos_ENABLE_SYCL=${Kokkos_ENABLE_SYCL}")
  message("  Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=${Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC}")
  message("  Kokkos_ENABLE_HIP=${Kokkos_ENABLE_HIP}")
  message("  Kokkos_ENABLE_TUNING=${Kokkos_ENABLE_TUNING}")
  message("  Kokkos_ARCH=${Kokkos_ARCH}")
endif(Kokkos_FOUND)
if (DYABLO_ENABLE_APEX)
  message("  APEX enabled")
  message("  APEX runtime library    : ${DYABLO_APEX_LIBRARY_PATH}")
endif (DYABLO_ENABLE_APEX)
if (DYABLO_USE_HDF5)
  if (HDF5_FOUND)
    message("  HDF5 found version    : ${HDF5_VERSION}")
    message("  HDF5 definitions      : ${HDF5_DEFINITIONS}")
    message("  HDF5 parallel         : ${HDF5_IS_PARALLEL}")
    message("  HDF5 includes dirs    : ${HDF5_INCLUDE_DIRS}")
    message("  HDF5 libraries        : ${HDF5_LIBRARIES}")
  endif(HDF5_FOUND)
endif(DYABLO_USE_HDF5)
message("")
if (DYABLO_BUILD_DOC)
  message("   Building doc of type ${DYABLO_DOC}")
endif (DYABLO_BUILD_DOC)
message("")
