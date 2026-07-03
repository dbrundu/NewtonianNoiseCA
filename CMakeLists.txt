### CMakeLists for NewtonianNoiseCA

cmake_minimum_required(VERSION 3.16)
project(NewtonianNoiseCA CXX)

if(NOT UNIX)
  message(FATAL_ERROR "This is an unsupported system.")
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# C++20, optimized host build
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-O3 -march=native -funroll-loops -fPIC)

#--- required dependencies (all header only except libconfig) ----------------
find_package(Hydra REQUIRED)                         # custom module, see cmake/
include_directories(${HYDRA_INCLUDE_DIR})

find_path(TCLAP_INCLUDE_DIR NAMES tclap/CmdLine.h REQUIRED)
include_directories(${TCLAP_INCLUDE_DIR})

find_package(PkgConfig REQUIRED)
pkg_check_modules(LIBCONFIGXX REQUIRED IMPORTED_TARGET libconfig++)

#--- optional parallel backends ----------------------------------------------
find_package(TBB CONFIG)                             # provides TBB::tbb, TBB_FOUND
find_package(OpenMP)                                 # provides OpenMP::OpenMP_CXX

find_package(CUDA 8.0)
if(CUDA_FOUND)
  link_directories(${CUDA_TOOLKIT_ROOT_DIR}/lib64/)
  set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -O3 --expt-relaxed-constexpr --expt-extended-lambda -fmad=true)
  set(CUDA_PROPAGATE_HOST_FLAGS ON)
  include(FindCudaArch)
  SELECT_NVCC_ARCH_FLAGS(NVCC_FLAGS_EXTRA)
  list(APPEND CUDA_NVCC_FLAGS ${NVCC_FLAGS_EXTRA})
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

message(STATUS "Backends -> CPP: ON  TBB: ${TBB_FOUND}  OpenMP: ${OpenMP_CXX_FOUND}  CUDA: ${CUDA_FOUND}")

#--- targets -----------------------------------------------------------------
include(AddTarget)
ADD_TARGET(Integrator_isotropic)
ADD_TARGET(Integrator_anisotropic)
