### CMakeList for NewtonianNoiseCA

cmake_minimum_required(VERSION 3.16)

# project name
project(NewtonianNoiseCA)

# warn user if system is not UNIX
if(NOT UNIX)
  message(FATAL_ERROR "This is an unsupported system.")
endif()

#cmake path dir
SET(TFBOOST_CMAKE_DIR "${PROJECT_SOURCE_DIR}/cmake")
SET(CMAKE_MODULE_PATH "${TFBOOST_CMAKE_DIR}" ${CMAKE_MODULE_PATH})
SET(CMAKE_VERBOSE_MAKEFILE  ON)

#check if compiler is C++20 compliant
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("--std=c++20" COMPILER_SUPPORTS_CXX20)
if(NOT COMPILER_SUPPORTS_CXX20)
  message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++20 support. Please use a different C++ compiler.")
endif()

#compiler flags
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

 MESSAGE(STATUS "Setting Clang flags")
 set(CMAKE_CXX_FLAGS " --std=c++20 -W -march=native -fPIC -O3 -funroll-loops -ldl" CACHE STRING "compile flags" FORCE)

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

 MESSAGE(STATUS "Setting GCC flags")
 set(CMAKE_CXX_FLAGS " --std=c++20 -W -march=native -fPIC -O3 -funroll-loops -ldl" CACHE STRING "compile flags" FORCE)

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")

 MESSAGE(STATUS "Setting ICC flags")
 set(CMAKE_CXX_FLAGS " --std=c++20 -W -march=native -fPIC -O3 -funroll-loops -ldl" CACHE STRING "compile flags" FORCE)

#elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
 #MESSAGE(STATUS "===> Compiling with MSVC")
endif()

# get Hydra (>= v4.0, header only)
find_package(Hydra REQUIRED)
include_directories(${HYDRA_INCLUDE_DIR})

# get TCLAP (header only)
find_package(TCLAP REQUIRED)
include_directories(${TCLAP_INCLUDE_DIR})

# get libconfig (used to read the .cfg configuration files)
find_package(LibConfig REQUIRED)
include_directories(${LIBCONFIG_CPP_INCLUDE_DIR})
link_directories(${LIBCONFIG_CPP_LIBRARIES})

#get CUDA
find_package(CUDA 8.0)
if(CUDA_FOUND)
  link_directories( ${CUDA_TOOLKIT_ROOT_DIR}/lib64/)

  #set cuda flags
  SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}  -ftemplate-backtrace-limit=0; --cudart;  static; -O4; --expt-relaxed-constexpr; -fmad=true; --expt-extended-lambda;--relocatable-device-code=false ;
-Xptxas -dlcm=ca; -Xptxas --opt-level=4 )

  SET(CUDA_PROPAGATE_HOST_FLAGS ON)
  SET(CUDA_SEPARABLE_COMPILATION OFF)
  SET(CUDA_VERBOSE_BUILD OFF)

  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.8)
    LIST(APPEND CUDA_NVCC_FLAGS " -Xcompiler -D_MWAITXINTRIN_H_INCLUDED ")
  endif()

  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.4)
    LIST(APPEND CUDA_NVCC_FLAGS " -Xcompiler -D__CORRECT_ISO_CPP11_MATH_H_PROTO ")
  endif()

  # Detect CUDA architecture and get best NVCC flags
  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindCudaArch.cmake)
  SELECT_NVCC_ARCH_FLAGS(NVCC_FLAGS_EXTRA)
  LIST(APPEND CUDA_NVCC_FLAGS ${NVCC_FLAGS_EXTRA})
endif(CUDA_FOUND)

#get TBB
find_package(TBB )
if(TBB_FOUND)
  include_directories(${TBB_INCLUDE_DIRS})
  link_directories(${TBB_LIBRARY})
endif(TBB_FOUND)

#get OpenMP
find_package(OpenMP)
if(OPENMP_CXX_FOUND OR OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_CXX_FOUND OR OPENMP_FOUND)

# examples etc
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#always on
SET(BUILD_CPP_TARGETS  ON)

#crucial for clang build, where openmp support is not embeded in the compiler
if( OPENMP_CXX_FOUND OR OPENMP_FOUND)
  SET(BUILD_OMP_TARGETS  ON)
else()
  SET(BUILD_OMP_TARGETS  OFF)
endif()

#not every computer has a GPU and/or cuda installed
if( CUDA_FOUND )
  SET(BUILD_CUDA_TARGETS  ON)
else()
  SET(BUILD_CUDA_TARGETS  OFF)
endif()

# to enable tbb builds
if(TBB_FOUND)
  SET(BUILD_TBB_TARGETS  ON)
else()
  SET(BUILD_TBB_TARGETS  OFF)
endif(TBB_FOUND)

# messages
MESSAGE(STATUS "TCLAP Include path: ${TCLAP_INCLUDE_DIR}")
MESSAGE(STATUS "Build CUDA/NVCC-based targets: ${BUILD_CUDA_TARGETS}")
MESSAGE(STATUS "CUDA include: ${CUDA_INCLUDE_DIRS}" )
MESSAGE(STATUS "CUDA libraries: ${CUDA_LIBRARIES}"  )
MESSAGE(STATUS "NVCC flags: ${CUDA_NVCC_FLAGS}" )
MESSAGE(STATUS "Build TBB targets: ${BUILD_TBB_TARGETS}")
MESSAGE(STATUS "TBB include: ${TBB_INCLUDE_DIRS}" )
MESSAGE(STATUS "TBB libraries: ${TBB_LIBRARIES}" )
MESSAGE(STATUS "Build OpenMP targets: ${BUILD_OMP_TARGETS}")
MESSAGE(STATUS "OpenMP flags: ${OpenMP_CXX_FLAGS}" )
MESSAGE(STATUS "OpenMP libraries:  ${OpenMP_CXX_LIBRARIES}" )
MESSAGE(STATUS "Host compiler flags: ${CMAKE_CXX_FLAGS}" )
MESSAGE(STATUS "Install diretory: ${CMAKE_INSTALL_PREFIX}"  )
MESSAGE(STATUS "Project source dir: ${PROJECT_SOURCE_DIR}" )
MESSAGE(STATUS "Project build dir: ${PROJECT_BINARY_DIR}" )

#+++++++++++++++++++++++++++
#   targets    +
#+++++++++++++++++++++++++++
include(${TFBOOST_CMAKE_DIR}/AddTarget.cmake)
include(${TFBOOST_CMAKE_DIR}/AddExample.cmake)

ADD_TARGET(Integrator_isotropic)
ADD_TARGET(Integrator_anisotropic)
