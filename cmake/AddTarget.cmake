# Generate one executable per available backend from src/<target>.{cpp,cu}.
# The host backends (CPP/TBB/OMP) share the same .cpp and differ only by the
# Hydra device system; the CUDA backend compiles the .cu with nvcc.
function(ADD_TARGET target_name)

  set(cpp_src "${CMAKE_CURRENT_SOURCE_DIR}/src/${target_name}.cpp")

  # single-thread host backend (always built)
  add_executable(${target_name}_cpp ${cpp_src})
  target_compile_definitions(${target_name}_cpp PRIVATE HYDRA_HOST_SYSTEM=CPP HYDRA_DEVICE_SYSTEM=CPP)
  target_link_libraries(${target_name}_cpp PkgConfig::LIBCONFIGXX)

  if(TBB_FOUND)
    add_executable(${target_name}_tbb ${cpp_src})
    target_compile_definitions(${target_name}_tbb PRIVATE HYDRA_HOST_SYSTEM=CPP HYDRA_DEVICE_SYSTEM=TBB)
    target_link_libraries(${target_name}_tbb PkgConfig::LIBCONFIGXX TBB::tbb)
  endif()

  if(OpenMP_CXX_FOUND)
    add_executable(${target_name}_omp ${cpp_src})
    target_compile_definitions(${target_name}_omp PRIVATE HYDRA_HOST_SYSTEM=CPP HYDRA_DEVICE_SYSTEM=OMP)
    target_link_libraries(${target_name}_omp PkgConfig::LIBCONFIGXX OpenMP::OpenMP_CXX)
  endif()

  if(CUDA_FOUND)
    cuda_add_executable(${target_name}_cuda "${CMAKE_CURRENT_SOURCE_DIR}/src/${target_name}.cu"
      OPTIONS -Xcompiler -DHYDRA_HOST_SYSTEM=CPP -DHYDRA_DEVICE_SYSTEM=CUDA)
    target_link_libraries(${target_name}_cuda PkgConfig::LIBCONFIGXX)
  endif()

endfunction()
