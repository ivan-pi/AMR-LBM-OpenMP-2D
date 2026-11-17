# OpenMP target-offload prototype of the collision kernel (Intel ifx, SPIR-V).
# Uses the solver's Kinds_Mod and the shared build settings
# (amr_options: flags, precision, block size).
#
# Run on the GPU:   ./collision_offload [nblocks] [nrepeat]
# Run on the host:  OMP_TARGET_OFFLOAD=DISABLED ./collision_offload ...
# FP64 on Arc A7xx (emulated):
#   OverrideDefaultFP64Settings=1 IGC_EnableDPEmulation=1 ./collision_offload

add_executable(collision_offload
  ${PROJECT_SOURCE_DIR}/src/Kinds_Mod.f90
  collision_offload.f90)

set(_offload_flags -fiopenmp -fopenmp-targets=spir64)
target_compile_options(collision_offload PRIVATE ${_offload_flags})
target_link_options(collision_offload PRIVATE ${_offload_flags})
target_link_libraries(collision_offload PRIVATE amr_options)
set_target_properties(collision_offload PROPERTIES
  Fortran_PREPROCESS ON
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
