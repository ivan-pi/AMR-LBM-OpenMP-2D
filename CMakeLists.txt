cmake_minimum_required(VERSION 3.20)

project(AMR_LBM_2D
  DESCRIPTION "2D AMR lattice Boltzmann solver for a rising bubble (OpenMP)"
  LANGUAGES Fortran)

# ---------------------------------------------------------------------------
# Options
#
#   cmake -B build                                   # defaults (Release, double)
#   cmake -B build -DCMAKE_Fortran_COMPILER=ifx      # Intel ifx
#   cmake -B build -DAMR_PRECISION=4                 # single precision
#   cmake -B build -DAMR_L0=1024 -DAMR_MIN_LEVEL=5   # bigger case
#   cmake -B build -DAMR_BLOCK_SIZE=16               # 16x16-cell blocks
#   cmake -B build -DAMR_WITH_FVTKHDF=ON -DCMAKE_PREFIX_PATH=<fvtkhdf prefix>
#   cmake -B build -DAMR_BUILD_OFFLOAD=ON            # GPU prototype (ifx only)
#   cmake --build build -j

set(AMR_PRECISION 8 CACHE STRING "Working precision wp in bytes: 8 (real64) or 4 (real32)")
set_property(CACHE AMR_PRECISION PROPERTY STRINGS 8 4)
set(AMR_L0 256 CACHE STRING "Domain width in finest-level lattice units (domain is L0 x 3*L0)")
set(AMR_MIN_LEVEL 1 CACHE STRING "Minimum refinement level (uniform background refinement)")
set(AMR_BLOCK_SIZE 8 CACHE STRING "Cells per block in each direction (power of two dividing AMR_L0)")
option(AMR_WITH_FVTKHDF "VTKHDF output through the fVTKHDF library" OFF)
option(AMR_BUILD_OFFLOAD "Build the OpenMP target-offload prototype in offload/ (Intel ifx)" OFF)
option(AMR_NATIVE "Optimize for the build machine (-march=native / -xHost)" ON)

if(NOT AMR_PRECISION MATCHES "^(4|8)$")
  message(FATAL_ERROR "AMR_PRECISION must be 4 or 8")
endif()

# The finest level is the L with AMR_BLOCK_SIZE * 2**L == AMR_L0 (none: L0 is
# not a power-of-two multiple of the block size); AMR_MIN_LEVEL must lie in
# [1, finest level]
unset(_amr_max_level)
foreach(_lv RANGE 30)
  math(EXPR _amr_width "${AMR_BLOCK_SIZE} << ${_lv}")
  if(_amr_width EQUAL AMR_L0)
    set(_amr_max_level ${_lv})
    break()
  endif()
endforeach()
if(AMR_BLOCK_SIZE LESS 4 OR NOT DEFINED _amr_max_level)
  message(FATAL_ERROR "AMR_L0/AMR_BLOCK_SIZE must be a power of two (block size >= 4)")
endif()
if(AMR_MIN_LEVEL LESS 1 OR AMR_MIN_LEVEL GREATER _amr_max_level)
  message(FATAL_ERROR "AMR_MIN_LEVEL must be between 1 and ${_amr_max_level} for this L0/AMR_BLOCK_SIZE")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo)
endif()

find_package(OpenMP REQUIRED COMPONENTS Fortran)

# ---------------------------------------------------------------------------
# Build settings shared by all executables (flags, precision, block size)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # -fstack-arrays: keep the per-block work arrays of the kernels off the heap
  set(AMR_FLAGS_RELEASE -O3 -fstack-arrays)
  set(AMR_FLAGS_DEBUG   -O0 -g -fcheck=all -fbacktrace -Wall)
  set(AMR_FLAG_NATIVE   -march=native)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  set(AMR_FLAGS_RELEASE -O3)
  set(AMR_FLAGS_DEBUG   -O0 -g -check all -traceback -warn all)
  set(AMR_FLAG_NATIVE   -xHost)
endif()

add_library(amr_options INTERFACE)
target_compile_definitions(amr_options INTERFACE
  AMR_BLOCK_SIZE=${AMR_BLOCK_SIZE}
  $<$<STREQUAL:${AMR_PRECISION},4>:WP_SINGLE>)
target_compile_options(amr_options INTERFACE
  $<$<CONFIG:Release,RelWithDebInfo>:${AMR_FLAGS_RELEASE}>
  $<$<CONFIG:Debug>:${AMR_FLAGS_DEBUG}>
  $<$<AND:$<BOOL:${AMR_NATIVE}>,$<NOT:$<CONFIG:Debug>>>:${AMR_FLAG_NATIVE}>)

# ---------------------------------------------------------------------------
# Solver: one module per file in src/, file name = module name.
# amr_core holds everything but the driver programs; each driver selects
# an LBM model for the rising-bubble case.

add_library(amr_core OBJECT
  src/Kinds_Mod.f90
  src/Config_Mod.f90
  src/Blocks_Mod.f90
  src/D2Q9_Mod.f90
  src/Physics_Mod.f90
  src/Model_Mod.f90
  src/Case_Mod.f90
  src/Interpolations_Mod.f90
  src/BC_Mod.f90
  src/LBM_Common_Mod.f90
  src/Phase_Field_Mod.f90
  src/Model_JCP2016_Mod.f90
  src/Model_PRE2017_Mod.f90
  src/Mesh_Mod.f90
  src/Ghost_C_Mod.f90
  src/Ghost_hg_Mod.f90
  src/Prolong_Restrict_Mod.f90
  src/Refine_Mod.f90
  src/TecPlot_Mod.f90
  src/VTK_Output_Mod.f90
  src/Output_Mod.f90
  src/Solver_Mod.f90
  src/Rising_Bubble_Mod.f90
)
target_compile_definitions(amr_core PRIVATE AMR_L0=${AMR_L0} AMR_MIN_LEVEL=${AMR_MIN_LEVEL})
target_link_libraries(amr_core PUBLIC amr_options OpenMP::OpenMP_Fortran)
set_target_properties(amr_core PROPERTIES
  Fortran_PREPROCESS ON
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(amr_core PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/mod)

if(AMR_WITH_FVTKHDF)
  enable_language(C)   # fVTKHDF's package config looks up HDF5, which tests with C
  find_package(fVTKHDF REQUIRED)
  target_link_libraries(amr_core PUBLIC fVTKHDF::fvtkhdf)
  target_compile_definitions(amr_core PRIVATE USE_FVTKHDF)
endif()

# drivers: amr_lbm (Fakhari, Geier & Lee, JCP 2016) and
#          amr_lbm_pre2017 (Fakhari et al., PRE 2017, improved locality)
add_executable(amr_lbm         src/main_jcp2016.f90)
add_executable(amr_lbm_pre2017 src/main_pre2017.f90)
foreach(_exe amr_lbm amr_lbm_pre2017)
  target_link_libraries(${_exe} PRIVATE amr_core)
  set_target_properties(${_exe} PROPERTIES Fortran_PREPROCESS ON)
endforeach()

install(TARGETS amr_lbm amr_lbm_pre2017)

# ---------------------------------------------------------------------------
# GPU offload prototype

if(AMR_BUILD_OFFLOAD)
  if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    message(FATAL_ERROR "AMR_BUILD_OFFLOAD requires Intel ifx (IntelLLVM)")
  endif()
  add_subdirectory(offload)
endif()

message(STATUS "AMR-LBM: precision=${AMR_PRECISION} bytes, L0=${AMR_L0}, "
               "min level=${AMR_MIN_LEVEL}, block=${AMR_BLOCK_SIZE}, fVTKHDF=${AMR_WITH_FVTKHDF}, "
               "offload=${AMR_BUILD_OFFLOAD}, build type=${CMAKE_BUILD_TYPE}")
