cmake_minimum_required(VERSION 3.20)

project(AMR_LBM_2D
  DESCRIPTION "2D AMR lattice Boltzmann solver for a rising bubble (OpenMP)"
  LANGUAGES Fortran)

# ---------------------------------------------------------------------------
# Options
#
#   cmake -B build                                   # defaults (Release, double)
#   cmake -B build -DCMAKE_Fortran_COMPILER=ifx      # Intel ifx
#   cmake -B build -DCMAKE_Fortran_COMPILER=nagfor   # NAG Fortran
#   cmake -B build -DAMR_PRECISION=4                 # single precision
#   cmake -B build -DAMR_L0=1024 -DAMR_MIN_LEVEL=5   # bigger case
#   cmake -B build -DAMR_BLOCK_SIZE=16               # 16x16-cell blocks
#   cmake -B build -DAMR_DOMAIN_HEIGHT=4             # L0 x 4*L0 domain (Rayleigh-Taylor)
#   cmake -B build -DAMR_WITH_FVTKHDF=ON -DCMAKE_PREFIX_PATH=<fvtkhdf prefix>
#   cmake -B build -DAMR_BUILD_OFFLOAD=ON            # GPU prototype (ifx only)
#   cmake --build build -j

set(AMR_PRECISION 8 CACHE STRING "Working precision wp in bytes: 8 (real64) or 4 (real32)")
set_property(CACHE AMR_PRECISION PROPERTY STRINGS 8 4)
set(AMR_L0 256 CACHE STRING "Domain width in finest-level lattice units (domain is L0 x AMR_DOMAIN_HEIGHT*L0)")
set(AMR_DOMAIN_HEIGHT 3 CACHE STRING "Domain height in base blocks of L0 x L0 (rising bubble: 3, Rayleigh-Taylor: 4)")
set(AMR_MIN_LEVEL 1 CACHE STRING "Minimum refinement level (uniform background refinement)")
set(AMR_BLOCK_SIZE 8 CACHE STRING "Cells per block in each direction (power of two dividing AMR_L0)")
option(AMR_PERIODIC_Y "Periodic in y as well (default: walls at the bottom and top)" OFF)
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
if(NOT AMR_DOMAIN_HEIGHT MATCHES "^[1-9][0-9]*$")
  message(FATAL_ERROR "AMR_DOMAIN_HEIGHT must be a positive integer")
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
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
  set(AMR_FLAGS_RELEASE -O4)
  set(AMR_FLAGS_DEBUG   -O0 -g -gline -C=all)
  set(AMR_FLAG_NATIVE   -target=native)
endif()

# OpenMP. For NAG, FindOpenMP's imported target lists nagfor's thread-safe
# runtime objects but omits -openmp at link time, so nagfor also links its
# serial runtime and the link fails with duplicate symbols; let the nagfor
# driver add the OpenMP runtime itself.
add_library(amr_openmp INTERFACE)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
  target_compile_options(amr_openmp INTERFACE ${OpenMP_Fortran_FLAGS})
  target_link_options(amr_openmp INTERFACE ${OpenMP_Fortran_FLAGS})
else()
  target_link_libraries(amr_openmp INTERFACE OpenMP::OpenMP_Fortran)
endif()

add_library(amr_options INTERFACE)
target_compile_options(amr_options INTERFACE
  $<$<CONFIG:Release,RelWithDebInfo>:${AMR_FLAGS_RELEASE}>
  $<$<CONFIG:Debug>:${AMR_FLAGS_DEBUG}>
  $<$<AND:$<BOOL:${AMR_NATIVE}>,$<NOT:$<CONFIG:Debug>>>:${AMR_FLAG_NATIVE}>)

# ---------------------------------------------------------------------------
# Fypp (https://fypp.readthedocs.io, "pip install fypp"): the sources src/*.fypp
# are turned into plain Fortran in the build tree. Fypp carries the
# configuration (block size, domain, precision, optional fVTKHDF) and
# generates the kernels that exist in several forms (strip shapes, numbers of
# ghost layers, direction lists). No source goes through the C preprocessor.

find_program(FYPP_EXECUTABLE fypp REQUIRED
  DOC "The Fypp preprocessor (pip install fypp)")

set(AMR_FYPP_FLAGS
  -DAMR_BLOCK_SIZE=${AMR_BLOCK_SIZE} -DAMR_L0=${AMR_L0} -DAMR_MIN_LEVEL=${AMR_MIN_LEVEL}
  -DAMR_DOMAIN_HEIGHT=${AMR_DOMAIN_HEIGHT})
if(AMR_PRECISION EQUAL 4)
  list(APPEND AMR_FYPP_FLAGS -DWP_SINGLE=True)
endif()
if(AMR_PERIODIC_Y)
  list(APPEND AMR_FYPP_FLAGS -DAMR_PERIODIC_Y=True)
endif()
if(AMR_WITH_FVTKHDF)
  list(APPEND AMR_FYPP_FLAGS -DUSE_FVTKHDF=True)
endif()

# amr_fypp(<variable> <module>...): src/<module>.fypp -> fypp/<module>.f90 in the
# current binary directory; <variable> receives the generated files
function(amr_fypp out_var)
  set(_generated)
  foreach(_module ${ARGN})
    set(_in  ${PROJECT_SOURCE_DIR}/src/${_module}.fypp)
    set(_out ${CMAKE_CURRENT_BINARY_DIR}/fypp/${_module}.f90)
    add_custom_command(
      OUTPUT  ${_out}
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/fypp
      COMMAND ${FYPP_EXECUTABLE} ${AMR_FYPP_FLAGS} -I ${PROJECT_SOURCE_DIR}/src ${_in} ${_out}
      DEPENDS ${_in} ${PROJECT_SOURCE_DIR}/src/Biquadratic.fypp
      COMMENT "Fypp ${_module}.fypp"
      VERBATIM)
    list(APPEND _generated ${_out})
  endforeach()
  set(${out_var} ${_generated} PARENT_SCOPE)
endfunction()

amr_fypp(AMR_FYPP_SOURCES
  Kinds_Mod Config_Mod Model_JCP2016_Mod Ghost_C_Mod Ghost_hg_Mod VTK_Output_Mod Output_Mod)

# ---------------------------------------------------------------------------
# Solver: one module per file in src/, file name = module name (.f90, or
# .fypp for the modules that Fypp generates code for).
# amr_core holds everything but the driver programs; each driver selects
# an LBM model for the rising-bubble case.

add_library(amr_core OBJECT
  ${AMR_FYPP_SOURCES}
  src/Blocks_Mod.f90
  src/D2Q9_Mod.f90
  src/Physics_Mod.f90
  src/Model_Mod.f90
  src/Case_Mod.f90
  src/Interpolations_Mod.f90
  src/BC_Mod.f90
  src/LBM_Common_Mod.f90
  src/Phase_Field_Mod.f90
  src/Model_PRE2017_Mod.f90
  src/Model_Parker2009_Mod.f90
  src/Mesh_Mod.f90
  src/Prolong_Restrict_Mod.f90
  src/Refine_Mod.f90
  src/TecPlot_Mod.f90
  src/Solver_Mod.f90
  src/Rising_Bubble_Mod.f90
  src/Rayleigh_Taylor_Mod.f90
  src/Parker_Cases_Mod.f90
)
target_link_libraries(amr_core PUBLIC amr_options amr_openmp)
set_target_properties(amr_core PROPERTIES
  Fortran_PREPROCESS OFF
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(amr_core PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/mod)

if(AMR_WITH_FVTKHDF)
  enable_language(C)   # fVTKHDF's package config looks up HDF5, which tests with C
  find_package(fVTKHDF REQUIRED)
  target_link_libraries(amr_core PUBLIC fVTKHDF::fvtkhdf)
endif()

# drivers, one per case and model:
#   amr_lbm,    amr_lbm_pre2017     rising bubble
#   amr_lbm_rt, amr_lbm_rt_pre2017  Rayleigh-Taylor instability
#   amr_lbm_parker                  two species (Parker & Jovanovic 2009): droplets, spinodal decomposition
# with the models of Fakhari, Geier & Lee (JCP 2016) and Fakhari et al.
# (PRE 2017, improved locality)
add_executable(amr_lbm            src/main_jcp2016.f90)
add_executable(amr_lbm_pre2017    src/main_pre2017.f90)
add_executable(amr_lbm_rt         src/main_rt_jcp2016.f90)
add_executable(amr_lbm_rt_pre2017 src/main_rt_pre2017.f90)
add_executable(amr_lbm_parker     src/main_parker2009.f90)
set(_amr_drivers amr_lbm amr_lbm_pre2017 amr_lbm_rt amr_lbm_rt_pre2017 amr_lbm_parker)
foreach(_exe ${_amr_drivers})
  target_link_libraries(${_exe} PRIVATE amr_core)
  set_target_properties(${_exe} PROPERTIES Fortran_PREPROCESS OFF)
endforeach()

install(TARGETS ${_amr_drivers})

# ---------------------------------------------------------------------------
# GPU offload prototype

if(AMR_BUILD_OFFLOAD)
  if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    message(FATAL_ERROR "AMR_BUILD_OFFLOAD requires Intel ifx (IntelLLVM)")
  endif()
  add_subdirectory(offload)
endif()

message(STATUS "AMR-LBM: precision=${AMR_PRECISION} bytes, L0=${AMR_L0}, height=${AMR_DOMAIN_HEIGHT}*L0, periodic y=${AMR_PERIODIC_Y}, "
               "min level=${AMR_MIN_LEVEL}, block=${AMR_BLOCK_SIZE}, fVTKHDF=${AMR_WITH_FVTKHDF}, "
               "offload=${AMR_BUILD_OFFLOAD}, build type=${CMAKE_BUILD_TYPE}")
