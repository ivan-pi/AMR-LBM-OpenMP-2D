# amrviz with fields in GPU memory: SYCL and OpenMP target-offload demos
# (Intel oneAPI; tested on an Arc A770)
#
#   source /opt/intel/oneapi/setvars.sh
#   cmake -S viewer/examples -B build-viewer-examples -DCMAKE_CXX_COMPILER=icpx -DCMAKE_Fortran_COMPILER=ifx
#   cmake --build build-viewer-examples -j
#   ./build-viewer-examples/sycl_demo         # or ./build-viewer-examples/offload_demo
cmake_minimum_required(VERSION 3.20)
project(amrviz_examples LANGUAGES C CXX Fortran)

find_package(SDL2 REQUIRED)
find_package(Threads REQUIRED)

add_library(amrviz STATIC ../amrviz.c ../amrviz.f90)
set_target_properties(amrviz PROPERTIES C_STANDARD 11 C_STANDARD_REQUIRED ON
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(amrviz PUBLIC .. ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_link_libraries(amrviz PRIVATE SDL2::SDL2 Threads::Threads m)

add_executable(sycl_demo sycl_demo.cpp)
target_compile_options(sycl_demo PRIVATE -fsycl)
target_link_options(sycl_demo PRIVATE -fsycl)
target_link_libraries(sycl_demo PRIVATE amrviz)

add_executable(offload_demo offload_demo.f90)
target_compile_options(offload_demo PRIVATE -fiopenmp -fopenmp-targets=spir64)
target_link_options(offload_demo PRIVATE -fiopenmp -fopenmp-targets=spir64)
target_link_libraries(offload_demo PRIVATE amrviz)
