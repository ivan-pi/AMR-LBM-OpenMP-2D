#!/bin/bash
# Build bench_ghost:  ./build.sh <tag> <block size> <compiler> <flags...>
#   ./build.sh gf14_8 8 gfortran-14 -O3 -fstack-arrays -march=native
# Mesh: rising bubble with L0 = 64*block size (four levels), background level 3 (override: L0=, MINLEV=)
set -e
here=$(cd "$(dirname "$0")" && pwd); src=$here/../../src
tag=$1; bs=$2; fc=$3; shift 3
b=$here/build/$tag; mkdir -p $b; cd $b
defs="-DAMR_BLOCK_SIZE=$bs -DAMR_L0=${L0:-$((64*bs))} -DAMR_MIN_LEVEL=${MINLEV:-3} -DAMR_DOMAIN_HEIGHT=3"
case $(basename $fc) in
  gfortran*) pp="-cpp -ffree-line-length-none"; omp=-fopenmp ;;
  ifx)       pp="-fpp";                         omp=-qopenmp ;;
  flang*)    pp="-cpp";                         omp=-fopenmp ;;
esac
for f in Kinds_Mod Config_Mod Blocks_Mod D2Q9_Mod Physics_Mod Model_Mod Case_Mod Interpolations_Mod BC_Mod \
         LBM_Common_Mod Phase_Field_Mod Model_JCP2016_Mod Mesh_Mod Ghost_C_Mod Ghost_hg_Mod Prolong_Restrict_Mod \
         Refine_Mod Rising_Bubble_Mod; do
  $fc $pp $omp $defs "$@" -I$src -c $src/$f.f90
done
$fc $pp $omp $defs "$@" -I$here -I$src -c $here/Ghost_Variants.F90
$fc $pp $omp $defs "$@" -I$here -I$src -S -o variants.s $here/Ghost_Variants.F90
$fc $pp $omp $defs "$@" -I$here -c $here/bench_ghost.F90
$fc $omp "$@" -o bench_ghost *.o ${GCC_XPACK:+-Wl,-rpath,$GCC_XPACK/lib64}
echo "built $b/bench_ghost"
