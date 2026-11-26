#!/bin/bash
# Build bench_collision:  ./build.sh <tag> <block size> <compiler> <flags...>
#   ./build.sh gf8 8 gfortran -O3 -fstack-arrays -march=native
# -> build/<tag>/bench_collision (and the assembly of the variants in variants.s)
set -e
here=$(cd "$(dirname "$0")" && pwd); src=$here/../../src
tag=$1; bs=$2; fc=$3; shift 3
b=$here/build/$tag; mkdir -p $b; cd $b
defs="-DAMR_BLOCK_SIZE=$bs -DAMR_L0=1024 -DAMR_MIN_LEVEL=1 -DAMR_DOMAIN_HEIGHT=3"
# the C compiler of the same family for the intrinsics kernel (override: CC_BENCH)
case $(basename $fc) in
  gfortran-*) pp="-cpp -ffree-line-length-none"; omp=-fopenmp;  cc=gcc-${fc##*-} ;;
  gfortran)   pp="-cpp -ffree-line-length-none"; omp=-fopenmp;  cc=gcc ;;
  ifx)        pp="-fpp";                         omp=-qopenmp;  cc=$(command -v icx || echo gcc) ;;
  flang*)     pp="-cpp";                         omp=-fopenmp;  cc=$(command -v clang || echo gcc) ;;
esac
cc=${CC_BENCH:-$cc}
$cc -O3 -march=native -DAMR_BLOCK_SIZE=$bs -c $here/collision_avx512.c
$cc -O3 -march=native -DAMR_BLOCK_SIZE=$bs -S -o c512.s $here/collision_avx512.c
for f in Kinds_Mod Config_Mod Blocks_Mod D2Q9_Mod Physics_Mod Model_Mod Phase_Field_Mod LBM_Common_Mod Model_JCP2016_Mod; do
  $fc $pp $omp $defs "$@" -I$src -c $src/$f.f90
done
$fc $pp $omp $defs "$@" -c $here/Model_JCP2016_Ref_Mod.f90
$fc $pp $omp $defs "$@" -I$here -c $here/Collision_Variants_Mod.F90
$fc $pp $omp $defs "$@" -I$here -S -o variants.s $here/Collision_Variants_Mod.F90
$fc $pp $omp $defs "$@" -I$here -c $here/bench_collision.F90
$fc $omp "$@" -o bench_collision *.o ${GCC_XPACK:+-Wl,-rpath,$GCC_XPACK/lib64}
echo "built $b/bench_collision"
