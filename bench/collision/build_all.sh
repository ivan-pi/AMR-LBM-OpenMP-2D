#!/bin/bash
# Build the configurations of the study for a block size:  ./build_all.sh 8
# tags: <compiler><block>[_avx2|_512]; default = 256-bit vectors with the AVX-512 register set
here=$(cd "$(dirname "$0")" && pwd); bs=${1:-8}
G="-O3 -fstack-arrays"
for gf in gfortran-14:gf14_ gfortran:gf15_; do
  fc=${gf%%:*}; t=${gf##*:}$bs
  command -v $fc >/dev/null || continue
  $here/build.sh ${t}_avx2 $bs $fc $G -march=x86-64-v3
  $here/build.sh ${t}      $bs $fc $G -march=native
  $here/build.sh ${t}_512  $bs $fc $G -march=native -mprefer-vector-width=512
done
if command -v ifx >/dev/null; then
  $here/build.sh ifx${bs}_avx2 $bs ifx -O3 -xCORE-AVX2
  $here/build.sh ifx${bs}      $bs ifx -O3 -xHost
  $here/build.sh ifx${bs}_512  $bs ifx -O3 -xHost -qopt-zmm-usage=high
fi
if command -v flang >/dev/null; then
  $here/build.sh fl${bs}_avx2 $bs flang -O3 -march=x86-64-v3
  $here/build.sh fl${bs}      $bs flang -O3 -march=native
  $here/build.sh fl${bs}_512  $bs flang -O3 -march=native -mllvm -force-vector-width=8
fi
