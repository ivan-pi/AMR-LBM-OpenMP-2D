# Performance study of the JCP 2016 collision kernel

`Collision_h_g_AMR` of `src/Model_JCP2016_Mod.f90` (with the moments of g and
the chemical potential computed before it) is the most expensive stage of a
time step. This directory holds the benchmark used to study it, the variants
that were tried, and the measurements. The outcome is the kernel now in
`src/Collision_JCP2016.inc`.

**Result in the solver** (rising bubble, `-DAMR_L0=512 -DAMR_MIN_LEVEL=3`, 3000
steps, 1 thread, ~470 leaf blocks of 8x8 = 8 MiB of block data, which is in L3;
stage "moments of g + collision", MLUPS = million cell updates per second):

| build (CMake defaults unless noted) | before | after | speed-up |
|---|---|---|---|
| gfortran 14.2, 8x8 | 49.0 | 63.7 | 1.30 |
| gfortran 15.2, 8x8 | 42.8 | 55.3 | 1.29 |
| ifx 2026.2, 8x8 | 47.9 | 64.9 | 1.35 |
| gfortran 14.2, 8x8, `-mprefer-vector-width=512` | 59.6 | 77.3 | 1.30 |
| ifx 2026.2, 8x8, `-qopt-zmm-usage=high` | 56.3 | 75.7 | 1.35 |
| gfortran 14.2, 16x16 | 48.9 | 72.4 | 1.48 |
| ifx 2026.2, 16x16 | 54.5 | 73.0 | 1.34 |

The mesh (blocks per level at every output) is identical before and after and
the diagnostics (`bubble_rise.dat`, density range, maximum velocity, mass)
agree to the printed 6 digits. The kernel itself agrees with the old one to
1e-14 relative (round-off; not bitwise). Round-off grows slowly with the flow:
after 20000 steps the bubble centroid agrees to 4e-6 and the rise velocity to
2e-4.

Machine: Intel Xeon (Emerald Rapids, family 6 model 207, AVX-512 with two
512-bit FMA units, 48 KiB L1d, 2 MiB L2 per core), 4 cores of a virtual machine,
~4.2 GHz with and without 512-bit code. Timings in a VM are noisy: all numbers
are the best of 11 samples and differences below ~4 % are not significant. No
hardware counters were available, so statements about the pipeline come from
the assembly and from probe variants, not from `perf`.

## The benchmark

`bench_collision` fills a pool of blocks (the solver's own `Blocks_Mod`) with a
synthetic two-phase state: a lattice of bubbles with about half of the cells
inside an interface, a vortex lattice, h and g 1 % off equilibrium. It then
times sweeps of a block kernel over all blocks, every variant starting from the
same saved state, samples interleaved between the variants. The first collision
of every variant is compared with the kernel as it was before the study
(`Model_JCP2016_Ref_Mod.f90`, a frozen copy): `err` is the largest difference
of the interior of h, g (and P, Ux, Uy) over the largest magnitude of the field.

Two kinds of kernels are timed, as the solver has two:

* `collide`: chemical potential + collision (`model%Collide`, used on the
  steps with diagnostics, when P and U are already known);
* `mom+col` (names `mc_*`): chemical potential + moments of g + collision
  (`model%Moments_and_Collide`, **every other step: this is the hot path**).

Two pool sizes: 1 MiB (in L2: the kernel alone) and 12 MiB (in L3 on this
machine and on a desktop CPU with 16 MiB of L3: as in the solver).

```bash
./build.sh gf8 8 gfortran -O3 -fstack-arrays -march=native   # tag, block size, compiler, flags
./build_all.sh 8                                             # all configurations of the study
./run.sh 1 gf14_8 gf14_8_512 ifx8                            # pool in MiB, tags: table of best MLUPS
BENCH_ONLY=new,mc_new,mc_c512 OMP_NUM_THREADS=4 ./run.sh 12 gf14_8
python3 asmstat.py build/gf14_8/variants.s k_sym             # instructions, spills, divisions per kernel
```

The variants are assembled by `gen_variants.py` from the `*.inc` fragments
(the cell body with the array accesses as macros) into
`Collision_Variants_Mod.F90`. Configurations (tags): compiler + block size,
then `_avx2` (AVX2 only: 16 vector registers), nothing (`-march=native` /
`-xHost`: 256-bit vectors, 32 registers) or `_512` (512-bit vectors:
`-mprefer-vector-width=512`, `-qopt-zmm-usage=high`, flang
`-mllvm -force-vector-width=8`; flang ignores `-mprefer-vector-width`).
The complete tables are in `results/study/`.

## Experiments

MLUPS, 1 thread, 8x8 blocks, 1 MiB pool unless noted.

### 1. Loop structure: nested, collapsed, flattened

| variant | gf14 | gf14_512 | gf15 | ifx | ifx_512 | flang |
|---|---|---|---|---|---|---|
| `2d_simd`: `DO j` / `!$OMP SIMD` `DO i` (the original) | 73.5 | 103.6 | 66.8 | 81.4 | 106.0 | 79.4 |
| `2d_collapse`: `!$OMP SIMD COLLAPSE(2)` | 73.6 | 104.2 | 67.0 | 42.9 | 50.0 | 43.3 |
| `1d_divmod`: one loop over the cells, `i`, `j` from the index | 39.4 | 55.0 | 35.0 | 43.2 | 51.4 | 43.2 |
| `1d_linear`: rank-1 dummies, index linearized by hand, nested loops | 73.7 | 105.1 | 67.1 | 83.5 | 103.9 | 77.3 |

Flattening cannot pay off here and is easy to get badly wrong. The row
length (8 or 16) is a multiple of the vector length, so the nested loop has no
remainder iterations: one or two full vectors per row. A flat loop over the
cells of a block needs either `i`, `j` recomputed from the flat index
(`1d_divmod`), or a collapsed loop; in both cases ifx and flang no longer see
unit-stride accesses (the rows are separated by ghost cells) and emit
**gathers**, at half the speed. gfortran recovers the nest from `COLLAPSE(2)`
but not from the div/mod form. A truly contiguous flat loop over the ghost
columns as well is not possible with this data layout: C has two ghost layers,
h and g one, and P, Ux, Uy none, so there is no common stride.

### 2. Array bounds, sequence association, compile-time sizes

`dflt_bounds` (dummies `C(nx+4,ny+4)`, `h(nx+2,ny+2,9)` instead of
`C(-1:nx+2,-1:ny+2)`, `h(0:nx+1,0:ny+1,9)`) and `1d_linear` generate the
same code as the original with all three compilers (identical instruction
counts with gfortran): custom lower bounds are free when they are constants.
What matters is that the sizes are compile-time constants: `rt_size` (the
same kernel with `nx`, `ny` as dummy arguments) is 6-20 % slower (gf14 68.9
vs 73.5, gf14_512 83.2 vs 103.6, ifx 73.5 vs 81.4, flang 68.9 vs 79.4).

### 3. `!$OMP SIMD`

`2d_nosimd` and `sym_nosimd` are within the noise of the versions with the
directive for all compilers and widths: every compiler vectorizes the loop on
its own (the loop is branch-free, the arguments cannot alias). The directive
is kept: it documents the intent and costs nothing. `COLLAPSE(2)` hurts (1.).

### 4. Symmetry of the MRT collision

The original already paired opposite directions. `sym` goes further: all
moments but the two stresses relax with rate 1, so the post-collision value is
`geq + forcing +- m` with `m` the relaxed stress moment of `g - geq`, and
only the *sums of geq over a pair* enter `m`. In these sums the odd parts
cancel, and so do the pressure term and the constant part of the equilibrium
between the two pairs of a moment: what is left are closed expressions in u,
grad C and the buoyancy (see the header of `src/Collision_JCP2016.inc`). The
individual `geq_k` and `g_k - geq_k` are never formed, the mixed and central
differences combine into one biased difference per direction, and the force
term folds into the coefficient of Gamma.

| | gf14 | gf14_512 | gf15 | ifx | ifx_512 | flang |
|---|---|---|---|---|---|---|
| `2d_simd` | 73.5 | 103.6 | 66.8 | 81.4 | 106.0 | 79.4 |
| `sym` | 92.7 | 133.5 | 83.0 | 105.8 | 133.2 | 95.7 |
| gain | +26 % | +29 % | +24 % | +30 % | +26 % | +21 % |

Vector instructions per pass of the loop (gfortran 15, 256-bit): 506 -> 378.

### 5. Loop fusion and fission

| variant | gf14 | gf14_512 | gf15 | ifx | ifx_512 | flang |
|---|---|---|---|---|---|---|
| `sym` (mu from a separate pass) | 92.7 | 133.5 | 83.0 | 105.8 | 133.2 | 95.7 |
| `sym_fuse_mu`: mu in the cell loop | 99.4 | 141.5 | 87.2 | 110.3 | 134.8 | 105.3 |
| `sym_fis_hg`: an h pass and a g pass over the block | 82.9 | 122.4 | 78.9 | 92.0 | 96.5 | 82.2 |
| `sym_fis_hg_row`: an h loop and a g loop per row | 82.2 | 116.7 | 79.6 | 88.4 | 119.5 | 80.9 |
| `pre_row`: square root and divisions of a row first, into row buffers | 90.6 | 126.3 | 83.6 | 94.8 | 129.8 | 92.6 |
| `mc_model`: mu pass, moments pass (+ gradient pass), collision (the original) | 57.6 | 80.0 | 54.8 | 62.7 | 73.6 | 59.1 |
| `mc_sym_split`: the same three passes with the `sym` collision | 68.1 | 93.3 | 63.4 | 75.9 | 86.4 | 69.1 |
| `mc_sym_fused`: moments in the cell loop | 79.8 | 107.8 | 69.9 | 90.1 | 116.5 | 81.5 |
| `mc_sym_fused_mu`: mu and moments in the cell loop | 85.1 | 113.6 | 71.1 | 90.3 | 124.2 | 87.0 |

Fusion wins, fission loses, every time. The fused loop shares the nine loads
of C and the gradient between the chemical potential, the moments and the
collision, loads g once, and needs no work arrays. Splitting h from g
recomputes the shared part (gradient, Gamma) and is 7-25 % slower. Computing
the slow operations first (`pre_*`, see 8.) does not help either.

### 6. Register pressure

Vector stores to / loads from the stack in the kernel (whole routine), and MLUPS:

| | AVX2, 16 registers | AVX-512VL, 32 registers, 256-bit | 512-bit |
|---|---|---|---|
| gfortran 15 `2d_simd` | 114 / 146, 55.0 | 72 / 70, 66.8 | 60 / 61, 82.5 |
| gfortran 15 `sym_fuse_mu` | 92 / 98, 69.6 | 52 / 49, 87.2 | 37 / 39, 124.9 |
| gfortran 14 `sym_fuse_mu` | 78.4 | 34 / 30, 99.4 | 20 / 21, 141.5 |
| ifx `sym` | 89.8 | 0 / 61, 105.8 | 0 / 70, 133.2 |
| flang `sym_fuse_mu` | 86.2 | 26 / 38, 105.3 | 30 / 30, 120.8 |

* The 32 registers of AVX-512VL are worth 10-25 % at the same 256-bit width
  (the `_avx2` column against the default one). `-march=native` on an
  AVX-512 CPU gives them even when the vectors stay 256 bits wide.
* ifx never stores a vector to the stack; its "loads" are constants used as
  memory operands. gfortran also spills per-block constants and intermediate
  values.
* Deliberately *raising* the pressure (`sym_late`: all eight g of a cell held
  to the end of the body, as the original did) changes nothing measurable; the
  compilers reorder the stores anyway. The fused moments (five more live
  values) cost nothing either. What reduces the spills is fewer live
  temporaries by algebra (4.), not statement order.

### 7. Block size

| mom+col | gf14 | gf14_512 | gf15 | ifx | ifx_512 | flang |
|---|---|---|---|---|---|---|
| 8x8, 1 MiB: original / new | 58.7 / 84.8 | 79.8 / 119.2 | 55.8 / 69.2 | 66.2 / 95.9 | 75.6 / 122.4 | 58.8 / 87.7 |
| 16x16, 1 MiB | 61.0 / 90.2 | 71.8 / 112.7 | 58.1 / 77.1 | 70.0 / 91.8 | 81.1 / 122.7 | 67.6 / 89.5 |
| 8x8, 12 MiB | 50.0 / 66.5 | 63.8 / 86.4 | 46.1 / 58.3 | 53.5 / 73.1 | 58.5 / 80.8 | 52.2 / 72.4 |
| 16x16, 12 MiB | 53.3 / 81.7 | 59.6 / 87.2 | 49.5 / 67.2 | 57.9 / 80.3 | 67.4 / 91.5 | 59.3 / 84.4 |

In L2 the block size hardly matters for the kernel (per-block overhead is
small: the loops are fully specialized for the size). Once the data comes from
L3, 16x16 blocks are 10-20 % faster: a block carries (nx+2)^2/nx^2 = 1.56
times its interior in h and g for nx = 8 but only 1.27 times for nx = 16, so
fewer bytes are moved per cell. (For the solver as a whole a larger block
also means more cells on fine levels; that trade-off is not part of this
study.)

A layout with rows padded to a multiple of 8 values and aligned to 64 bytes
(`al_*`, so that no 512-bit load or store straddles two cache lines) gains
2-5 % in L2 and loses up to 25 % in L3 for 8x8 blocks, because the padded
block is 60 % larger. Not pursued.

### 8. Divisions and the square root

A cell needs `1/|grad C|` (interface normal), `1/rho` (velocity) and
`(2-t)/(2+t)` (relaxation rate): three divisions and a square root. Probe
variants with these replaced by a multiply-add (`x_nodiv`, wrong results)
give the cost of the divider:

| mom+col | gf14 | gf14_512 | ifx | ifx_512 |
|---|---|---|---|---|
| `mc_sym_fused_mu` | 86.0 | 115.3 | 95.6 | 125.1 |
| `mc_x_nodiv` (no division at all) | 87.6 | 135.1 | 98.7 | 132.8 |

With 256-bit vectors the divisions cost at most 2-3 %: they are not the
bottleneck. With 512-bit vectors they cost 6-17 % (one 512-bit division
occupies the divider for about 16 cycles).

Can they be avoided? Not exactly. The normal is zero only where `|grad C|` is
below 1e-16, which after a few steps is nowhere: in a run of the rising bubble
all leaf cells have a non-zero normal from the first output on (the bulk
carries noise of 1e-6 in C). `1/rho` and the relaxation rate are needed in
every cell anyway. A physical threshold (blocks with `|grad C|` < 1e-3: all
blocks below the finest level and a quarter of the finest) would skip the
normal there, but it changes the model beyond round-off, and at 256 bits
there is nothing to gain. Other attempts:

* one division for all three quotients (`*_1div`: `x = 1/(a b c)`, then
  products): 7-14 % *slower*, as the three quotients then wait for one
  another; gfortran with AVX2 does not even vectorize it;
* square root and divisions of a row or block first, into buffers (`pre_row`,
  `pre_blk`): 3-20 % slower;
* `-ffast-math -mrecip=all` (gfortran), `-no-prec-div -no-prec-sqrt` (ifx):
  no change; neither compiler uses the reciprocal estimates for doubles here.

### 9. AVX-512 intrinsics

`collision_avx512.c` is the final kernel written by hand: one row of 8 cells
per vector, every multiply-add an explicit FMA, zero-masking for the normal,
instantiated with macros for moments or not and for exact divisions (`c512`)
or the 14-bit reciprocal estimates `rcp14`/`rsqrt14` with two Newton steps
(`c512_rcp`, ~1 ulp). It is exported under C names and called through
`BIND(C)` interfaces. Compiled with the C compiler of the same family (gcc 14,
gcc 15, icx 2026.2, clang 18), `-O3 -march=native`.

| 8x8, 1 MiB | gf14 (256-bit) | gf14_512 | ifx (256-bit) | ifx_512 | flang_512 |
|---|---|---|---|---|---|
| `new` (Fortran) | 97.5 | 140.5 | 113.8 | 134.0 | 119.9 |
| `c512` | 136.5 | 130.0 | 123.5 | 120.8 | 116.5 |
| `c512_rcp` | 148.4 | 144.3 | 146.4 | 133.3 | 140.0 |
| `mc_new` (Fortran) | 84.8 | 119.2 | 95.9 | 122.4 | 110.2 |
| `mc_c512` | 119.0 | 117.2 | 122.4 | 113.3 | 111.4 |
| `mc_c512_rcp` | 127.0 | 128.8 | 125.8 | 118.6 | 122.7 |

(The C kernel is the same in both columns of a compiler; the spread is noise.)

* Against the Fortran kernel compiled for 512-bit vectors, intrinsics with
  exact divisions are **not faster** (-7 % to +2 %): the compilers' 512-bit
  code is as good as the hand-written one. gcc even spills more in the C
  version (38 vector stores against 20).
* What intrinsics add is what Fortran cannot express: the reciprocal
  estimates. They recover most of the cost of the 512-bit divider (8.):
  +3 to +10 % over the best Fortran build, results still within 1e-14.
* Against the *default* Fortran build (256-bit vectors) the C kernel is 25-45 %
  faster, but so is the Fortran kernel with the 512-bit flag (+25-40 %). On
  this CPU (two 512-bit FMA units, no frequency drop) the flag is the cheaper
  way to get there. On CPUs with one 512-bit FMA unit or a frequency penalty
  (the Rocket Lake of the main README) neither pays off.
* With 16x16 blocks the gcc builds of the C kernel reach 126 (`mc_c512`) and
  139 (`mc_c512_rcp`) against 113 for Fortran at 512 bits.
* Once the data comes from L3 (12 MiB) everything at 512 bits converges to
  80-87 MLUPS for 8x8 blocks: memory, not arithmetic, sets the pace.

The intrinsics kernel stays in the benchmark; the solver uses the Fortran one.

### 10. Threads

4 threads, 8x8, `mom+col`: 4 MiB pool (L2) gf14 212.6 -> 313.6, gf14_512
290.5 -> 436.2, ifx 236.6 -> 342.6, ifx_512 270.3 -> 441.9; 48 MiB pool gf14
186.2 -> 257.9, ifx 206.4 -> 283.5. The gain of the new kernel holds
(1.35-1.6 in L2, 1.37-1.4 beyond), scaling is 3.5-3.7 on 4 cores.

### 11. The moments kernel on its own (`Macroscopic_Values_g_AMR`)

Since the moments are fused into the collision loop, this kernel is off the
hot path: `Moments_g_Block` only runs on the steps with diagnostics (about 4 %
of the block updates) and for newly refined blocks. In the 3000-step run above
it takes 0.024 s of 5.1 s (0.5 %); before, it ran in every step. The same
levers applied to it (kind `moments`, variants `mg_*`, MLUPS, 1 thread):

| | gf14 | gf14_512 | gf15 | ifx | ifx_512 | flang |
|---|---|---|---|---|---|---|
| 8x8, 1 MiB: `mg_model` (mu pass, gradient pass, moments loop) | 235.5 | 260.8 | 293.0 | 237.6 | 239.6 | 213.2 |
| `mg_fuse_grad`: gradient in the loop | 265.6 | 293.0 | 303.5 | 301.6 | 397.9 | 286.1 |
| `mg_fused`: mu and gradient in the loop, one pass | 342.7 | 326.5 | 343.0 | 296.6 | 512.0 | 406.2 |
| `mg_fused_nosimd` | 340.1 | 330.7 | 344.2 | 331.3 | 413.9 | 409.2 |
| `mg_fused_sym`: pair sums and differences shared | 337.3 | 331.2 | 334.8 | 300.3 | 507.6 | 413.4 |
| `mg_x_nodiv`: no division (wrong results) | 341.0 | 328.1 | 340.3 | 307.1 | 405.3 | 403.8 |
| 8x8, 12 MiB: `mg_model` | 182.6 | 185.2 | 168.3 | 143.9 | 160.0 | 144.0 |
| `mg_fused` | 194.1 | 190.7 | 195.1 | 176.1 | 193.2 | 189.7 |

16x16, 12 MiB: `mg_model` 161.4 (gf14), 190.0 (ifx); `mg_fused` 229.9, 229.8.

* Fusion is again the only lever: one pass instead of three is 1.2-2.1 times
  faster in L2 and 1.05-1.3 times (8x8) to 1.2-1.4 times (16x16) from L3.
* Everything else is neutral: the directive, the shared pair sums (the
  compilers find the common subexpressions themselves), the division (the
  no-division probe is not faster), 512-bit vectors once the data comes from
  L3. The loop has no spills (about 40 vector operations per pass).
* From L3 all fused forms and all compilers end at 175-195 MLUPS for 8x8
  blocks: the kernel reads 9 directions of g and C and writes three fields,
  and that data movement is the cost. There is nothing left in the arithmetic.

The one-pass form is not adopted in `src/`: on the 0.5 % of the run it would
save about 0.1 %.

## What the compilers get right and wrong

* **Is gfortran the fastest?** No. On the final kernel (`mc_new`, 8x8, L2):
  256-bit ifx 95.9 > flang 87.7 > gfortran 14 84.8 > gfortran 15 69.2;
  512-bit ifx 122.4 ~ gfortran 14 119.2 > flang 110.2 > gfortran 15 104.4.
  Only for the collide-only kernel at 512 bits is gfortran 14 ahead (140.5
  against 134.0 for ifx).
* **gfortran 15.2 is 8-20 % slower than gfortran 13.3 and 14.2** on every
  form of this loop. Its code keeps about 34 array pointers (one per row of C
  and per direction of h and g) on the stack and reloads 30 of them in each
  pass (14.2: none), and it has about twice the vector spills (69 against 38
  in `mc_sym_fused_mu`). `-fno-ivopts` makes it worse (57 reloads),
  `-fira-loop-pressure -fira-region=all` and `-ffast-math` change nothing.
  (gfortran 16 could not be installed in the sandbox and was not measured.)
* All compilers vectorize the loop without the directive, fuse multiply-adds,
  turn the `MERGE` into a mask, and specialize for the constant block size.
* ifx and flang turn a collapsed or flattened loop into gathers (1.).
* flang needs `-mllvm -force-vector-width=8` for 512-bit vectors.
* None of them uses reciprocal estimates for double precision, with any flag (8.).
* Where the time goes: the kernel sustains about 1.3 of the 2 possible 512-bit
  floating-point operations per cycle (about 230 cycles for the ~290
  operations of a row), so it is not limited by the execution ports alone. The
  probes rule out the divider at 256 bits and spills as main causes; without
  hardware counters the rest could not be attributed.

## What went into the solver

`src/Collision_JCP2016.inc`, instantiated twice in `Model_JCP2016_Mod`
(`Collision_h_g_AMR` and `Moments_Collision_h_g_AMR`): the symmetric algebra
(4.), chemical potential and moments fused into the cell loop (5.), nested
loops with `!$OMP SIMD` on the inner one and explicit-shape dummies with the
block size as constants (1.-3.), exact divisions (8.). `Moments_g_Block`
(diagnostic steps, new blocks after refinement) is unchanged.

## Files

| file | contents |
|---|---|
| `bench_collision.F90` | the benchmark driver |
| `gen_variants.py`, `*.inc` | generator and fragments of the Fortran variants -> `Collision_Variants_Mod.F90`, `variant_*.inc` |
| `Model_JCP2016_Ref_Mod.f90` | the model before the study (reference results and baseline) |
| `moments_variants.inc`, `moments_kernel.inc` | the variants of the moments kernel (11.) |
| `collision_avx512.c` | the kernel with AVX-512 intrinsics |
| `build.sh`, `build_all.sh`, `run.sh` | build one / all configurations, run and tabulate |
| `asmstat.py` | instruction statistics per routine from the assembly |
| `results/study/` | the tables of the study |
