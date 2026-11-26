# Kernel performance study, 2026-11-26

The kernels of a time step were studied one by one and rewritten where a
measured gain justified it. Whole runs became 1.33-1.48 times faster on one
thread and 1.22-1.32 times on four; results are equal to round-off.

## Setup

* Case: rising bubble, JCP 2016 model, `-DAMR_L0=512 -DAMR_MIN_LEVEL=3`, 3000
  steps (about 470 leaf blocks on levels 3-6, 8 MiB of block data).
* Before: commit 209db98. After: this change.
* Machine: Intel Xeon (Emerald Rapids, AVX-512), 4 cores of a virtual
  machine, 2 MiB L2 per core. Timings in a VM are noisy: every number is the
  best of 3 runs and differences below about 4 % are not significant. No
  hardware counters were available.
* Compilers: gfortran 14.2 and 15.2, ifx 2026.2, with the flags of the CMake
  build (`-O3 -march=native` / `-O3 -xHost`; 256-bit vectors).
* MLUPS: million leaf-cell updates per second of wall time for the whole time
  step; every leaf cell counts once per step whatever its level. It is used
  here only to compare the same simulation before and after.

## Whole run

| build | threads | before (MLUPS) | after | speed-up |
|---|---|---|---|---|
| gfortran 14.2, 8x8 blocks | 1 | 17.2 | 22.8 | 1.33 |
| | 4 | 46.1 | 57.5 | 1.25 |
| gfortran 15.2, 8x8 | 1 | 14.9 | 22.1 | 1.48 |
| | 4 | 42.6 | 53.6 | 1.26 |
| ifx 2026.2, 8x8 | 1 | 16.5 | 23.0 | 1.39 |
| | 4 | 43.2 | 57.1 | 1.32 |
| gfortran 14.2, 16x16 | 1 | 21.5 | 29.9 | 1.39 |
| | 4 | 67.7 | 82.6 | 1.22 |
| ifx 2026.2, 16x16 | 1 | 22.1 | 30.5 | 1.38 |
| | 4 | 68.3 | 83.7 | 1.23 |

A full run to t* = 8 (36000 steps, gfortran 14.2, 1 thread) took 88.6 s before
and 64.6 s after. Over that run the leaves on the finest level grow from 192
to 328 and the speed in MLUPS falls by about 5 %; the speed-up is 1.37 early
and 1.29 late.

## By stage

Seconds for 3000 steps, 1 thread, 8x8 blocks. Collision and ghost cells of h
and g are given together because part of the ghost work moved into the
collision loop.

| stage | gfortran 14.2 | gfortran 15.2 | ifx 2026.2 |
|---|---|---|---|
| collision + ghost cells of h, g | 3.74 -> 2.73 | 3.75 -> 2.83 | 3.76 -> 2.69 |
| streaming and C | 0.97 -> 0.89 | 1.72 -> 0.90 | 1.04 -> 0.90 |
| ghost layers of C | 0.36 -> 0.17 | 0.37 -> 0.16 | 0.51 -> 0.19 |
| mesh adaptation, output | 0.14 -> 0.13 | 0.15 -> 0.15 | 0.11 -> 0.12 |
| total | 5.21 -> 3.92 | 6.00 -> 4.04 | 5.41 -> 3.89 |

## What changed and why

**Collision of the JCP 2016 model** (`Model_JCP2016_Mod.fypp`, generated in
two forms). The kernel alone became 1.3-1.5 times faster.

* The chemical potential, the moments of g and the collision are computed in
  one pass over the cells instead of three (plus a gradient pass). The
  differences of C are shared, g is loaded once and no work array is passed
  between loops. Fusion was worth 25-50 % on the path of a normal step;
  splitting the loop (an h pass and a g pass, or the divisions first) was
  7-25 % slower.
* All moments of g but the two stresses relax with rate 1, and only the sums
  of the equilibrium over a pair of opposite directions enter the stress
  moments. These sums have closed forms, so the individual `g - geq` are
  never formed: 21-30 % faster, about a quarter fewer operations.
* It is generated in two forms because the solver needs it with the moments computed in
  the loop (normal steps) and with P, U given (first step, steps with
  diagnostics). A single routine with the choice inside the cell loop, or the
  cell body as an internal procedure, is not vectorized by gfortran and three
  times slower.
* Tried and neutral: loop flattening (the row length is a multiple of the
  vector length; `COLLAPSE(2)` or a flat index makes ifx and flang emit
  gathers, at half the speed), default lower bounds and sequence association
  (identical code), the `!$OMP SIMD` directive (all compilers vectorize the
  loops anyway, so the collision and the moments kernel carry none). Run-time block sizes cost 6-20 %.
* The three divisions and the square root per cell cost 2-3 % with 256-bit
  vectors. They are needed in every cell and stay exact.

**Streaming** (`LBM_Common_Mod`). Both kernels compute a direction into a
temporary of one plane and copy it back before the next direction, instead of
one temporary for all directions and a single array assignment. Results are
bitwise identical. gfortran 15 had compiled the single copy into several
hundred shuffles, 2.4 times slower than gfortran 14. Unrolling the loop over
the directions by directive, linear indexing and flat loops gained nothing.

**9-point streaming as an option** (environment variable
`AMR_COARSE_STENCIL=9` at run time, default 3;
`Streaming_LaxWendroff9Pt_AMR`, selected through a procedure pointer). On
the levels below the finest the diagonal directions can use the tensor
product of the 1-D Lax-Wendroff schemes in x and y (all nine neighbours)
instead of three points along the direction; the axis directions are the
same. It is the 2-D Lax-Wendroff scheme with its 4-point cross derivative
plus terms of third and fourth order.

* Cost: done as a pass in x and a pass in y it costs the same as the 3-point
  scheme (streaming stage 0.90 s against 0.93 s for 3000 steps; a full run
  to t* = 8 on 4 threads 23.1 s against 24.2 s). A general 9-point stencil
  with nine weights per direction would cost 1.5 times as much.
* It reads all four diagonal directions in the ghost corners, so the corner
  lists of the ghost cells grow from 2 to 4 directions, and the wall values
  of the end columns are redone after a side has been interpolated from a
  coarser block.
* Behaviour on the rising bubble to t* = 8 (36000 steps): stable, total mass
  changes by -1.3e-4 (3 points: -1.7e-4), bubble centroid within 0.2 % of the
  3-point run. The late flow differs: at t* = 8 the rise velocity is 0.0177
  against 0.0141 and the largest velocity 0.060 against 0.029. Which of the
  two is closer to the uniform-grid solution was not examined.
* Known approximation: at a wall, the corner ghost cell of a block takes the
  bounce-back value from the end column instead of the neighbour's column
  (the clamp of `Bounce_Back_3/4`). The 3-point scheme never reads that
  value; the 9-point scheme reads it with a weight of about 0.05.

**Ghost cells of h and g** (`Ghost_hg_Mod`, the collision loop of
`Solver_Mod`). The same-level copies were 56 % of this stage and are
bound by cache misses: each value of a left or right edge is in its own cache
line of a block long evicted.

* Every block now writes its edges into the ghost cells of its same-level
  neighbours right after its own collision (`Push_Same_Level`), while it is
  in L1. This is valid in any block order, as the collision of a neighbour
  does not touch its ghost cells.
* From a coarser neighbour only the directions the block's streaming reads
  are interpolated (3 or 6 of 8 on a side), unless finer blocks read these
  ghost cells.
* The timer `t_collision` now includes the same-level copies and `t_Ghost_hg`
  only the coarse/fine part and the walls.

**Ghost layers of C** (`Ghost_C_Mod`). The strip routines are generated for
one and for two layers, with the interpolation
fixed to each strip shape, instead of one set of routines with the sizes at
run time: 2-2.5 times faster. Pushing the copies from the streaming loop was
tried and gains nothing more on 8x8 blocks. The run-time-size strip
interpolations of `Interpolations_Mod` had no caller left and were removed.

**Fypp instead of the C preprocessor.** The routines that exist in several
forms (the collision with and without the moments, the copies per list of
directions, the strips per number of ghost layers, the interpolations per
strip shape) are generated by Fypp loops and macros in the `.fypp` sources,
which CMake turns into plain Fortran in the build tree. The configuration
(every CMake option, under its own name) is passed to Fypp as well; no source
goes through the C preprocessor. The generated code is the same as with the
C-preprocessor templates it replaces: identical diagnostics, same timings.

## Validation

* The mesh (blocks and leaves at every output) is the same before and after.
* The diagnostics (`bubble_rise.dat`, density range, maximum velocity, mass)
  agree to the last printed digits (1e-6 to 7e-6) after 3000 steps: JCP 2016
  with gfortran 14.2 (1 and 4 threads, 8x8 and 16x16 blocks), gfortran 15.2
  and ifx; PRE 2017 and Parker 2009 (1500 steps) with gfortran 14.2. After
  20000 steps of an intermediate state of this change the bubble centroid
  agreed to 4e-6 and the rise velocity to 2e-4 (growth of round-off).
* A single collision agrees with the previous kernel to 1e-14, a fill of the
  ghost layers of C to 4e-16; the streaming and the same-level copies are
  bitwise identical.
* A Debug build (`-fcheck=all`) runs without a run-time error.

## Tried and not adopted

* AVX-512 intrinsics for the collision (a C kernel): no faster than the
  Fortran kernel compiled for 512-bit vectors; only reciprocal estimates with
  Newton steps gained 3-10 %.
* A layout with padded, cache-line aligned rows: 2-5 % in L2, slower from L3.
* `TASKLOOP` instead of the static loops over the levels: 14-29 % slower at
  the best grain size.
* A cap on the number of threads for small meshes: no case found on 4 cores
  where more threads were slower for the whole step.

## Compilers

* On the final collision kernel ifx is ahead at 256 bits (96 MLUPS in L2
  against 88 for flang, 85 for gfortran 14 and 69 for gfortran 15).
* gfortran 15.2 is 8-20 % slower than 13.3 and 14.2 on the collision loop
  (about 30 array pointers reloaded from the stack per pass, twice the vector
  spills). gfortran 16 was not measured.
