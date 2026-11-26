# Performance study of the ghost cells of h and g

`Fill_Ghost_Cells_hg` (`src/Ghost_hg_Mod.f90`) was the most expensive stage
of a time step once the collision and the streaming had been dealt with
(`bench/collision`). This directory holds the benchmark and the variants that
were tried. The outcome is in `src/Ghost_hg_Mod.f90`, `src/Ghost_Copy.inc`
and the collision loop of `src/Solver_Mod.f90`.

**Result in the solver** (rising bubble, `-DAMR_L0=512 -DAMR_MIN_LEVEL=3`, 3000
steps; seconds for collision + ghost cells of h and g, and for the run):

| build | before | after | | run before | after |
|---|---|---|---|---|---|
| gfortran 14.2, 8x8, 1 thread | 3.27 | 2.85 | 1.15 | 4.8 | 4.4 |
| gfortran 14.2, 8x8, 4 threads | 1.13 | 0.88 | 1.30 | 1.9 | 1.4 |
| gfortran 15.2, 8x8, 1 thread | 3.44 | 3.00 | 1.15 | 5.2 | 4.5 |
| ifx 2026.2, 8x8, 1 thread | 3.25 | 2.51 | 1.30 | 4.7 | 3.9 |
| gfortran 14.2, 16x16, 1 thread | 5.24 | 4.64 | 1.13 | 8.0 | 7.3 |
| PRE 2017 model, gfortran 14.2, 1500 steps | 1.71 | 1.35 | 1.26 | 2.5 | 2.2 |
| Parker 2009 model, gfortran 14.2, 1500 steps | 4.31 | 4.12 | 1.05 | 5.0 | 5.0 |

With gfortran the printed log and the diagnostics file are identical before
and after for all three models and for 1 and 4 threads; with ifx the rise
velocity differs in the last printed digit (2e-6).

## The benchmark

`bench_ghost` sets up the solver's own mesh (the rising bubble, refined level
by level as the solver does before its first step: 444 leaf blocks on levels
3-6, i.e. 176 / 32 / 92 / 144 leaves, 2792 same-level neighbour directions,
712 coarse/fine ones, 48 at walls) and times the fill of all leaf blocks.
Before every fill all blocks are collided, not timed, so that the fill finds
the caches as in a time step. A second column counts the collision as well:
it is the one that matters for the variants that move work into the
collision loop. Check: from the same state, fill, then stream every leaf
block; `err` is the largest difference of the streamed interior against the
module as it was before the study.

```bash
./build.sh gf14_8 8 gfortran-14 -O3 -fstack-arrays -march=native
./run.sh gf14_8 ifx8                       # 1 thread, pinned
OMP_NUM_THREADS=4 ./run.sh gf14_8
```

`ghost_template.inc` is the module before the study with options;
`variants.json` lists the variants (one module each, `gen_ghost.py`).

## Where the time goes

MLUPS of the fill, 1 thread, 8x8, gfortran 14 (ifx and gfortran 15 within 5 %):

| | fill | share of the fill |
|---|---|---|
| `ref`: everything | 47.8 | |
| `x_nocopy`: without the same-level copies | 108.6 | copies 56 % |
| `x_nocf`: without the coarse/fine fills | 85.2 | coarse/fine 44 % |
| `x_nocoarse`: without the interpolation from coarser blocks | 64.7 | from coarse 26 % |
| `x_nofine`: without the restriction of finer blocks | 58.6 | from fine 18 % |
| `x_nowall`: without the walls | 49.3 | walls ~0 |

A same-level copy moves about 2.5 ns per value (10 cycles): far more than a
load and a store. The values of a left or right edge are each in a cache line
of their own, in a block that the collision sweep has long pushed out of L1
and L2. The copies are bound by data movement, not by instructions.

## Experiments

Fill MLUPS | collision + fill MLUPS, 1 thread:

| 8x8 | gf14 | gf15 | ifx | flang |
|---|---|---|---|---|
| `ref` | 47.8 / 29.9 | 51.6 / 30.1 | 47.7 / 29.9 | 36.8 / 24.9 |
| `copy_k`: `SELECT CASE( k )` outside the loops, blocks as dummy arguments | 51.2 / 30.9 | 53.5 / 30.0 | 48.9 / 29.9 | 45.0 / 28.2 |
| `copy_spec`: one routine per list of directions (constants) | 53.0 / 31.9 | 54.3 / 30.5 | 53.2 / 32.2 | 48.8 / 29.6 |
| `copy_cell`: columns cell by cell, directions inside | 53.2 / 32.1 | 54.0 / 30.3 | 52.4 / 31.6 | 48.3 / 29.1 |
| `push`: a block writes its edges into its neighbours' ghost cells | 52.6 / 31.9 | 53.2 / 30.2 | 51.2 / 31.4 | 46.9 / 28.7 |
| `cfdirs`: from a coarser block only the directions read | 54.2 / 32.4 | 55.0 / 30.8 | 55.2 / 32.9 | 40.1 / 26.3 |
| `spec_cfdirs` | 59.2 / 34.0 | 61.1 / 32.8 | 58.7 / 33.6 | 54.8 / 31.8 |
| `fused_push`: the push right after the collision of the block | 105.6 / 33.8 | 112.0 / 32.4 | 98.5 / 33.2 | 88.0 / 32.0 |
| `fused_push_cfdirs` (= `src` now) | 141.1 / **37.6** | 144.0 / **34.7** | 131.0 / **37.2** | 112.7 / **34.6** |
| gain in collision + fill | 1.26 | 1.15 | 1.24 | 1.39 |

16x16 (gf14 / gf15 / ifx): `ref` 84.5 | 44.0, 85.9 | 41.3, 83.5 | 44.2;
`fused_push_cfdirs` 230.5 | 51.3, 228.2 | 47.5, 233.5 | 52.8 (1.15-1.19).
4 threads, 8x8: `ref` 155.4 | 93.1 (gf14), 158.1 | 80.8 (gf15), 155.4 | 96.5
(ifx), 111.9 | 74.6 (flang); `fused_push_cfdirs` 388.7 | 111.4, 400.5 | 104.6,
375.4 | 113.3, 332.3 | 111.2 (1.17-1.49).

* **Code shape of the copies** (the levers of the kernel studies: the
  selection outside the loops, explicit-shape dummies, compile-time
  constant directions, loop order): 5-10 % of the fill with gfortran and ifx,
  30 % with flang. The copies were already written with constant sections;
  what is left is the dispatch per direction. Column by column or cell by
  cell makes no difference.
* **Only the directions that are read** (`cfdirs`): a block whose ghost cells
  no finer block reads needs from its coarser neighbour only what its own
  streaming reads: 3 of 8 directions on a side on the finest level (1 at a
  corner), 6 (2) with Lax-Wendroff. The same-level copies already did this;
  the interpolation gathered and interpolated all 8. 13 % of the fill
  (interpolation from coarse: -35 %), results equal to round-off (the
  interpolation is instantiated for 3, 6 and 8 directions).
* **Fusion with the collision** (`fused_push`): the copies are pushed instead
  of pulled, which by itself changes nothing (`push`), but a push can run
  right after the collision of its block: the collision of a neighbour does
  not touch that neighbour's ghost cells, so the order of the blocks does not
  matter. The source edges are then in L1 and only the neighbours' ghost
  cells miss. Collision + fill: +12-13 % (flang +29 %), with `cfdirs` the
  totals above. On 4 threads the gain is larger (1.2-1.5): less traffic for
  the shared L3.
* **Not tried again**: gather / scatter of the columns (AVX-512 and AVX2, C
  kernels) and streaming fused with the fill through a second pool were
  measured before and are slower (main README).
* **Block size**: per cell the fill is 1.8 times cheaper with 16x16 blocks
  (the surface to volume ratio halves); the gain of the new form is
  1.15-1.19 there.
* **Compilers**: gfortran 14, 15 and ifx are within 5 % of one another on
  every variant: the stage is bound by memory. flang is 25 % slower on the
  original form and equal on the final one.

## What went into the solver

* `Push_Same_Level( n )` (`Ghost_hg_Mod`), called in the collision loop of
  `Solver_Mod` right after the collision of block n; the copies are the
  template `Ghost_Copy.inc`, instantiated for the three lists of directions.
* `Fill_Ghost_Cells_hg` does the rest: restrictions of refined neighbours and
  walls in the first pass, interpolation from coarser neighbours in the
  second (`Fill_From_Coarse`), for the directions read.
* The timer `t_collision` now includes the same-level ghost cells and
  `t_Ghost_hg` only the coarse/fine part and the walls.

What remains in the fill is the coarse/fine part: restriction of finer
blocks (all 8 directions are needed there, as the finer blocks interpolate
from these ghost cells) and the interpolation from coarser ones, 3 to 8 ns
per leaf cell on this mesh.
