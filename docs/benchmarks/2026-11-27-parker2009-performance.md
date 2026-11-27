# Parker 2009 model: collision kernel, 2026-11-27

The collision of the two-species model was rewritten with the methods of
the kernel study (`2026-11-26-kernel-performance.md`). A spinodal run gets
2.3 to 3.7 times faster; the collision kernel alone 5.6 to 7.9 times.
Results agree with the previous code to round-off.

All throughput numbers are MLUPS: million leaf-cell updates per second of
wall time. *Before* is commit 654e836, *after* is this change.

## Setup

* Machine: Intel Xeon (Emerald Rapids, AVX-512), 4 cores of a virtual
  machine, 2 MiB L2 per core, no hardware counters. The timing noise was up
  to 10 % in these sessions; numbers are the best of the stated runs.
* Compilers: gfortran 14.2, gfortran 15.2, ifx 2026.2, flang 23, with the
  flags of the CMake build (flang: `-O3 -march=native` in the kernel
  benchmark only).
* Case: `cases/parker2009_spinodal.nml`, `-DAMR_L0=256 -DAMR_MIN_LEVEL=2`,
  8x8 blocks, 1000 steps, about 75 000 leaf cells (1200 leaf blocks) on
  levels 2 to 5, output off.
* Kernel benchmark: the solver's mesh and initial state of that case; times
  the *moments + collision* path (`model%Moments_and_Collide`, every time
  step) and the *collision only* path (`model%Collide`, steps with
  diagnostics) over all leaf blocks, 1 thread.

## Whole runs

Best of 3 interleaved runs. MLUPS for the whole time step.

| compiler      | threads | before | after | speed-up |
|---------------|--------:|-------:|------:|---------:|
| gfortran 14.2 | 1 | 5.9  | 21.6 | 3.7 |
| gfortran 14.2 | 4 | 22.5 | 58.6 | 2.6 |
| gfortran 15.2 | 1 | 6.3  | 22.6 | 3.6 |
| gfortran 15.2 | 4 | 21.2 | 64.6 | 3.0 |
| ifx 2026.2    | 1 | 8.5  | 22.7 | 2.7 |
| ifx 2026.2    | 4 | 28.3 | 66.2 | 2.3 |

Seconds per stage, 1 thread:

| stage | gfortran 14.2 | gfortran 15.2 | ifx 2026.2 |
|-------|--------------:|--------------:|-----------:|
| collision + same-level ghost cells of h, g | 10.15 -> 1.98 | 9.64 -> 1.88 | 6.90 -> 1.73 |
| streaming and C                            | 1.30 -> 0.87  | 1.13 -> 0.83 | 0.88 -> 0.84 |
| ghost layers of C                          | 0.26 -> 0.18  | 0.22 -> 0.17 | 0.20 -> 0.20 |
| coarse/fine ghost cells of h, g, walls     | 0.23 -> 0.14  | 0.18 -> 0.13 | 0.15 -> 0.14 |
| mesh adaptation                            | 0.14 -> 0.13  | 0.14 -> 0.13 | 0.23 -> 0.23 |
| total                                      | 12.07 -> 3.30 | 11.30 -> 3.14 | 8.36 -> 3.14 |

The collision was 83 to 85 % of a run and is now 55 to 60 %. Only the
collision changed; the other stages differ between the two columns by the
state of the machine (gfortran) or not at all (ifx).

## The kernel

Kernel benchmark, *moments + collision*, best of 3 runs:

| compiler      | before | after | speed-up |
|---------------|-------:|------:|---------:|
| gfortran 14.2 | 9.0  | 58.0 | 6.4 |
| gfortran 15.2 | 8.3  | 54.0 | 6.5 |
| ifx 2026.2    | 11.5 | 64.4 | 5.6 |
| flang 23      | 4.2  | 33.0 | 7.9 |

The populations after a collision differ from those of the previous code
by at most 4e-16 (values up to 0.4) with every compiler.

What changed:

* **A loop the compiler can vectorize.** Before, `Collision` called
  `Collide_Species` once per cell and species with the array section
  `h(i,j,:)`, and each call ran a loop over the nine directions. Now the
  body of one species is a Fypp macro, expanded for h and for g inside one
  loop over the cells of a block (function inlining by the template, not
  left to the compiler).
* **Raw moments from pairs.** Opposite directions enter the moments only
  as their sum (even moments) or their difference (odd moments): 16
  additions per species instead of nine directions times eight products.
* **Central moments by two one-dimensional shifts**, along x and then
  along y, each of the form `k1 = m1 - u m0`, `k2 = m2 - u (m1 + k1)`: 12
  multiplications instead of the expanded binomial formulas.
* **Only the changes are transformed back.** The transformations are
  linear, so the relaxation changes of the central Hermite moments are
  shifted back and added to the populations. The density is conserved, so
  its change is zero and all its terms drop out.
* **One gradient for both species.** With
  `G = -D grad( ln gamma_1 - ln gamma_2 )` the diffusion velocities relative
  to the mass average are `u'_1 = x2 G` and `u'_2 = -x1 G`. One gradient of
  one field replaces two gradients, and `c1 u'_1 = -c2 u'_2`.
* **No division in the loop.** `1/(c1 + c2)` is stored by the pass that
  evaluates the activity coefficients on the block and its ghost layer.
* **Moments fused with the collision.** The first moments of the two
  species, summed, are the mixture momentum, so the path of a normal time
  step gets the velocity from sums it needs anyway. `Moments_Collision` and
  `Collision` are generated from one template.

## Variants measured

Kernel benchmark, *moments + collision*, change against the adopted form,
best of 3 or 4 runs:

| variant | gfortran 14.2 | gfortran 15.2 | ifx 2026.2 | flang 23 |
|---------|--------------:|--------------:|-----------:|---------:|
| species locals declared in a `BLOCK` inside the macro (adopted) | -2 % | +1 % | -2 % | -3 % |
| full back-transform: the relaxed moments are transformed and the populations assigned | +3 % | +2 % | -6 % | -5 % |
| one loop per species | -2 % | within noise | -13 % | -13 % |
| activity terms in the loop: no pass before it, `ln gamma` of the eight neighbours evaluated per cell | -18 % | -19 % | -21 % | -24 % |
| species body as a module procedure, called per cell | -62 % | -61 % | -69 % | -52 % |
| species body as an internal procedure of the kernel | -64 % | -64 % | -67 % | -52 % |

* The `BLOCK` only moves the declarations of the per-species temporaries
  next to their use; the timings are within the noise of the form without.
* The full back-transform is within the noise of the adopted form. A
  second session of 6 runs gave -1 % (gfortran 14.2), +2 % (gfortran 15.2)
  and -5 % (ifx).
* Splitting the species into two loops is not faster with any compiler.
* The pass before the loop evaluates each activity term and each division
  once per cell (100 for a block of 64 cells with its ghost layer). In the
  loop, every cell evaluates its eight neighbours: 9 divisions per cell.
* A procedure call per cell and species, module or internal, costs a factor
  of 2.1 to 3.2 with every compiler: the loop over the cells is then not
  vectorized. Hence the macro.

## Not done

* Timings on Arm (the model is also run on Apple Silicon) and with 16x16
  blocks.
* The remaining stages. After this change streaming is about a quarter of
  a run; it handles two sets of populations and two scalars, as before.
