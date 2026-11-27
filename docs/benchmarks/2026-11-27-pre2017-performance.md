# PRE 2017 model: performance after the kernel study, 2026-11-27

The findings of the kernel study (`2026-11-26-kernel-performance.md`) were
applied to the PRE 2017 model. Outcome:

* The changes that are shared by all models (streaming, ghost cells of h and
  g, ghost layers of C) make the PRE 2017 rising bubble 1.21-1.41 times
  faster than before the study.
* The two levers of the JCP 2016 collision, one fused pass and stress moments
  from pair sums, do **not** carry over to the PRE 2017 collision with
  gfortran. They do with ifx. The model is left as it is.

## Whole runs

![Time per cell update of the PRE 2017 driver](2026-11-27-pre2017-whole-run.svg)

Rising bubble with `amr_lbm_pre2017`, `-DAMR_L0=512 -DAMR_MIN_LEVEL=3`, 3000
steps, 8x8 blocks, MLUPS for the whole time step, best of 3. `original` is
commit 209db98, `master` the merged kernel study (7f1bdeb), `two loops` the
candidate described below, which was not adopted.

| build | threads | original | master | | two loops | against master |
|---|---|---|---|---|---|---|
| gfortran 14.2 | 1 | 19.5 | 23.6 | 1.21 | 23.3 | -1 % |
| | 4 | 52.2 | 63.1 | 1.21 | 67.0 | +6 % |
| gfortran 15.2 | 1 | 17.2 | 24.1 | 1.40 | 22.8 | -5 % |
| | 4 | 46.0 | 64.7 | 1.41 | 60.3 | -7 % |
| ifx 2026.2 | 1 | 18.1 | 22.3 | 1.23 | 23.5 | +6 % |
| | 4 | 46.8 | 63.2 | 1.35 | 67.1 | +6 % |

The diagnostics of the three states agree to the last printed digit (6e-7);
the mesh is the same.

## The collision kernel

On master a block is processed in four passes: chemical potential, gradient
of C, moments of g, collision. The collision already pairs opposite
directions. Three rearrangements were tried, all with results equal to
round-off:

* **one loop**: gradient, chemical potential, moments and collision in one
  pass over the cells, the stress and the stress moments from the sums of g
  over the pairs (no equilibrium per direction), every g stored as soon as it
  is computed;
* **two loops**: a moments loop and a collision loop, each with the gradient
  and the chemical potential computed in the loop;
* **hybrid**: the three moments passes of master and the new collision loop.

![The collision kernel in the arrangements tried](2026-11-27-pre2017-kernel-variants.svg)

Kernel only, on the solver's mesh (492 leaf blocks, 8.4 MiB, in L3), 1 thread,
MLUPS, best of 6 interleaved runs:

| moments + collision | gfortran 14.2 | gfortran 15.2 | ifx 2026.2 |
|---|---|---|---|
| master: four passes | 64.4 | 62.7 | 64.2 |
| one loop | 63.5 (-1 %) | 59.1 (-6 %) | 75.7 (+18 %) |
| two loops | 64.8 (+1 %) | 62.7 (0 %) | 70.4 (+10 %) |

| collision only | gfortran 14.2 | gfortran 15.2 | ifx 2026.2 |
|---|---|---|---|
| master | 79.8 | 75.8 | 78.0 |
| one loop | 88.7 (+11 %) | 82.0 (+8 %) | 94.5 (+21 %) |
| two loops | 87.0 (+9 %) | 85.2 (+12 %) | 93.8 (+20 %) |

The hybrid, measured in a run of its own, was 2-6 % slower than master on the
path with the moments for all three compilers.

* The collision alone gets 8-21 % faster with every compiler, but that path
  only runs on the steps with diagnostics.
* On the path of a normal step gfortran gains nothing from any arrangement,
  and gfortran 15 loses 6 % with one loop. The same holds on a small mesh
  that fits in L2 (gfortran 14: 73 MLUPS on master, 71 with one or two
  loops), so memory is not what hides a gain.
* Computing the gradient and the chemical potential once into two small
  arrays and reading them in two loops, as master does, is as cheap as
  computing them in each loop: the hybrid and the two-loop form do not beat
  it.
* Neutral, as for JCP 2016: the `!$OMP SIMD` directive on the loops (within
  2 % either way) and merging the divisions of the relaxation time and the
  density into one (within noise).
* In the solver the two-loop form is 5-7 % slower than master with gfortran
  15, where the kernel benchmark showed no difference.

Why the model is left unchanged: no arrangement is at least as fast as master
with all three compilers on the path that matters. One fused loop would be
worth about 6 % of a run with ifx; keeping it next to the four passes for
gfortran would mean two implementations of the collision.

## `DO CONCURRENT` and the pair-sum collision

Two further ideas, measured with the kernel benchmark and not adopted.
MLUPS of moments + collision, 1 thread, best of 5-6 interleaved runs; the
numbers of one compiler were taken together, those of different tables on
different occasions.

### `DO CONCURRENT` in place of `DO j` / `DO i`

Applied to the collision loops of both models, without the `!$OMP SIMD`
directive and without a locality clause.

| compiler      | model    | `DO j` / `DO i` (master) | `DO j` / `DO CONCURRENT (i)` | `DO CONCURRENT (j, i)` |
|---------------|----------|-----:|-----:|-----:|
| gfortran 14.2 | PRE 2017 | 66.9 | 66.7 | 67.6 |
| gfortran 14.2 | JCP 2016 | 67.8 | 68.6 | 68.5 |
| gfortran 15.2 | PRE 2017 | 64.8 | 63.8 | 64.7 |
| gfortran 15.2 | JCP 2016 | 57.0 | 56.1 | 53.7 |
| ifx 2026.2    | PRE 2017 | 64.3 | 14.4 | 26.5 |
| ifx 2026.2    | JCP 2016 | 68.8 | 20.9 | 31.2 |

* gfortran: no difference beyond the timing noise, except the two-index
  form of JCP 2016 with 15.2 (-6 %). The plain loops already vectorize, so
  `DO CONCURRENT` has nothing to add.
* ifx with `-qopenmp`: 2 to 4.5 times slower. Under that flag ifx turns
  every `DO CONCURRENT` into an OpenMP work-sharing loop; its optimization
  report lists the loops as `OMP LOOP` and the enclosing loop as not
  vectorized because of a function call. Here that loop sits inside the
  parallel loop over blocks.
* A `LOCAL( ... )` clause with all scalars of the loop body does not change
  this (ifx: PRE 2017 17.4 against 63.2, JCP 2016 24.8 against 67.3;
  gfortran 15.2 unchanged). gfortran 14 rejects the clause.
* ifx has no flag that maps `DO CONCURRENT` on the host to SIMD only; its
  options for the construct concern offload
  (`-fopenmp-target-do-concurrent`). What works is to compile the two model
  files with `-qopenmp-simd` instead of `-qopenmp`, which is possible
  because they hold no OpenMP construct other than `SIMD`:

| ifx, model files with `-qopenmp-simd` | master | `DO CONCURRENT (i)` | `DO CONCURRENT (j, i)` |
|----------|-----:|-----:|-----:|
| PRE 2017 | 63.6 | 65.1 | 66.0 |
| JCP 2016 | 68.0 | 68.4 | 69.0 |

  With that, `DO CONCURRENT` is as fast as the plain loops (+1 to +4 %,
  inside the noise), at the price of per-file compiler flags.

The plain loops stay: `DO CONCURRENT` gains nothing with any compiler and
is a trap with ifx under the flags of this build.

### Pair sums in the PRE 2017 collision alone

The collision of master already treats opposite directions in pairs. The
further step of the JCP 2016 kernel is to take the stress and the two
relaxed moments from sums over pairs, without forming the equilibrium of
each direction. Above it was measured only together with fusing the
gradient and the chemical potential into the loop; here it is applied
alone, keeping the four passes of master.

| compiler      | moments + collision: master | pair sums | collision only: master | pair sums |
|---------------|-----:|-----:|-----:|-----:|
| gfortran 14.2 | 66.9 | 68.7 | 82.9 | 87.0 |
| gfortran 15.2 | 64.8 | 62.5 | 79.4 | 78.7 |
| ifx 2026.2    | 64.3 | 65.0 | 78.2 | 77.0 |

Between -3 % and +3 % on the path the solver takes at every step, +5 % to
-2 % on the collision-only path: within the noise. The algebra is not what
limits this kernel.

## What limits `Collision_h_g`

The kernel alone in a standalone benchmark (gfortran 14.2, 1 thread, 8x8
blocks, MLUPS, range of two runs), with the pool of blocks sized to fit in
each cache level. "No div" replaces the four divisions and the square root
by multiplications: wrong results, timing only.

| pool of blocks  | 256-bit | 256-bit, no div | 512-bit | 512-bit, no div |
|-----------------|--------:|--------:|--------:|--------:|
| 18 KiB (L1)     | 142     | 146-149 | 171-187 | 240-245 |
| 1.2 MiB (L2)    | 127-131 | 125-130 | 157-161 | 174-178 |
| 37 MiB (L3)     | 87-91   | 82-85   | 88-90   | 90-92   |

* In the solver the pool is in L3, and there the kernel is bound by memory
  traffic: every variant gives about 90 MLUPS, which is also what the
  solver's collision-only path reaches (about 83 with the gradient and the
  chemical potential). A block of 8x8 cells occupies 18 KiB, 284 bytes per
  cell with the ghost layers, half of it written back.
* In cache it is bound by the throughput of the floating-point ports, not
  by a stall. The inner loop has 222 instructions per 4 cells, 157 of them
  vector arithmetic (69 FMA), 45 memory operands and 20 stores, and runs at
  about 3 instructions per cycle.
* The divisions do not limit it at 256 bits (3 % in L1, nothing beyond).
  They do at 512 bits in L1 only.
* This is why the rearrangements above change nothing: they shorten
  arithmetic the solver is not waiting for.

Block size, 256 bits, pool in L3:

| block | bytes per cell | MLUPS   |
|-------|---------------:|--------:|
| 8x8   | 284            | 87-91   |
| 16x16 | 237            | 105     |

Larger blocks load fewer ghost cells per interior cell and give the
prefetcher longer runs. This is the kernel only; a run with larger blocks
also refines more cells around the interface.

## Setup

Machine and compilers as in `2026-11-26-kernel-performance.md` (4 cores of a
virtual machine, timing noise of a few percent; no hardware counters). The
kernel benchmark sets up the solver's mesh for the rising bubble, restores
the same state before every sweep, and times
`model%Moments_and_Collide` and `model%Collide` over all leaf blocks.
