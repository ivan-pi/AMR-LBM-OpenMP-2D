# AMR-LBM-OpenMP-2D
Adaptive Mesh Refinement (AMR) + Lattice-Boltzmann Method (LBM) + OpenMP

Based on the following papers:

A. Fakhari, D. Bolster, L.-S. Luo, “A weighted multiple-relaxation-time lattice Boltzmann method for multiphase flows and its application in partial coalescence cascades”, Journal of Computational Physics 341:22–43 (2017) [https://doi.org/10.1016/j.jcp.2017.03.062].

A. Fakhari, M. Geier, T. Lee, “A mass-conserving lattice Boltzmann method with dynamic grid refinement for immiscible two-phase flows”, Journal of Computational Physics 315:434–457 (2016) [https://doi.org/10.1016/j.jcp.2016.03.058].

A. Fakhari, T. Mitchell, C. Leonardi, D. Bolster, “Improved locality of the phase-field lattice-Boltzmann model for immiscible fluids at high density ratios”, Physical Review E 96:053301 (2017) [https://doi.org/10.1103/PhysRevE.96.053301] (second model; uniform-grid reference code: https://github.com/AbbasHub/Improved_Phase-Field_LBM_2D).

This package is for numerical simulation of multiphase flows (rising bubble).
It employs OpenMP parallelization for share-memory machines.

The computational setup is a 2D bubble rising due to buoyancy force.
The computational domain is periodic in the x-direction and no-slip (bounce-back) at the bottom and top boundaries.

The geometry, boundary conditions, density ratio 1000, viscosity ratio 100,
body force and refinement criterion (|∇C| with ε_r = 0.002, ε_d = 0.001) are
those of the rising-bubble case of Fakhari, Geier & Lee (2016), Sec. 5.1.

### Case parameters

The dimensionless numbers are read at start-up from the namelist
`&rising_bubble` in the file named by the environment variable `AMR_CASE`:

| variable | default | meaning |
|---|---|---|
| `Eo` | 40 | Eötvös number Gy (ρh−ρl) D² / σ |
| `Re` | 40 | gravity Reynolds number √(Gy ρh (ρh−ρl) D³) / μh |
| `rho_ratio`, `mu_ratio` | 1000, 100 | ρh/ρl, μh/μl |
| `Gy_lu` | 5e-6 | gravity in lattice units; sets the Mach number and the compressibility error (the bubble swells as it rises, in proportion to `Gy_lu`) |
| `Pe` / `Mobility` | – / 0.03 | Péclet number √(Gy D) W / M if `Pe > 0`, else the mobility M |
| `Cn` / `Width` | – / 5 | Cahn number W / L0 if `Cn > 0`, else the interface width W |
| `t_end` | 8 | end time t* = t √(Gy/D) |

Without `AMR_CASE` the defaults give the original case (Eo = 40, Re_h = 40,
Pe ≈ 2.7, Cn ≈ 0.02 with L0 = 256). `cases/paper2016_Eo{5,20,200}.nml` are the
runs of Figs. 3–5 of the paper (Re_h = 100, Pe = 1.1, Cn = 0.010, t* = 12);
they need L0 = 512 (l_max = 6):

```
cmake -B build-paper -DAMR_L0=512 && cmake --build build-paper -j
AMR_CASE=cases/paper2016_Eo20.nml ./build-paper/amr_lbm
```

Every run writes `bubble_rise.dat` (120 samples and every output): the
centroid height y_c and rise velocity u_r of the light phase (cells with
C < 0.9, weight 1 − C) and the rise Reynolds number Re_r = ρh u_r D / μh.
The light-phase average u_r is noisy (pressure waves in the light phase);
d y_c/dt is the better rise velocity. The paper cases use `Gy_lu = 2.5e-6`,
which reproduces the bubbles of Fig. 4:

```
python3 tools/plot_rise.py jcp/bubble_rise.dat:JCP2016 pre/bubble_rise.dat:PRE2017 -o rise.png
```

## Building and running

The build uses CMake (3.20 or later) and needs a Fortran compiler with OpenMP
(tested with gfortran 13 and Intel ifx 2026.1):

```
cmake -B build                                  # Release, double precision
cmake --build build -j
OMP_NUM_THREADS=8 ./build/amr_lbm [number_of_time_steps [output_interval_in_steps]]
OMP_NUM_THREADS=8 ./build/amr_lbm_pre2017 [...]    # same case, PRE 2017 model
```

Two drivers are built, one per LBM model (see *Models*): `amr_lbm` (JCP 2016)
and `amr_lbm_pre2017` (PRE 2017).

The block fields live in one pool array per field, allocated at start-up for
the blocks of the fully refined quadtree (`4**(maxLevels+1) - 1`, e.g. 4095
blocks or 67 MiB for the default case). Only the slots that get used are
written, so the operating system commits only those pages. The environment
variable `AMR_BLOCK_CAPACITY` sets a smaller start capacity; the pool then
doubles whenever a refinement needs more blocks.

Configuration options (`cmake -B build -D<option>=<value>`):

| option | default | meaning |
|---|---|---|
| `CMAKE_Fortran_COMPILER` | system default | e.g. `ifx` |
| `CMAKE_BUILD_TYPE` | `Release` | `Debug` adds bounds checking and warnings |
| `AMR_PRECISION` | `8` | working precision `wp`: `8` (real64) or `4` (real32) |
| `AMR_L0` | `256` | base block width in finest lattice units (domain `L0 x 3*L0`) |
| `AMR_MIN_LEVEL` | `1` | minimum (uniform background) refinement level |
| `AMR_BLOCK_SIZE` | `8` | cells per block in each direction (`16` is faster for larger cases) |
| `AMR_WITH_FVTKHDF` | `OFF` | VTKHDF output; set `CMAKE_PREFIX_PATH` to the fVTKHDF install |
| `AMR_BUILD_OFFLOAD` | `OFF` | build the GPU offload prototype in `offload/` (ifx only) |
| `AMR_NATIVE` | `ON` | optimize for the build machine (`-march=native` / `-xHost`) |

The number of threads is taken from `OMP_NUM_THREADS`. The optional arguments
override the number of time steps (default: 8 bubble-rise time units) and the
output interval in time steps (default: one twentieth of the run; use a smaller
interval for smooth animations). The initial and the final state are always
written.

All floating-point data use the kind parameter `wp` from `Kinds_Mod`
(`real64` by default, `real32` with `-DAMR_PRECISION=4`). The total-mass diagnostic is
always accumulated in double precision.

## Source layout

| file (`src/`) | contents |
|---|---|
| `main_jcp2016.f90`, `main_pre2017.f90` | driver programs: select the model and the case |
| `Solver_Mod.f90` | the solver: initialization and time loop for any model and case, progress and timers |
| `Case_Mod.f90`, `Rising_Bubble_Mod.f90` | a test case as a set of procedures (setup, initial state, diagnostics, finish); the rising bubble |
| `Model_Mod.f90` | an LBM model as a set of block kernels (initialize, collide, moments of g) and its refinement criterion |
| `Model_JCP2016_Mod.f90`, `Model_PRE2017_Mod.f90` | the two models' kernels and criteria |
| `LBM_Common_Mod.f90`, `Phase_Field_Mod.f90` | kernels shared by the models: streaming and C = sum(h); finite differences of C (chemical potential, gradient, largest gradient) |
| `Kinds_Mod.f90`, `Config_Mod.f90`, `D2Q9_Mod.f90`, `Physics_Mod.f90` | kinds, compile-time block size and domain, lattice, physical parameters (set by the case) |
| `Blocks_Mod.f90`, `Mesh_Mod.f90` | block storage: integer metadata `Blk(:)` and one pool array per field (`Blk_h(:,:,:,n)`, `Blk_C(:,:,n)`, ...), leaf lists, free-slot queue; the domain and child blocks, the neighbour table (rebuilt from the block positions after every mesh change) |
| `Refine_Mod.f90`, `Prolong_Restrict_Mod.f90` | flagging, refinement/derefinement; data transfer to children and parents |
| `Ghost_C_Mod.f90`, `Ghost_hg_Mod.f90`, `BC_Mod.f90`, `Interpolations_Mod.f90` | ghost cells of C (1 or 2 layers) and of h/g, wall conditions, coarse/fine interpolation |
| `Output_Mod.f90`, `VTK_Output_Mod.f90`, `TecPlot_Mod.f90` | writer selection (`AMR_OUTPUT`), VTK/VTKHDF and TecPlot writers |

`code_structure.mmd` is a Mermaid diagram of the module dependencies.

Every file holds one module of the same name, so all calls are checked against
explicit interfaces at compile time. Modules are `PRIVATE` by default and list
their public entities explicitly; helpers used by a single routine are internal
procedures of that routine.

The solver knows the model and the case only through the procedure pointers
of `LBM_Model` (`Model_Mod`) and `AMR_Case` (`Case_Mod`): a new case is a new
module with a constructor like `Rising_Bubble_Case()`, a new model one like
`JCP2016_Model()`. Each model brings its refinement criterion (`Block_Error`
interface) and tolerances, since the blocks it needs on the finest level depend
on the model; mesh adaptation takes the model (`Refine_Derefine_Mesh(model)`)
and also uses its moments kernel for new blocks. A case returns its end time
and its time unit per lattice time step from `Setup`; the solver passes times in
that unit to the output and to the case's diagnostics. The domain shape
(`nBaseX` x `nBaseY` base blocks of `L0` x `L0`) is set in `Config_Mod`.
All sources go through the preprocessor
(CMake `Fortran_PREPROCESS`); the case size, precision and optional fVTKHDF
support are set with the CMake options above.

## Models

Both models solve the conservative Allen–Cahn equation for the phase field C
(distribution h) and a pressure-based hydrodynamic equation (distribution g)
on the same AMR mesh, ghost-cell, streaming and output machinery.

| | `amr_lbm` – JCP 2016 | `amr_lbm_pre2017` – PRE 2017 |
|---|---|---|
| interface sharpening | in the equilibrium of h | source term of the h equation |
| hydrodynamic forcing | central and mixed directional differences of C (two ghost layers of C) | surface tension, pressure and viscous forces from the isotropic gradient and Laplacian (one ghost layer) |
| g | physical pressure | normalised by the density (sum g = 3p/ρ) |
| collision of g | reduced MRT: stress moments relax with 1/τ(C), all others with 1 | the same MRT as in the paper, Eqs. (27)–(29) (the reference code uses BGK); viscous force from the stress relaxed with the same rates, Eq. (31) |
| relaxation time τ(C) | 1/τ interpolated linearly in C (harmonic) | dynamic viscosity interpolated linearly in C, τ = μ/(ρc²ₛ), Eqs. (24)–(25) |
| refinement criterion (in the model module) | `Interface_Error`: maximum gradient of C in a block | `Interface_Light_Phase_Error`: the interface, and additionally the light phase (max(1−C) of a block) |

The PRE 2017 model needs the bubble interior on the finest level: its g
carries 3p/ρ, a thousand times the physical pressure in the light phase, and
on a coarser level (Lax–Wendroff streaming, coarse/fine interpolation) the
buoyancy-driven pressure gradient inside the bubble becomes unstable after
about one time unit. With the light-phase criterion the adaptive run agrees
with a uniform finest-level run (`-DAMR_MIN_LEVEL=5`) to within 1% in the
maximum velocity up to t = 2. The kernels were checked term by term against
Eqs. (5)–(35) of the paper. The two models give visibly different bubbles at
t = 8 (the PRE 2017 bubble rises less and sheds no satellite droplets).

`Blk_P` holds the physical pressure for both models, so the output, restriction
and prolongation are the same.

## Output

One output writer is selected at start-up with the environment variable
`AMR_OUTPUT`: `vtk` (default), `vtkhdf`, `tecplot` or `none`, e.g.
`AMR_OUTPUT=tecplot ./build/amr_lbm`. All writers have the same interface
(`Write_Step(Time)`); `Output_Mod` points a procedure pointer at the chosen one.

* `vtk` – VTK XML (binary) files in `vtk/`, with ParaView time series
  `Rising_Bubble_2D.pvd` (all leaf cells with `C`, `P`, `Velocity`, `Level`,
  `Block`) and `AMR_Blocks.pvd` (one quad per leaf block with `Level`, `Block`,
  `Parent`, `dx`: the adaptive mesh structure; show it as *Wireframe* coloured
  by `Level`).
* `vtkhdf` – the same two data sets as single time-dependent VTKHDF files
  `Rising_Bubble_2D.vtkhdf` and `AMR_Blocks.vtkhdf`, written with the
  [fVTKHDF](https://github.com/nncarlson/fvtkhdf) library (requires
  `-DAMR_WITH_FVTKHDF=ON -DCMAKE_PREFIX_PATH=<fVTKHDF install prefix>`, a
  serial fVTKHDF build made with the same compiler, and ParaView 5.13 or later
  to read).
* `tecplot` – the original TecPlot ASCII files `Rising_Bubble_2D_NN.dat`.

## Plotting and animation

`tools/animate_vtk.py` plots or animates the VTK/VTKHDF output with only
numpy and matplotlib (h5py for `.vtkhdf`):

```
python3 tools/animate_vtk.py Rising_Bubble_2D.pvd -o bubble.mp4          # C + AMR blocks
python3 tools/animate_vtk.py Rising_Bubble_2D.pvd -f speed -o speed.gif   # |u|
python3 tools/animate_vtk.py Rising_Bubble_2D.vtkhdf --frame -1 -o last.png
```

Run `python3 tools/animate_vtk.py -h` for all options.

## Conservation across levels

The phase field is conserved exactly on a uniform grid, but not across
refinement levels: all levels advance with the same time step, and on a
coarse/fine block edge the two sides compute the streaming flux (Lax–Wendroff
on coarse levels, a perfect shift on the finest) from their own interpolated
ghost cells, so they disagree. The error is largest where the coarse cells
are large, and grows with the slow compression flow that the weakly
compressible model sets up as the bubble rises (proportional to `Gy_lu`).

A sweep over `AMR_MIN_LEVEL` (8x8 blocks, 2 threads each; uniform = all
cells on the finest level) shows two different cases:

| case | min level | sum C change | y_c − uniform | run time |
|---|---|---|---|---|
| JCP 2016, paper Eo = 5 (L0 = 512) | 1 / 3 / 4 / 5 | −6.4 / −4.4 / −2.4 / −0.3 % | +13 / +13 / +9 / +1.6 | 181 / 279 / 608 / 2079 s |
| PRE 2017, paper Eo = 20 (L0 = 512) | 1 / 3 / 4 / 5 | +0.7 / +0.7 / +0.8 / +0.8 % | +1.1 / +0.8 / +0.6 / +0.4 | 202 / 303 / 599 / 2086 s |
| JCP 2016, default case (L0 = 256, t* = 8) | 1 / 3 / 4 | −0.4 / −0.1 / +0.3 % | +1.5 / +1.4 / +1.2 | 20 / 36 / 163 s |
| PRE 2017, default case (L0 = 256, t* = 36) | 1 / 3 / 4 | +3.2 / +3.2 / +3.2 % | −0.1 / −0.4 / −0.4 | 76 / 150 / 665 s |

(sum C change as a fraction of the bubble area; the uniform runs took 4277,
4352, 346 and 1170 s.) Only the JCP 2016 paper Eo = 5 case depends on the
coarsest level: there the error sits in the far field, where the heavy phase
is compressed (C = 1 + 3P), and the AMR bubble swells less than the uniform
one (width 187 vs 199 lu). It takes level 5, at half the cost of the uniform
run, to bring it within 2 lu; level 4 halves the error at 3.4 times the cost
of level 1. With PRE 2017 the light phase shrinks steadily, by the same amount
at every coarsest level. Accounting the streaming fluxes through every block
side shows where: same-level edges and walls contribute nothing, edges between
the coarser levels cancel to below 1 lu², and almost all of the change (+72 of
+66 lu² net in the default case; mesh adaptation gives −5) comes from the edges
between the two finest levels, next to the bubble. There the fine side shifts
values interpolated from the coarse side, the coarse side applies
Lax–Wendroff to values restricted from the fine side, and the two fluxes
differ by an amount that follows the small light-phase tail (1 − C) left at
that edge: with `refinTol` = 0.0005 / 0.002 / 0.008 the loss is 1.4 / 3.2 /
11 % of the bubble by t* = 36. The bubble's position is not affected (y_c
within 0.2 lu of the uniform run in all three). The default stays at 1; build with `-DAMR_MIN_LEVEL=5` (L0 = 512)
for a JCP 2016 run that should match the uniform grid. (A refluxing correction
was tried: per-direction corrections at the corners of refined regions are
inconsistent for the diagonal directions, and none of the variants was free of
side effects.)

## Progress output

Twenty times per run the solver prints a progress line: step, t*, number of
leaf cells, the speed over the last interval in MLUPS (million lattice-cell
updates per second, total and per thread; every leaf cell counts once per
step, whatever its level), the elapsed time and the estimated time left from
the average time per step so far. The estimate follows the number of refined
blocks, so it drifts as the mesh adapts.

## Time step and OpenMP

One time step, starting from h, g and C after streaming:

1. fill the ghost layers of C (as many as the model's collision reads: two for
   JCP 2016, one for PRE 2017); they also serve the refinement criterion and
   the moments of g
2. adapt the mesh; only if it changed are the ghost layers filled again
3. on output steps, compute the moments of g first and write the state
4. in one parallel region: moments of g and collision in one block loop
   (sharing the chemical potential, and for PRE 2017 the gradient of C),
   ghost cells of h and g, then streaming together with C = sum(h)

The stages of step 4 share work with orphaned `!$OMP DO` constructs whose
implicit barriers order the stages. Leaf blocks are kept sorted by level, so
each ghost-cell pass only visits the blocks of its level. Mesh adaptation
(flagging is parallel) runs between the parallel regions. The block kernels take
the block size from `Config_Mod` (compile-time constants) instead of as
dummy arguments, so their loop bounds and array shapes are known to the
compiler even when they are not inlined across modules: this made streaming
plus C = sum(h) 42–48 % faster with gfortran and about 25 % with ifx (8–12 %
of the whole step). LTO (`-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON`) gives
gfortran the same gain but nothing more, and makes the PRE 2017 collision
slower with ifx, so it is not enabled. Three cell loops carry an `!$OMP SIMD`
hint (with `PRIVATE` for their scalar temporaries): the two direction sums in
PRE 2017's `Moments_g` and the per-cell and equilibrium loops of the JCP 2016
collision. They make the ifx collision 14–20 % faster and do not change
gfortran's results or speed; on the other kernels (streaming, finite
differences, the other model loops) such hints were neutral or slower with
one of the two compilers, so they are left out.

The ghost cells of h and g are interpolated from coarser and finer
neighbours in small strips (3 x 8 cells, 8 directions). The interpolation is
written once, as a template (`src/Biquadratic.inc`), and instantiated with
the strip shape fixed at compile time for the three shapes used by
`Ghost_hg_Mod` (left/right edges, bottom/top edges, corners) and with
run-time sizes for the ghost cells of C. Both ghost modules have one routine
per strip shape for the restriction and interpolation, with the side as a
sign (the strips across an edge run outwards from a fixed start with that
stride); the same-level copies of h and g, the hot path, are written out per
direction with constant sections. A block whose neighbours are all leaves
or coarser gets only the directions its own streaming reads (3 of 8 per side
on the finest level, 6 of 8 with Lax-Wendroff); all 8 only next to a refined
neighbour, whose children interpolate from its ghost cells (ghost cells of h
and g 7 % faster on 1 thread, 12 % on 6, results unchanged).

Things that did not help (measured, gfortran 16, i7-11700K): 512-bit vectors
for the cell kernels (`-mprefer-vector-width=512`; a row of 8 cells per
vector) make the collision ~20 % faster on 1 thread, but on 6 threads the
cores clock down from 4.5 to 4.05 GHz and the whole step gets slower; fusing
the moments of g into the PRE 2017 collision loop makes it 14 % slower
(register pressure in the already long loop body); and dropping the barrier
between the levels of the ghost-cell fills does not make them faster (the
waiting is load imbalance), so a two-phase fill would not pay off. The
JCP 2016 collision loop spills registers (about 35-50 vector stores to the
stack per pass of ~450 instructions; PRE 2017: 7), but splitting it into an
h and a g pass made it 2 % slower: the spills stay in L1.

The ghost cells of h and g are limited by data movement, not by
instructions: per perf they take ~19 % of the cycles but ~12 % of the
instructions, ~40 % of the L1-miss stall cycles and ~45 % of the stalls on a
full store buffer (each left/right ghost value is in its own cache line),
and on 6 threads a third of the reads of lines modified by another core.
Two attempts at that did not pay off: filling the ghost cells and streaming
block by block while the block is in L1 (the streaming then writes into a
second h/g pool, as the neighbours still read the unstreamed values; no
copy-back) was 8 % slower on 1 thread, as the doubled data no longer fits
the caches as well; and copying the left/right columns with an AVX-512
gather and scatter (a C kernel with intrinsics) was neutral on 1 thread and
10 % slower on 6, the cores clocking down for the 512-bit instructions.
This CPU (Rocket Lake) executes 512-bit FMAs on one fused port, so AVX-512
has no more flops per cycle than AVX2, which the build uses throughout. An
AVX2 version (two 4-wide gathers, scalar stores as AVX2 has no scatter) made
the ghost cells 8-10 % slower on 1 thread: gfortran's scalar moves are as fast
and inline, the C kernel does not.

The refinement criteria reduce over a block with a running maximum
(minimum) per cell of a row and one short reduction at the end; gfortran's
MAXVAL/MINVAL follow the IEEE NaN rules and do not vectorize (the PRE 2017
criterion went from 4 % to 2 % of the run). The final division by 64 or 1024 is
folded into the weights (powers of two, so the results are unchanged). This
made the ghost cells of h and g about 27 % faster with both compilers (10 %
of the whole step on one thread); a gprofng line profile showed the generic
version spending a fifth of its time in a separate division pass and most of
the rest in loop overhead of the short loops. Most of the remaining ghost
time is copying the left/right edges, which touches one cache line per value
in the current layout.

## GPU solver (SYCL)

`gpu/` contains a C++/SYCL version of the solver for the JCP 2016 model and
the rising bubble, in single precision (Intel Arc A-series GPUs have no FP64
hardware):

    cmake -B build-gpu -DAMR_BUILD_SYCL=ON -DCMAKE_CXX_COMPILER=icpx -DAMR_PRECISION=4
    cmake --build build-gpu
    ./build-gpu/amr_lbm_sycl [steps [output_interval]]     # AMR_CASE as above

The host keeps the quadtree and runs the same mesh logic as the Fortran code
(slot allocation, Z-order, neighbour table, flagging); all field work runs
on the device (`gpu/lbm_kernels.hpp`, one kernel per Fortran routine). The
same build directory has the Fortran programs in FP32, for comparison: on
the default case the diagnostics (block counts, U_max) are identical and
mass and bubble position agree to about 1e-5. The host thread waits for the
GPU several times per step; on a loaded machine pin it to an idle CPU.
`AMR_GPU_SYNC=1` gives per-stage times. No field output yet.

## GPU offload prototype

`offload/` contains a stand-alone OpenMP target-offload version of the
collision kernel (about half of the time step) on a flat block-pool data
layout, which compares device results and timings against the host. See the
header of `offload/collision_offload.f90`. It uses the solver's `Kinds_Mod`
(with local copies of the D2Q9 constants) and follows `AMR_PRECISION` and
`AMR_BLOCK_SIZE`. Intel Arc A-series GPUs have no
FP64 hardware, so double precision there relies on driver emulation.
