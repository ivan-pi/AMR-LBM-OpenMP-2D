# AMR-LBM-OpenMP-2D
Adaptive Mesh Refinement (AMR) + Lattice-Boltzmann Method (LBM) + OpenMP

Based on the following papers:

A. Fakhari, D. Bolster, L.-S. Luo, “A weighted multiple-relaxation-time lattice Boltzmann method for multiphase flows and its application in partial coalescence cascades”, Journal of Computational Physics 341:22–43 (2017) [https://doi.org/10.1016/j.jcp.2017.03.062].

A. Fakhari, M. Geier, T. Lee, “A mass-conserving lattice Boltzmann method with dynamic grid refinement for immiscible two-phase flows”, Journal of Computational Physics 315:434–457 (2016) [https://doi.org/10.1016/j.jcp.2016.03.058].

A. Fakhari, T. Mitchell, C. Leonardi, D. Bolster, “Improved locality of the phase-field lattice-Boltzmann model for immiscible fluids at high density ratios”, Physical Review E 96:053301 (2017) [https://doi.org/10.1103/PhysRevE.96.053301] (second model; uniform-grid reference code: https://github.com/AbbasHub/Improved_Phase-Field_LBM_2D).

This package is for numerical simulation of multiphase flows (rising bubble,
Rayleigh–Taylor instability).
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
(tested with gfortran 13, Intel ifx 2026.1 and NAG Fortran 7.2):

```
cmake -B build                                  # Release, double precision
cmake --build build -j
OMP_NUM_THREADS=8 ./build/amr_lbm [number_of_time_steps [output_interval_in_steps]]
OMP_NUM_THREADS=8 ./build/amr_lbm_pre2017 [...]    # same case, PRE 2017 model
AMR_CASE=drop.nml ./build/amr_lbm_parker [...]     # two species (see Two-species model)
```

Five drivers are built, one per case and LBM model (see *Models*): `amr_lbm`
(JCP 2016) and `amr_lbm_pre2017` (PRE 2017) for the rising bubble,
`amr_lbm_rt` and `amr_lbm_rt_pre2017` for the Rayleigh–Taylor instability, and
`amr_lbm_parker` for the two-species model of Parker & Jovanovic (2009) with
its own cases.

The block fields live in one pool array per field, allocated at start-up for
the blocks of the fully refined quadtree (`4**(maxLevels+1) - 1`, e.g. 4095
blocks or 67 MiB for the default case). Only the slots that get used are
written, so the operating system commits only those pages. The environment
variable `AMR_BLOCK_CAPACITY` sets a smaller start capacity; the pool then
doubles whenever a refinement needs more blocks.

Configuration options (`cmake -B build -D<option>=<value>`):

| option | default | meaning |
|---|---|---|
| `CMAKE_Fortran_COMPILER` | system default | e.g. `ifx`, `nagfor` |
| `CMAKE_BUILD_TYPE` | `Release` | `Debug` adds bounds checking and warnings |
| `AMR_PRECISION` | `8` | working precision `wp`: `8` (real64) or `4` (real32) |
| `AMR_L0` | `256` | base block width in finest lattice units |
| `AMR_DOMAIN_HEIGHT` | `3` | domain height in base blocks (domain `L0 x AMR_DOMAIN_HEIGHT*L0`; 4 for the Rayleigh–Taylor paper case) |
| `AMR_MIN_LEVEL` | `1` | minimum (uniform background) refinement level |
| `AMR_BLOCK_SIZE` | `8` | cells per block in each direction (`16` is faster for larger cases) |
| `AMR_PERIODIC_Y` | `OFF` | periodic in y as well (default: walls at the bottom and top) |
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

### Rayleigh–Taylor instability

A heavy fluid above a light one, the interface at mid-height perturbed by
y = 0.1 L0 cos(2πx/L0): the single-mode case of Fakhari, Mitchell, Leonardi &
Bolster (PRE 2017), Sec. IV B, Figs. 5–6. The drivers are `amr_lbm_rt` (JCP 2016
model) and `amr_lbm_rt_pre2017`; the paper's L × 4L domain needs
`-DAMR_DOMAIN_HEIGHT=4`:

```
cmake -B build-rt -DAMR_L0=256 -DAMR_DOMAIN_HEIGHT=4 && cmake --build build-rt -j
AMR_CASE=cases/paper2017_RT_At05.nml ./build-rt/amr_lbm_rt
python3 tools/plot_rt.py rt_fronts.dat:JCP2016 --paper -o rt.png
```

The namelist `&rayleigh_taylor` sets `At` (Atwood number), `Re` = ρh U0 L0/μh,
`mu_ratio`, `Ca` = μh U0/σ, `Pe` = U0 L0/M (or `Mobility`), `Cn` / `Width`,
`Gy_lu`, `Amplitude` (/L0), `t_end` and `refine_light_phase`, with U0 = √(g L0) and the time
t* = t/t0, t0 = √(L0/(g At)). The defaults are the paper's At = 0.5 case
(Re = 3000, μh/μl = 1, Ca = 0.26, Pe = 1000, W = 5, t0 = 16000 steps at
L0 = 256) but end at t* = 2.5, when the spike nears the bottom wall of the
default L0 × 3L0 domain. The light fluid starts with the hydrostatic pressure
of the models' buoyancy force (ρh − ρ) g. `rt_fronts.dat` records the bubble
front (x = 0) and the spike tip (x = L0/2), measured from the unperturbed
interface in units of L0.

Paper case (L0 = 256, l_max = 5, min level 1, 3 threads each on a loaded machine):

| t* | 1 | 1.5 | 2 | 2.5 | 3 |
|---|---|---|---|---|---|
| bubble, PRE 2017 Fig. 6 | 0.29 | 0.425 | 0.555 | 0.685 | 0.835 |
| bubble, `amr_lbm_rt` | 0.289 | 0.428 | 0.561 | 0.696 | 0.841 |
| bubble, `amr_lbm_rt_pre2017` | 0.287 | 0.423 | 0.554 | 0.686 | 0.828 |
| spike, PRE 2017 Fig. 6 | −0.36 | −0.61 | −0.86 | −1.11 | −1.43 |
| spike, `amr_lbm_rt` | −0.347 | −0.593 | −0.825 | −1.066 | −1.372 |
| spike, `amr_lbm_rt_pre2017` | −0.350 | −0.600 | −0.839 | −1.085 | −1.387 |

The paper values are read off Fig. 6 (its "current" curve, about ±0.01). Both
models reproduce the bubble front within 0.01 L0 and the shapes of Fig. 5 (the
mushroom, the roll-ups and the satellite drops at t* = 3). The spike tip trails
the paper by 0.01–0.06 L0 (up to 4 %) late in the run, about the spread between
the references in Fig. 6. With the JCP 2016 model, a uniform finest-level run
(`-DAMR_MIN_LEVEL=5`) gives a spike of −0.836, −1.083, −1.394 at t* = 2, 2.5, 3:
the AMR run trails the uniform one by up to 0.02 L0, and min level 3 is
the same as min level 1 to 0.001 L0. The total mass changes by −0.015 % (AMR)
and 0 (uniform) by t* = 3.2.

`refine_light_phase` (default `.false.`) selects the refinement criterion of
the PRE 2017 model in `amr_lbm_rt_pre2017`. The rising bubble needs the light
phase on the finest level (see *Models*), but at At = 0.5 the interface alone
suffices: the fronts agree with the light-phase criterion to within
0.004 L0 (bubble) and 0.001 L0 (spike), at about half the run time (196 vs 413 s
side by side, 3 threads each; 280 instead of 2188 leaf blocks at t* = 0). The
PRE 2017 row of the table above is from the light-phase criterion. Set it to
`.true.` for high density ratios (e.g. the paper's At = 0.998 case).

## Source layout

| file (`src/`) | contents |
|---|---|
| `main_jcp2016.f90`, `main_pre2017.f90`, `main_rt_*.f90`, `main_parker2009.f90` | driver programs: select the model and the case |
| `Solver_Mod.f90` | the solver: initialization and time loop for any model and case, progress and timers |
| `Case_Mod.f90`, `Rising_Bubble_Mod.f90`, `Rayleigh_Taylor_Mod.f90`, `Parker_Cases_Mod.f90` | a test case as a set of procedures (setup, initial state, diagnostics, finish); the rising bubble, the Rayleigh–Taylor instability, droplets and spinodal decomposition of the two-species model |
| `Model_Mod.f90` | an LBM model as a set of block kernels (initialize, collide, moments of g) and its refinement criterion |
| `Model_JCP2016_Mod.f90`, `Model_PRE2017_Mod.f90`, `Model_Parker2009_Mod.f90` | the models' kernels and criteria; `Collision_JCP2016.inc` is the template of the JCP 2016 collision kernel |
| `LBM_Common_Mod.f90`, `Phase_Field_Mod.f90` | kernels shared by the models: streaming and C = sum(h); finite differences of C (chemical potential, gradient, largest gradient) |
| `Kinds_Mod.f90`, `Config_Mod.f90`, `D2Q9_Mod.f90`, `Physics_Mod.f90` | kinds, compile-time block size and domain, lattice, physical parameters (set by the case) |
| `Blocks_Mod.f90`, `Mesh_Mod.f90` | block storage: integer metadata `Blk(:)` and one pool array per field (`Blk_h(:,:,:,n)`, `Blk_C(:,:,n)`, ...), leaf lists, free-slot queue; the domain and child blocks, the neighbour table (rebuilt from the block positions after every mesh change) |
| `Refine_Mod.f90`, `Prolong_Restrict_Mod.f90` | flagging, refinement/derefinement (optionally refining coarse islands enclosed by finer blocks); data transfer to children and parents |
| `Ghost_C_Mod.f90`, `Ghost_hg_Mod.f90`, `BC_Mod.f90`, `Interpolations_Mod.f90` | ghost cells of C (1 or 2 layers; and of C2 for two-scalar models) and of h/g, wall conditions, coarse/fine interpolation |
| `Output_Mod.f90`, `VTK_Output_Mod.f90`, `TecPlot_Mod.f90` | writer selection (`AMR_OUTPUT`), VTK/VTKHDF and TecPlot writers |

`code_structure.mmd` is a Mermaid diagram of the module dependencies.
`bench/collision/` holds the benchmark and the performance study of the JCP 2016
collision kernel and of the streaming, `bench/ghost/` those of the ghost cells
of h and g.

Every file holds one module of the same name, so all calls are checked against
explicit interfaces at compile time. Modules are `PRIVATE` by default and list
their public entities explicitly; helpers used by a single routine are internal
procedures of that routine.

The solver knows the model and the case only through the procedure pointers
of `LBM_Model` (`Model_Mod`) and `AMR_Case` (`Case_Mod`). A case may adjust
the model's mesh policy (tolerances, island refinement) through the optional
`Adjust_Model`, called after its setup: a new case is a new
module with a constructor like `Rising_Bubble_Case()`, a new model one like
`JCP2016_Model()`. Each model brings its refinement criterion (`Block_Error`
interface) and tolerances, since the blocks it needs on the finest level depend
on the model; mesh adaptation takes the model (`Refine_Derefine_Mesh(model)`)
and also uses its moments kernel for new blocks. A case returns its end time
and its time unit per lattice time step from `Setup`; the solver passes times in
that unit to the output and to the case's diagnostics. The domain shape
(`nBaseX` x `nBaseY` base blocks of `L0` x `L0`, `nBaseY` from
`AMR_DOMAIN_HEIGHT`) is set in `Config_Mod`.
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
| refinement criterion (chosen by the model) | `Interface_Error`: maximum gradient of C in a block | `Interface_Light_Phase_Error`: the interface, and additionally the light phase (max(1−C) of a block); `PRE2017_Model(refine_light_phase=.false.)` gives `Interface_Error` |

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

### Two-species model (Parker & Jovanovic 2009)

`amr_lbm_parker` runs the multicomponent model of J. M. Parker (PhD thesis,
Oregon State University, 2008; Parker & Jovanovic, Chem. Biochem. Eng. Q. 23,
2009): species 1 in h (`Blk_C` = c₁), species 2 in g (`Blk_C2` = c₂, a second
ghost-cell scalar that the solver allocates, fills and transfers only for
models with `n_scalars = 2`). Each species relaxes towards the equilibrium of
a fluid moving with u + u′ᵢ, where u′ᵢ = −D ∇ln γᵢ (mass-average reference) and
ln γᵢ comes from Margules' equation. Diffusion is driven by the chemical
potential, and the phases separate where the mixture is unstable. The u′u′
terms of the second moments are scaled by κ, so the mixture carries the
capillary stress T = κ Σᵢ cᵢ u′ᵢu′ᵢ and the surface tension is
σ = κ ∫ Σᵢ cᵢ u′²ᵢ,ₙ dn ≈ κ D² c I (I ≈ 0.65–0.81, depending on κD²).
The collision is in central Hermite moments about u (Galilean invariant). The
refinement criterion is max |∇ln γᵢ|: zero in the bulk, largest in
interfaces (∇μᵢ itself vanishes across an equilibrium interface). The model
also refines small coarse islands enclosed by finer blocks (`max_island_blocks`), as the
interior of a large droplet would otherwise be coarser than its interface and
lose mass at the coarse/fine edges inside it.

Parameters come from the namelist `&parker` in the file `AMR_CASE`:

| variable | default | meaning |
|---|---|---|
| `kind` | `'droplet'` | `'droplet'` (Laplace test, or n = 2 oscillation with `eps > 0`; `droplet.dat`), `'spinodal'` (mixture square of L0 × L0; `spinodal.dat`) or `'box'` (the whole domain mixed, e.g. a periodic box; `spinodal.dat`) |
| `D_lu`, `kappa`, `nu_lu` | 0.2, 1, 0.02 | diffusivity, surface-tension scaling, viscosity |
| `A_12`, `A_21` | 2.7627, 2.4474 | Margules parameters |
| `Radius`, `eps` | 20, 0 | droplet radius and elongation |
| `x_mix`, `noise` | 0.5, 0.05 | spinodal/box: mean mole fraction of species 1 and perturbation (`noise = 0.5`: uniform in [0, 1]) |
| `n_steps` | 12000 | number of time steps |
| `refine_tol`, `deref_tol`, `max_island` | model: 0.02, 0.01, 16 | mesh policy overrides (negative: keep the model's value); the refine tolerance must lie below the gradient of ln γ that the initial state shows on the coarsest level, or the mesh never resolves it |

Spinodal decomposition in a periodic box as in Parker's `periodic.f90`
(128 × 128, uniform random x₁, D = ν = 1/3, κ = 1/(36 D²)):

```
cmake -B build-box -DAMR_L0=128 -DAMR_DOMAIN_HEIGHT=1 -DAMR_PERIODIC_Y=ON && cmake --build build-box -j
printf "&parker kind='box', noise=0.5, D_lu=0.3333333, nu_lu=0.3333333, kappa=0.25, n_steps=60000 /\n" > box.nml
AMR_CASE=box.nml ./build-box/amr_lbm_parker 60000 300
python3 tools/animate_vtk.py Rising_Bubble_2D.pvd --no-blocks --clim 0 1 --cmap gnuplot -o box.mp4
```

Checked against the exact linear viscous theory of 2D droplet oscillations
(`docs/parker_model/`): adaptive and uniform finest-level runs agree, the
oscillation frequency is within ±1 % of theory for R = 10, 20 and 40, and the
Laplace pressure converges to the flat-interface σ (−0.3 % at R = 40).
The derivation, the numerical checks and the open points (the composition
error of the grid-scale interface, the stencil mismatch of the drift) are in
`docs/parker_model/README.md`.

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
   each block then writing its edges into the ghost cells of h and g of its
   same-level neighbours; the remaining ghost cells of h and g (coarse/fine
   interpolation, walls); then streaming together with C = sum(h)

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
PRE 2017's `Moments_g` and the cell loop of the JCP 2016 collision. On the
other kernels (streaming, finite differences, the other model loops) such
hints were neutral or slower with one of the compilers, so they are left out.

The JCP 2016 collision (`src/Collision_JCP2016.inc`) computes the chemical
potential, the moments of g and the collision of h and g in one pass over the
cells, and relaxes the stress moments from closed expressions for the pair
sums of the equilibrium instead of forming every g - geq. This form is the
outcome of a study of the kernel (loop structure, bounds, SIMD, symmetry,
fusion and fission, register pressure, block size, divisions, AVX-512
intrinsics) with gfortran, ifx and flang, documented with all measurements in
[`bench/collision/README.md`](bench/collision/README.md): the stage "moments
of g and collision" became 1.3-1.5 times faster, with results equal to
round-off. The same study covers streaming (section 12 there): both streaming
kernels use one plane temporary per direction, 13 % faster with gfortran 14
and ifx and 1.9 times with gfortran 15, results unchanged. The ghost cells of
h and g were studied the same way ([`bench/ghost/README.md`](bench/ghost/README.md)):
the same-level copies are pushed by each block right after its collision,
while it is in the cache, and the interpolation from a coarser block covers
only the directions that are read; collision plus ghost cells became
1.15-1.3 times faster.

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

## GPU offload prototype

`offload/` contains a stand-alone OpenMP target-offload version of the
collision kernel (about half of the time step) on a flat block-pool data
layout, which compares device results and timings against the host. See the
header of `offload/collision_offload.f90`. It uses the solver's `Kinds_Mod`
(with local copies of the D2Q9 constants) and follows `AMR_PRECISION` and
`AMR_BLOCK_SIZE`. Intel Arc A-series GPUs have no
FP64 hardware, so double precision there relies on driver emulation.
