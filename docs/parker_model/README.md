# Parker's multicomponent LBM on the AMR code: analysis

Source: `parker-multiphase-lbm-20261118.zip` (transcribed thesis codes of
J. M. Parker 2008 / Parker & Jovanovic 2009, the Matlab moment scripts and the
central-moment paper draft `latex/moments.tex`). Target: the refactored AMR
code of `claude/block-pool-prototype` (this branch was fast-forwarded to it).

Files here:

| file | contents |
|---|---|
| `parker_surface_tension.py` | sympy: basis check, scaled moments, second-order analysis, surface tension, Gibbs–Duhem, Galilean-invariance defects (sections 1–6) |
| `parker_check.py` | NumPy implementation of the 2D model (periodic, MRT as in the thesis code, Margules) with flat-interface, Laplace and oscillating-droplet tests; drift variants A/B/C |
| `linkcheck.py` | check of the exact discrete link balance (section 3 below) |
| `droplet_dispersion.py` | exact linear viscous theory of 2D droplet shape oscillations (two fluids, equal ρ and ν) |
| `fit_oscillation.py` | fits `parker_check.py osc` histories; σ from the frequency through the exact theory; dynamic σ_eff |

## 1. The model in moment form

Each species *i* has a D2Q9 distribution; its equilibrium is that of a fluid
moving with the species velocity **v**ᵢ = **u** + **u**′ᵢ, where **u** is the
mixture velocity and **u**′ᵢ the diffusion velocity. In the code,

 **u**′ᵢ = −D ∇ln γᵢ + D ∇c − Σⱼ xⱼ(−D ∇ln γⱼ + D ∇c),  D = c_s²(1/ω_D − ½).

The species momentum (and the energy flux) relaxes with ω_D, the stress with
ω_ν. The model's κ multiplies the u′u′ terms of the second-order equilibrium.

sympy (section 1) shows:

* Parker's equilibria e₂…e₉ are exactly the scaled raw moments of the draft
  with k₄ = k₅ = k₆ = κ. The third- and fourth-order rows are those of the
  second-order Lallemand–Luo equilibrium (m₂₁ = cᵢvᵧ/3, no cubic terms), so
  there is no k₇…k₉ in the original model.
* The ∇c term (comment: ∇ln c; code: ∇c, with a typo `c(x-1,y-1)` for
  `c(x+1,y-1)` in the x component) is the same for every species and **cancels
  exactly** in the mass-average correction: u′₁ = −D(1−x₁)(g₁ − g₂). It is dead
  code, so the typo is harmless.
* For Margules (and any Gibbs–Duhem-consistent γ), Σ xⱼ∇ln γⱼ = 0, so
  **u**′ᵢ = −D∇ln γᵢ exactly.
* γ = xᵃ (the default `a = 2` of the samples) has a monotonic μ: it does not
  phase-separate. It is only a non-ideal diffusion test.

## 2. Surface tension in terms of the moments

Summed over species with Σcᵢ**u**′ᵢ = 0, the second moments are

 Σᵢ Π^eqᵢ = c c_s² **I** + c **uu** + **T**,  **T** = κ Σᵢ cᵢ **u**′ᵢ**u**′ᵢ.

**T** is a Korteweg-type capillary stress (u′ ∝ ∇ composition). Second-order
expansion of the steady LB for a flat interface at rest (normal *n*) gives
(sympy, section 3; no ω_ν, ω_e, ω_ε survive):

 Jᵢ = cᵢu′ᵢ − (1/ω_D − ½) ∂ₙ[cᵢ(c_s² + κu′ᵢ²)] (species flux)
 Π_nn − Π_tt = κ Σᵢ cᵢ u′ᵢ,ₙ², ∂ₙΠ_nn = 0

so the mechanical surface tension is

 **σ = κ ∫ Σᵢ cᵢ u′ᵢ,ₙ² dn = κ D² ∫ Σᵢ cᵢ (∂ₙ ln γᵢ)² dn**,

and with the zero-flux condition at leading order (u′ᵢ = D∂ₙ ln cᵢ) for a binary
mixture with c ≈ const and x₁ = sin²θ:

 σ = κ D² c ∫ x₁′² / (x₁x₂) dn = 4 κ D² c ∫ θ′² dn ≈ 4κD²c (Δθ)²/ℓ.

Consequences, all confirmed numerically (`parker_check.py flat`):

* σ_mech = ∫(Π_xx − Π_yy) dx equals κΣcᵢu′ᵢ² to all printed digits. Π_xx is
  constant to 1e−15, and κ = 0 gives σ = 0 (yet the phases still separate).
* **The steady state depends on κD² only.** (D, κ) = (0.05, 1) and (0.1, 0.25)
  give identical profiles and σ. This is the group behind
  `kappa = E_sigma/36/diffco**2` in `periodic.f90`: E_σ = 36κD².
* **Control relation:** σ = κD²c·I. The interface is grid-scale (≈3 cells for
  Margules A₁₂ = 2.76, A₂₁ = 2.45), and I ≈ 0.81 (κD² = 0.0025), 0.79 (0.01),
  0.72 (0.04), 0.64 (0.09). The drop comes from the κ term in the species
  flux, which moves the profile.
* The leading-order closed form κD²Σ(cᵢ′)²/cᵢ is 15–20% low: the interface is
  too thin for the leading-order zero-flux relation (see §3).

| test (D = 0.2, κ = 1, ν = 0.1) | σ |
|---|---|
| flat interface, ∫(Π_xx − Π_yy) | 0.02899 |
| Laplace, R = 11.1 / 17.5 / 23.7 | 0.0298 / 0.0294 / 0.0293 (+2.8 / +1.5 / +1.2 %) |
| oscillating droplet, n = 2, R = 20, ν = 0.02, against the exact viscous theory | 0.0297 (+2.5 %) |

### Oscillating droplet

The reference is the exact linear normal mode of a 2D droplet in an unbounded
fluid of the same density and viscosity (`droplet_dispersion.py`: potential
flow plus I_n/K_n vortical parts; continuity of velocity and shear stress; a
normal-stress jump σ(n²−1)η/R²). It recovers the inviscid ω₀ as ν → 0. With
two fluids the damping comes from the interfacial boundary layer
(δ = √(2ν/ω) ≈ 3.6 lu here), so the textbook pair ω₀² = 6σ/(2ρR³),
γ = 4ν/R² is wrong by much more than the model error: at ν = 0.02 the
exact frequency is 9 % below ω₀ and the damping is 1.7× 4ν/R². That
accounted for most of the "−12.5 %" first reported.

`fit_oscillation.py` turns each fitted frequency into σ through the exact
theory. It also prints the dynamic surface tension
σ_eff(t) = κΣ∫cᵢ|u′ᵢ|²dA / perimeter, which equals σ_mech for a static interface.
All runs: R = 20, ε = 0.08, Margules, drift A, static σ_flat = 0.02899 (same
κD² = 0.04 except where noted):

| run | ω vs exact theory | σ_osc | σ_eff (mean, swing) | γ / γ_theory |
|---|---|---|---|---|
| D = 0.2, ν = 0.02 (box 112) | +1.3 % | +2.5 % | +0.8 %, 0.8 % | 1.41 |
| D = 0.2, ν = 0.02 (box 176) | +1.7 % | +3.2 % | −0.4 %, 1.0 % | 1.44 |
| D = 0.2, ν = 0.005 | +1.2 % | +2.5 % | +0.6 %, 0.8 % | 2.46 |
| D = 0.2, ν = 0.05 | +1.2 % | +2.2 % | +1.1 %, 0.8 % | 1.19 |
| D = 0.4, κ = 0.25 | +4.4 % | +8.5 % | +2.1 %, 1.4 % | 1.64 |
| D = 0.025, κ = 64 | −3.7 % | −6.9 % | −1.2 %, **8.0 %** | 1.33 |

Diagnosis:

1. **For D = 0.2 the surface tension is right.** σ_osc is 2.2–3.2 % above
   σ_flat at every viscosity and both box sizes, the same as the Laplace test at
   this radius (+1.2 to +2.8 %): a small curvature/finite-width effect, not a dynamic
   one.
2. **Small D: the interface lags the flow.** Its composition relaxes at
   about D/w² (w ≈ 3 cells). At D = 0.025 that is 2.8e−3, the same as ω, and
   σ_eff swings by 8 % over a period instead of < 1.5 %. The interface
   behaves viscoelastically and the frequency drops. Rule:
   **D/(w²ω) ≫ 1** (≈ 7 at D = 0.2, ≈ 1 at D = 0.025). With a grid-scale interface
   this ties D to the capillary time scale, which is another reason for the
   resolved width of §3.1.
3. **Large D: the frequency rises (+4.4 % at D = 0.4)** more than σ_eff
   does (+2.1 %). **This part is the Galilean-invariance error of Parker's
   equilibrium.** With the collision in central Hermite moments
   (`parker_check.py ... coll=chm`) it drops to +2.6 %, which matches the
   +2.0 % of σ_eff. The damping drops from 1.64× to 1.46× theory. The static
   state is bit-identical to the MRT version, and the spurious velocities fall
   from about 1e−11 to 1e−16. The source term for the aliased xxx/yyy
   moments (`coll=chmg`) changes nothing measurable at these velocities
   (|u| ≲ 5e−3): it can wait for flows at higher Mach numbers.

   | collision | D = 0.4 | D = 0.2 | D = 0.025 |
   |---|---|---|---|
   | Parker MRT | +4.4 % / ×1.64 | +1.3 % / ×1.41 | −3.7 % / ×1.33 |
   | central Hermite | +2.6 % / ×1.46 | +1.0 % / ×1.36 | – |
   | central Hermite + xxx/yyy source | +2.6 % / ×1.46 | +1.0 % / ×1.36 | −3.8 % / ×1.35 |

   (frequency relative to exact theory / damping relative to exact theory)
4. **Extra damping** of 1.1–2.3e−4 on top of the exact viscous value. It grows
   with D (diffusive dissipation in the interface, as in any
   diffuse-interface model) and is largest at ν = 0.005, where τ_ν = 0.515
   and the boundary layer (δ ≈ 1.8 lu) is not resolved.

## 3. The discrete balance and the stencils

For a steady 1D state at rest, the zero link-flux condition of the lattice
reduces **exactly** to

 ½[(cᵢu′ᵢ)ₓ + (cᵢu′ᵢ)ₓ₊₁] = (D/c_s²)(Π^eqᵢ,ₓ₊₁ − Π^eqᵢ,ₓ),  Π^eqᵢ = cᵢ(c_s² + κu′ᵢ²)

(a link average of the nodal drift against a compact two-point difference).
`linkcheck.py` confirms it to 2e−10, independent of ω_ν. The nodal
central-difference form misses by 18% of the signal. Parker's drift uses the
wide isotropic G stencil of ln γ at the nodes, so on a 3-cell interface the
discrete equilibrium is not the thermodynamic one:

* coexistence (0.089, 0.870) instead of the Margules values (0.1005, 0.8564)
* |∇ln a₁| ≈ 0.09 left at the interface (|∇ln γ₁| ≈ 0.78)

Drift variants tested (`parker_check.py flat D kappa A|B|C`):

| drift | u′ᵢ before the mass-average correction | result |
|---|---|---|
| A (Parker) | −D G[ln γᵢ] | coexistence off; I falls 0.81 → 0.64 as κD² grows |
| B | −D G[ln aᵢ] + D G[cᵢ]/cᵢ | coexistence worse (0.068, 0.884) |
| C | −D G[ln aᵢ] + D/(c_s²cᵢ) G·Π^eqᵢ (lagged) | **I = 0.78–0.81 over κD² = 0.0025…0.09**: σ ≈ linear in κD²; coexistence still off |

So the κ-compensation of C does what it was designed for. Writing the drift in
terms of activities with wide nodal stencils does not fix the composition
error, because the error is the link-average/compact-difference mismatch
above. Ways to remove it (not yet tested):

1. **Resolved interface:** add a square-gradient term to the chemical
   potential (ln aᵢ − λ∇²xᵢ, Cahn–Hilliard-like) so that W ≈ 4–6 cells. The
   mismatch then drops to O(Δx²/W²), and σ ≈ 4κD²c(Δθ)²/W becomes grid-convergent
   and set by λ. This is the prerequisite for meaningful AMR (see §5).
2. **Link-wise drift:** inject the drift as a flux on the links (antisymmetric
   post-collision corrections of the populations crossing a link) with the
   logarithmic-mean mobility c̄ = Δc/Δln c. Then cᵢΔln cᵢ = Δcᵢ holds exactly
   on every link, and zero flux ⇔ Δln aᵢ = 0: exact discrete phase equilibrium.
   The link balance above has to be redone with the added link flux (sympy, 1D,
   as in section 3) before implementing.
3. Keep C's compensation term (or put it into the link drift) so that σ stays
   linear in κ.

## 4. Galilean invariance and the collision space

sympy section 6 compares the mixture's third moments with those of the target
(per species a Gaussian with mean **v**ᵢ and covariance c_s²**I** + (κ−1)**u**′ᵢ**u**′ᵢ):

* Parker's equilibrium: u-dependent errors in all of Q_xxx, Q_xxy, Q_xyy, Q_yyy.
* **Central Hermite moment equilibrium** about **u** (a₀ = cᵢ, a₁ = cᵢ**u**′ᵢ,
  a₂ = κcᵢ**u**′ᵢ**u**′ᵢ, a_xxy = κ₃cᵢu′ₓ²u′ᵧ, …, all free of **u**): Q_xxy, Q_xyy
  become exact. Only the aliased Q_xxx, Q_yyy (D2Q9: c_x³ = c_x) remain:

 δQ_xxx = c uₓ³ + 3 uₓ T_xx

  i.e. the standard cubic defect **plus an interface term**: the capillary
  stress advected by **u**. It gives a viscous stress error
  −(1/ω_ν − ½)∂ₓδQ_xxx (moving or oscillating droplets).
  Correction: a source S_xx = −∂ₓδQ_xxx (and S_yy) in the normal-stress moments,
  applied with (1 − ω_ν/2) in the collision. ∂ₓuₓ can come from the local
  non-equilibrium stress; ∂ₓT_xx needs a second derivative of ln γ (two ghost
  layers of C, as the JCP 2016 model already has).

Recommended collision: central Hermite moments about the mixture velocity,
one set of rates for both species (ω_D for a₁, ω_ν/ω_b for a₂, free rates for
a₃, a₄), plus the xxx/yyy correction. The draft paper's κ₄…κ₉ reduce to κ (a₂),
κ₃ (a₃) and κ₄ (a₄). At rest only κ enters σ (section 3).

## 5. Refinement indicator

* ∇μᵢ = ∇ln aᵢ is ~0 at an equilibrium interface (0.09 here, from the
  discretization alone, and it vanishes in a consistent scheme). As an indicator
  it would **derefine converged interfaces**. It only flags interfaces out of
  equilibrium.
* ∇ln γᵢ is the right quantity: zero in the bulk, maximal at the interface, and
  ∝ u′ᵢ, the quantity σ integrates. At equilibrium it equals −∇ln xᵢ (a
  log-weighted composition gradient, more sensitive on the dilute side than |∇C|).
  Suggested block error (finest-level units, as `Max_Gradient`):
  max over cells and species of |∇ln γᵢ|. Optionally also max |∇ln aᵢ|, to keep
  diffusion fronts (dissolution, spinodal decomposition) refined.
  Alternative with a physical tolerance: T/(c c_s²), the capillary stress
  relative to the pressure.
* Because the Margules interface is grid-scale, its width in cells is the same
  on every level. With the same Δt on all levels (D, ν, κ level-independent in
  finest units) σ on level L would be σ_finest·Δx_f/Δx_L. **The interface
  must always sit on the finest level**, or the model needs the resolved width
  of §3.1.

## 6. Adding the model to the AMR code

The solver already separates models (`LBM_Model` procedure pointers) and cases
(`AMR_Case`), so the model fits in as a new module with little infrastructure
change:

| piece | approach |
|---|---|
| distributions | species 1 in `Blk_h`, species 2 in `Blk_g` (same shapes, ghosts, streaming, prolongation, restriction) |
| scalars | `Blk_C` = c₁ (already sum(h) after streaming). New pool `Blk_C2` = c₂ = sum(g) with the same ghost layers, allocated only for two-scalar models; the solver sums g after streaming only then |
| ghost cells | `Fill_Ghost_Cells_C` takes the pool array as an argument (`Blk_C` or `Blk_C2`) instead of using `Blk_C` from the module; same for the prolongation/restriction of C |
| model module | `Model_Parker2009_Mod` with `Parker_Model()`: `Moments_g` → c₂, **u** (mixture), `Blk_P` = c c_s²; `Collide` → ln γ (with clamping x ∈ [ε, 1−ε]), G stencil scaled by `inv_dx_level`, u′ᵢ, central Hermite collision of both species; `Refinement_Error` → max |∇ln γᵢ| |
| activity models | a small `Activity_Mod` (Margules, UNIQUAC from `uniquac_mod.f90`) selected by the case |
| case | `Droplet_Mod`: static droplet (Laplace), n = 2 oscillation, later the rising droplet of `TwoDSampleCode_gravity.f90` |
| walls | `Symmetric_C` (zero normal gradient, as Parker) and bounce-back for both species |
| levels | same Δt, Lax–Wendroff on coarse levels and unscaled populations work here too: τ and D are level-independent in finest units. Interfaces stay on the finest level (§5) |
| performance | the JCP/PRE models are unchanged; `Blk_C2` and its ghost pass exist only for the new model. Check with the usual 1+6-thread runs of both models |

Suggested order: (1) the generic ghost/prolong routines + `Blk_C2`, with timings
of the two existing models; (2) Parker's original MRT model on a uniform grid
(`AMR_MIN_LEVEL = maxLevels`), reproducing the flat/Laplace numbers above;
(3) AMR with the ln γ indicator; (4) central Hermite collision and the xxx/yyy
correction; (5) the stencil fix of §3 (resolved interface or link-wise drift).

## 7. The AMR implementation and its checks

Implemented as planned in §6 (`src/Model_Parker2009_Mod.f90`,
`src/Parker_Cases_Mod.f90`, driver `amr_lbm_parker`). The collision is in central
Hermite moments without the xxx/yyy source (§4 showed no effect at these
velocities). The criterion is max |∇ln γᵢ| (refine 0.02, derefine 0.01).
Changes to the shared code:
`Blk_C2` with ghost layers, prolongation and restriction (only for
`n_scalars = 2`); `Fill_Ghost_Cells_C` with the scalar pool as an argument; the
optional refinement of coarse islands (`fill_holes`).

| check (L0 = 128 unless noted, D = 0.2, ν = 0.02) | result |
|---|---|
| Laplace, R = 20: adaptive vs uniform finest level | Δp 1.510e−3 vs 1.516e−3 (0.4 %) |
| oscillation, R = 20: adaptive vs uniform | identical; ω +0.5 % vs exact theory |
| cost, R = 20 oscillation (12000 steps) | 5.7 s on 2 threads vs 16.4 s on 4 threads uniform |
| JCP 2016 / PRE 2017 timings vs 283bacf (1 and 6 threads, best of 3) | ratio 0.96–0.99: no slowdown |

Convergence series: R = 10, 20, 40 at L0 = 64, 128, 256, κ = 2, 1, 0.5 (σR
nearly constant). σ from flat interfaces at the same κD². Oscillation minus
the static run with the same R:

| R | ω vs exact theory | γ / theory | static Laplace σ vs flat |
|---|---|---|---|
| 10 | +0.4 % | 1.52 | +11.6 % (the droplet dissolves: R 10.2 → 9.1) |
| 20 | +0.3 % | 1.24 | +2.7 % |
| 40 | −0.9 % (island refinement) | 1.11 | −0.3 % |

Without island refinement the R = 40 droplet had a 2×2 island of level-4
blocks inside it (its interior has no gradient to flag). Its flow crossed
coarse/fine edges inside the dense phase, where the same-Δt scheme is not
conservative: +2.3e−3 of the species-1 mass over the run, the droplet grew from
R = 40.04 to 40.27, and ω was +2.25 % off. Raising the coarsest level to 3 or 4
changed nothing. A uniform finest-level run was too slow to wait for. With the
islands refined, the drift drops tenfold (−2.3e−4) and ω is −0.9 %.

So the grid-scale interface converges both statically and dynamically, and the
resolved interface of §3.1 is not needed for that. It remains the way to fix
the composition error (coexistence 0.108 / 0.853 instead of 0.101 / 0.855).

Kinetic sensor (refining on |f_neq|): from §3, at a steady interface the
species momentum is j = −ω_D/(2 − ω_D)·cᵢu′ᵢ, so the species f_neq in the
first moment is ∝ cᵢu′ᵢ: largest in equilibrium interfaces, unlike ∇μ. The
mixture's deviatoric second moment minus c(**uu**)_dev is **T** plus the
viscous stress, so a sensor on it would flag interfaces and shear layers from
the local populations alone (no stencil, no ghost cells). Not implemented yet.
