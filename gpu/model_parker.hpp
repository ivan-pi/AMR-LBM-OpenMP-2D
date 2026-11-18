//==============================================================================
//  Multicomponent LBM with chemical-potential driven diffusion on the GPU
//  (src/Model_Parker2009_Mod.f90), for NS = 2 or 3 species
//
//  J. M. Parker, PhD thesis, Oregon State University (2008);
//  J. M. Parker, G. N. Jovanovic, Chem. Biochem. Eng. Q. 23 (2009) 399-409.
//  Analysis and the surface-tension relation: docs/parker_model/README.md
//
//  Species s in the distribution f[s], its concentration c_s = sum f[s] in
//  the scalar C[s] (one ghost layer). Each species relaxes towards the
//  equilibrium of a fluid moving with u + u'_s, u the mixture velocity and
//      u'_s = -D grad ln gamma_s  -  sum_r x_r (-D grad ln gamma_r)
//  (mass-average reference). The u'u' terms of the second moments are scaled
//  by kappa: the mixture carries the capillary stress
//  T = kappa sum_s c_s u'_s u'_s. Collision in central Hermite moments about
//  u (Galilean invariant), the rates of the Fortran model: a1, a3 omega_D,
//  deviatoric a2 omega_nu, trace of a2 and a4 omega_E.
//
//  Activity coefficients: the excess Gibbs energy (units of RT) of the
//  three-suffix Margules (Wohl) form
//      g = sum_{r<s} x_r x_s (A_sr x_r + A_rs x_s) + C123 x_1 x_2 x_3,
//  a homogeneous cubic, so that ln gamma_s = dg/dx_s - 2 g. For two species
//  this is Margules' two-parameter equation, ln gamma_1 =
//  (A12 + 2 (A21 - A12) x1) x2^2, evaluated as in the Fortran code (with x1
//  clamped to [x_eps, 1 - x_eps]); for three species the mole fractions are
//  clamped to [0, 1] and renormalized. tools/ternary_phase.py computes the
//  binodal and the spinodal of the ternary form.
//
//  P holds the pressure c cs^2, Ux/Uy the mixture velocity. Derivatives are
//  scaled with the block spacing (finest-level units); all levels use the
//  same time step, so D, nu and kappa are level-independent.
//==============================================================================
#pragma once
#include "lbm_kernels.hpp"

namespace amr {

template <int NS>
struct ParkerParams : Levels {
    real Diff, kappa, nu, omega_E;      // D, kappa, nu (lattice units), rate of the trace and a4
    real wD, wV;                        // 1/(D/cs2 + 1/2), 1/(nu/cs2 + 1/2)
    real A[NS][NS];                     // Margules A_rs (A[r][r] unused)
    real C123;                          // ternary Margules constant (NS = 3)
    real x_eps;                         // NS = 2: x1 clamped to [x_eps, 1 - x_eps]
};

namespace parker {

constexpr real cs2 = real(1) / real(3), cs4 = real(1) / real(9);
//== weights of the isotropic gradient: straight and diagonal neighbours
constexpr real g1 = real(1) / real(3), g2 = real(1) / real(12);

//== ln gamma_s of the concentrations c (Model_Parker2009_Mod: Ln_Gamma)
template <int NS>
inline void ln_gamma(const ParkerParams<NS> &p, const real (&c)[NS], real (&lg)[NS]) {
    if constexpr (NS == 2) {
        const real x1 = c[0] / (c[0] + c[1]);
        const real xa = sycl::fmin(sycl::fmax(x1, p.x_eps), real(1) - p.x_eps);
        const real xb = real(1) - xa;
        lg[0] = (p.A[0][1] + real(2) * (p.A[1][0] - p.A[0][1]) * xa) * xb * xb;
        lg[1] = (p.A[1][0] + real(2) * (p.A[0][1] - p.A[1][0]) * xb) * xa * xa;
    } else {
        real x[NS], sum = 0;
        for (int s = 0; s < NS; ++s) {
            x[s] = sycl::fmax(c[s], real(0));
            sum += x[s];
        }
        const real r = real(1) / sycl::fmax(sum, real(1e-30));
        for (int s = 0; s < NS; ++s) x[s] *= r;
        real g = 0, dg[NS] = {};
        for (int a = 0; a < NS; ++a)
            for (int b = a + 1; b < NS; ++b) {
                const real Aba = p.A[b][a], Aab = p.A[a][b];
                g += x[a] * x[b] * (Aba * x[a] + Aab * x[b]);
                dg[a] += x[b] * (real(2) * Aba * x[a] + Aab * x[b]);
                dg[b] += x[a] * (Aba * x[a] + real(2) * Aab * x[b]);
            }
        if constexpr (NS == 3) {
            g += p.C123 * x[0] * x[1] * x[2];
            dg[0] += p.C123 * x[1] * x[2];
            dg[1] += p.C123 * x[0] * x[2];
            dg[2] += p.C123 * x[0] * x[1];
        }
        for (int s = 0; s < NS; ++s) lg[s] = dg[s] - real(2) * g;
    }
}

//== ln gamma of all species on block n with its ghost layer, into local
//   memory: lg[s * SP + i + SX * j], i, j = 0..nx+1 (SX, SP: lbm_kernels.hpp)
template <int NS, class F, class Acc>
inline void stage_ln_gamma(const F &fields, const ParkerParams<NS> &p, int n, int t, int WG, Acc &lg) {
    for (int w = t; w < SP; w += WG) {
        const int i = w % SX, j = w / SX;
        real c[NS], l[NS];
        for (int s = 0; s < NS; ++s) c[s] = fields.C[s][cix(i, j, n)];
        ln_gamma<NS>(p, c, l);
        for (int s = 0; s < NS; ++s) lg[s * SP + w] = l[s];
    }
}

//== isotropic gradient of a staged field at interior cell (i, j), finest-level
//   units (Model_Parker2009_Mod: Gradient)
template <class Acc>
inline void gradient(const Acc &q, int s, int i, int j, real inv_dx, real &qx, real &qy) {
    auto Q = [&](int a, int b) { return q[s * SP + a + SX * b]; };
    qx = inv_dx * (g1 * (Q(i + 1, j) - Q(i - 1, j)) + g2 * (Q(i + 1, j + 1) + Q(i + 1, j - 1) - Q(i - 1, j + 1) - Q(i - 1, j - 1)));
    qy = inv_dx * (g1 * (Q(i, j + 1) - Q(i, j - 1)) + g2 * (Q(i + 1, j + 1) + Q(i - 1, j + 1) - Q(i + 1, j - 1) - Q(i - 1, j - 1)));
}

//== mixture velocity and pressure c cs^2 of cell (i, j) of block n
//   (Model_Parker2009_Mod: Mixture_Moments)
template <int NS, class F>
inline void mixture_moments(const F &fields, int n, int i, int j) {
    real c = 0, jx = 0, jy = 0;
    for (int s = 0; s < NS; ++s) c += fields.C[s][cix(i, j, n)];
    for (int k = 2; k <= 9; ++k) {
        real fk = 0;
        for (int s = 0; s < NS; ++s) fk += fields.f[s].in[iix(i, j, k, n)];
        jx += real(ex[k]) * fk;
        jy += real(ey[k]) * fk;
    }
    const size_t pc = pix(i, j, n);
    fields.Ux[pc] = jx / c;
    fields.Uy[pc] = jy / c;
    fields.P[pc]  = c * cs2;
}

//== central Hermite moment collision of one species at one cell
//   (Model_Parker2009_Mod: Collide_Species)
template <int NS>
inline void collide_species(real (&f)[10], real c, real ux, real uy, real px, real py, const ParkerParams<NS> &p) {
    const real wD = p.wD, wV = p.wV, kappa = p.kappa, omega_E = p.omega_E;
    //== central moments about u
    real k00 = 0, k10 = 0, k01 = 0, k20 = 0, k02 = 0, k11 = 0, k21 = 0, k12 = 0, k22 = 0;
    for (int k = 1; k <= 9; ++k) {
        const real fk = f[k];
        const real dx = real(ex[k]) - ux, dy = real(ey[k]) - uy;
        const real dx2 = dx * dx, dy2 = dy * dy;
        k00 = k00 + fk;
        k10 = k10 + fk * dx;
        k01 = k01 + fk * dy;
        k20 = k20 + fk * dx2;
        k02 = k02 + fk * dy2;
        k11 = k11 + fk * dx * dy;
        k21 = k21 + fk * dx2 * dy;
        k12 = k12 + fk * dx * dy2;
        k22 = k22 + fk * dx2 * dy2;
    }
    (void)k00;

    //== central Hermite moments, relaxed (c = k00 is conserved)
    const real h10 = k10 + wD * (c * px - k10);
    const real h01 = k01 + wD * (c * py - k01);
    const real h21 = (real(1) - wD) * (k21 - cs2 * k01);
    const real h12 = (real(1) - wD) * (k12 - cs2 * k10);
    const real h11 = k11 + wV * (kappa * c * px * py - k11);
    real tr = k20 + k02 - real(2) * cs2 * c;
    real dv = k20 - k02;
    const real trq = kappa * c * (px * px + py * py);
    const real dvq = kappa * c * (px * px - py * py);
    tr = tr + omega_E * (trq - tr);
    dv = dv + wV * (dvq - dv);
    const real h20 = real(0.5) * (tr + dv);
    const real h02 = real(0.5) * (tr - dv);
    const real h22 = (real(1) - omega_E) * (k22 - cs2 * (k20 + k02) + cs4 * c);

    //== back to central moments
    k10 = h10;
    k01 = h01;
    k20 = h20 + cs2 * c;
    k02 = h02 + cs2 * c;
    k11 = h11;
    k21 = h21 + cs2 * k01;
    k12 = h12 + cs2 * k10;
    k22 = h22 + cs2 * (k20 + k02) - cs4 * c;

    //== raw moments (binomial shift by u)
    const real m10 = k10 + ux * c;
    const real m01 = k01 + uy * c;
    const real m20 = k20 + real(2) * ux * k10 + ux * ux * c;
    const real m02 = k02 + real(2) * uy * k01 + uy * uy * c;
    const real m11 = k11 + ux * k01 + uy * k10 + ux * uy * c;
    const real m21 = k21 + real(2) * ux * k11 + uy * k20 + ux * ux * k01 + real(2) * ux * uy * k10 + ux * ux * uy * c;
    const real m12 = k12 + real(2) * uy * k11 + ux * k02 + uy * uy * k10 + real(2) * ux * uy * k01 + ux * uy * uy * c;
    const real m22 = k22 + real(2) * ux * k12 + real(2) * uy * k21 + ux * ux * k02 + uy * uy * k20 + real(4) * ux * uy * k11
                   + real(2) * ux * uy * uy * k10 + real(2) * ux * ux * uy * k01 + ux * ux * uy * uy * c;

    //== populations
    f[1] = c - m20 - m02 + m22;
    f[2] = real(0.5) * (m10 + m20 - m12 - m22);
    f[3] = real(0.5) * (m01 + m02 - m21 - m22);
    f[4] = real(0.5) * (-m10 + m20 + m12 - m22);
    f[5] = real(0.5) * (-m01 + m02 + m21 - m22);
    f[6] = real(0.25) * (m11 + m21 + m12 + m22);
    f[7] = real(0.25) * (-m11 + m21 - m12 + m22);
    f[8] = real(0.25) * (m11 - m21 - m12 + m22);
    f[9] = real(0.25) * (-m11 - m21 + m12 + m22);
}

//------------------------------------------------------------------------------
//  mixture moments (optional) and collision of all species, one cell per
//  work-item; ln gamma staged in local memory (Model_Parker2009_Mod:
//  Moments_and_Collide_Block / Collide_Block)
//------------------------------------------------------------------------------
template <int NS>
void collide(sycl::queue &q, FieldsT<NS, NS> F, Mesh M, ParkerParams<NS> p, int nLeaf, bool moments) {
    constexpr int WG = NX * NY;
    q.submit([&](sycl::handler &cgh) {
        sycl::local_accessor<real, 1> lg(sycl::range<1>(NS * SP), cgh);
        cgh.parallel_for(sycl::nd_range<1>(size_t(nLeaf) * WG, WG), [=](sycl::nd_item<1> it) {
            const int n = M.leaves[it.get_group(0)];
            const int t = int(it.get_local_id(0)), i = t % NX + 1, j = t / NX + 1;
            const real inv_dx = p.inv_dx[M.level[n - 1]];
            stage_ln_gamma<NS>(F, p, n, t, WG, lg);
            if (moments) mixture_moments<NS>(F, n, i, j);
            sycl::group_barrier(it.get_group());

            const size_t pc = pix(i, j, n);
            const real ux = F.Ux[pc], uy = F.Uy[pc];
            real c[NS], ctot = 0;
            for (int s = 0; s < NS; ++s) {
                c[s] = F.C[s][cix(i, j, n)];
                ctot += c[s];
            }
            //== diffusion velocities, mass-average reference
            real x[NS], px[NS], py[NS];
            if constexpr (NS == 2) {
                x[0] = c[0] / (c[0] + c[1]);
                x[1] = real(1) - x[0];
            } else {
                const real r = real(1) / ctot;
                for (int s = 0; s < NS; ++s) x[s] = c[s] * r;
            }
            real qx = 0, qy = 0;
            for (int s = 0; s < NS; ++s) {
                gradient(lg, s, i, j, inv_dx, px[s], py[s]);
                px[s] = -p.Diff * px[s];
                py[s] = -p.Diff * py[s];
                qx += x[s] * px[s];
                qy += x[s] * py[s];
            }
            for (int s = 0; s < NS; ++s) {
                real f[10];
                for (int k = 1; k <= 9; ++k) f[k] = F.f[s].in[iix(i, j, k, n)];
                collide_species<NS>(f, c[s], ux, uy, px[s] - qx, py[s] - qy, p);
                for (int k = 1; k <= 9; ++k) F.f[s].in[iix(i, j, k, n)] = f[k];
            }
        });
    });
}

//------------------------------------------------------------------------------
//  mixture moments of the listed blocks (Model_Parker2009_Mod: Moments_Block)
//------------------------------------------------------------------------------
template <int NS>
void moments(sycl::queue &q, FieldsT<NS, NS> F, const BlockRef *list, int count) {
    constexpr int WG = NX * NY;
    q.parallel_for(sycl::nd_range<1>(size_t(count) * WG, WG), [=](sycl::nd_item<1> it) {
        const BlockRef b = list[it.get_group(0)];
        const int t = int(it.get_local_id(0));
        mixture_moments<NS>(F, b.n, t % NX + 1, t / NX + 1);
    });
}

//------------------------------------------------------------------------------
//  refinement criterion: the largest |grad ln gamma_s| of each leaf, over its
//  cells and the species (Model_Parker2009_Mod: Activity_Error). It is zero
//  in the bulk phases and largest in the interfaces; grad mu_s is not used,
//  as it vanishes across an equilibrium interface
//------------------------------------------------------------------------------
template <int NS>
void criterion(sycl::queue &q, FieldsT<NS, NS> F, Mesh M, ParkerParams<NS> p, int nLeaf, real refine_tol,
               real deref_tol, real *err, std::uint8_t *code, int *changed) {
    constexpr int WG = NX * NY;
    q.submit([&](sycl::handler &cgh) {
        sycl::local_accessor<real, 1> lg(sycl::range<1>(NS * SP), cgh);
        cgh.parallel_for(sycl::nd_range<1>(size_t(nLeaf) * WG, WG), [=](sycl::nd_item<1> it) {
            const int m = int(it.get_group(0)), n = M.leaves[m];
            const int t = int(it.get_local_id(0)), i = t % NX + 1, j = t / NX + 1;
            const real inv_dx = p.inv_dx[M.level[n - 1]];
            stage_ln_gamma<NS>(F, p, n, t, WG, lg);
            sycl::group_barrier(it.get_group());
            real e2 = 0;
            for (int s = 0; s < NS; ++s) {
                real ax, ay;
                gradient(lg, s, i, j, inv_dx, ax, ay);
                e2 = sycl::fmax(e2, ax * ax + ay * ay);
            }
            e2 = sycl::reduce_over_group(it.get_group(), e2, sycl::maximum<real>());
            store_criterion(t, m, sycl::sqrt(e2), refine_tol, deref_tol, err, code, changed);
        });
    });
}

//------------------------------------------------------------------------------
//  initial state of the listed blocks: the mixture at rest, c_s = c_total x_s
//  with the ghost layers (Init::x: the case's mole fractions at a point,
//  finest-level lattice units), then all species at equilibrium
//  (Model_Parker2009_Mod: Initialize_Block)
//------------------------------------------------------------------------------
template <int NS, class Init>
void initial_state(sycl::queue &q, FieldsT<NS, NS> F, Init init, real c_total, int L0, const BlockRef *list,
                   int count) {
    constexpr int WG = NX * NY;
    q.parallel_for(sycl::nd_range<1>(size_t(count) * WG, WG), [=](sycl::nd_item<1> it) {
        const BlockRef b = list[it.get_group(0)];
        const int t = int(it.get_local_id(0));
        const int dx = (L0 / NX) >> b.L;
        auto xc = [&](int i) { return real(b.I * dx * NX) + real(dx) * (real(i) - real(0.5)); };
        auto yc = [&](int j) { return real(b.J * dx * NY) + real(dx) * (real(j) - real(0.5)); };
        for (int w = t; w < CX * CY; w += WG) {             // concentrations with the ghost layers
            const int i = w % CX - 1, j = w / CX - 1;
            real x[NS];
            init.x(xc(i), yc(j), x);
            for (int s = 0; s < NS; ++s) F.C[s][cix(i, j, b.n)] = c_total * x[s];
        }
        sycl::group_barrier(it.get_group());
        const int i = t % NX + 1, j = t / NX + 1;
        const size_t pc = pix(i, j, b.n);
        real c = 0;
        for (int s = 0; s < NS; ++s) {
            const real cs = F.C[s][cix(i, j, b.n)];
            for (int k = 1; k <= 9; ++k) F.f[s].in[iix(i, j, k, b.n)] = cs * Wa[k];
            c += cs;
        }
        F.Ux[pc] = 0;
        F.Uy[pc] = 0;
        F.P[pc]  = cs2 * c;
    });
}

//------------------------------------------------------------------------------
//  initial mole fractions of the cases (Parker_Cases_Mod: X1_Initial)
//------------------------------------------------------------------------------

//== a reproducible pseudo-random number in [-1, 1] of the integer position
//   (i, j) (Parker_Cases_Mod: Hash_Noise)
//   (64-bit arithmetic modulo 2^64, as the Fortran integer(8) code)
inline real hash_noise(int i, int j) {
    std::uint64_t h = std::uint64_t(std::int64_t(i) * 73856093) ^ std::uint64_t(std::int64_t(j) * 19349663);
    h = (h ^ (h >> 13)) * std::uint64_t(1274126177);
    h = (h ^ (h >> 16)) & std::uint64_t(2147483647);
    return real(2) * real(std::int64_t(h)) / real(2147483647) - real(1);
}

//== two species (namelist &parker): kind 0 droplet of radius R, elongated by
//   eps along x, rich in species 1 (xb) in the phase xa; 1 the whole domain
//   mixed, x1 = x_mix + noise; 2 the mixture in the middle L0 x L0 square
struct BinaryInit {
    int kind;
    real xa, xb, X0, Y0, R, eps, x_mix, noise, L0;
    void x(real px, real py, real (&xs)[2]) const {
        real x1;
        if (kind == 0) {
            const real a = R * (real(1) + eps), b = R / (real(1) + eps);
            const real r = sycl::sqrt(((px - X0) / a) * ((px - X0) / a) + ((py - Y0) / b) * ((py - Y0) / b));
            const real s = (r - real(1)) * sycl::sqrt(a * b);         // approximate distance to the interface
            x1 = xa + (xb - xa) * real(0.5) * (real(1) - sycl::tanh(real(0.5) * s));
        } else if (kind == 1) {
            x1 = x_mix + noise * hash_noise(int(sycl::floor(px)), int(sycl::floor(py)));
        } else {
            const real s = sycl::fmax(sycl::fabs(px - X0), sycl::fabs(py - Y0)) - (real(0.5) * L0 - real(4));
            const real a = real(0.5) * (real(1) - sycl::tanh(real(0.5) * s));
            x1 = xa + a * (x_mix + noise * hash_noise(int(sycl::floor(px)), int(sycl::floor(py))) - xa);
        }
        xs[0] = x1;
        xs[1] = real(1) - x1;
    }
};

//== three species (namelist &ouzo): the ouzo stock xo (species 1 oil, 2
//   ethanol, 3 water) and pure water xw. kind 0: a drop of radius R at
//   (X0, Y0); 1: a layer of the stock below y = Y0, its edge perturbed by
//   amp cos(2 pi x / L0); 2: the whole domain diluted, the fraction dil of
//   stock. A smooth edge of width W (tanh), and on the oil and ethanol
//   fractions a relative perturbation of amplitude noise (finest cells)
struct OuzoInit {
    int kind;
    real xo[3], xw[3];
    real X0, Y0, R, W, amp, dil, noise, L0;
    void x(real px, real py, real (&xs)[3]) const {
        real a;                                             // fraction of stock
        if (kind == 0) {
            const real r = sycl::sqrt((px - X0) * (px - X0) + (py - Y0) * (py - Y0));
            a = real(0.5) * (real(1) - sycl::tanh(real(2) * (r - R) / W));
        } else if (kind == 1) {
            const real yi = Y0 + amp * sycl::cos(real(2) * real(3.14159265358979) * px / L0);
            a = real(0.5) * (real(1) - sycl::tanh(real(2) * (py - yi) / W));
        } else {
            a = dil;
        }
        const int ix = int(sycl::floor(px)), iy = int(sycl::floor(py));
        const real n1 = real(1) + noise * hash_noise(ix, iy), n2 = real(1) + noise * hash_noise(iy + 7919, ix);
        xs[0] = (a * xo[0] + (real(1) - a) * xw[0]) * n1;
        xs[1] = (a * xo[1] + (real(1) - a) * xw[1]) * n2;
        xs[2] = real(1) - xs[0] - xs[1];
    }
};

} // namespace parker

//== the model as the solver sees it (amr_lbm_sycl.cpp: Solver)
template <int NS>
struct Parker {
    static_assert(NS == 2 || NS == 3, "two or three species");
    using Params = ParkerParams<NS>;
    using Fields = FieldsT<NS, NS>;         // f[s], C[s] = c_s = sum f[s]
    static constexpr int NG = 1;
    static constexpr real refine_tol = real(0.02), deref_tol = real(0.01);
    static constexpr int max_island = 16;   // droplet interiors: no coarse/fine edges inside a phase
    static const char *name() {
        return NS == 2 ? "two-species LBM, chemical-potential diffusion (Parker & Jovanovic 2009)"
                       : "three-species LBM, chemical-potential diffusion (Parker & Jovanovic 2009; ternary Margules)";
    }

    static void collide(sycl::queue &q, Fields F, Mesh M, const Params &p, int nLeaf, bool moments) {
        parker::collide<NS>(q, F, M, p, nLeaf, moments);
    }
    static void moments(sycl::queue &q, Fields F, const Params &, const BlockRef *list, int count) {
        parker::moments<NS>(q, F, list, count);
    }
    static void criterion(sycl::queue &q, Fields F, Mesh M, const Params &p, int nLeaf, real refine_tol,
                          real deref_tol, real *err, std::uint8_t *code, int *changed) {
        parker::criterion<NS>(q, F, M, p, nLeaf, refine_tol, deref_tol, err, code, changed);
    }
};

} // namespace amr
