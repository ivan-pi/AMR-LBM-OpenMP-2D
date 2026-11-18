//==============================================================================
//  The mass-conserving phase-field model of Fakhari, Geier & Lee (JCP 2016)
//  on the GPU (src/Model_JCP2016_Mod.f90): the phase field C in h, the
//  pressure-based hydrodynamics in g; the initial states of its cases (rising
//  bubble, Rayleigh-Taylor instability)
//==============================================================================
#pragma once
#include "lbm_kernels.hpp"

namespace amr {

struct JCPParams : Levels {
    real W, M0, wc, Gy, taul, tauh, Rhol, Rhoh, Beta, kappa;
};

namespace jcp {

using Params = JCPParams;
using Fields = FieldsT<2, 1>;       // f[0] = h, f[1] = g; C[0] = C = sum(h)

//== initial states of the cases at rest, for initial_state(): the phase
//   field C(x, y) and the pressure P(x, y, C) at a point (finest-level
//   lattice units)

// a bubble of radius R at (X0, Y0) with the Laplace pressure jump
// (Rising_Bubble_Mod: Initial_State)
struct BubbleInit {
    real R, X0, Y0, Sigma, W;
    real C(real x, real y) const {
        const real Ri = sycl::sqrt((x - X0) * (x - X0) + (y - Y0) * (y - Y0));
        return real(0.5) + real(0.5) * sycl::tanh(real(2) * (Ri - R) / W);
    }
    real P(real, real, real c) const { return -Sigma / R * c; }
};

// heavy over light fluid, interface y = Y0 + A cos(2 pi x / L0), hydrostatic
// pressure of the buoyancy (rho_h - rho) g, zero in the heavy fluid
// (Rayleigh_Taylor_Mod: Initial_State; ln(1+exp(s)) without overflow)
struct RTInit {
    real Y0, A, L0, W, dRhoG;
    real yi(real x) const { return Y0 + A * sycl::cos(real(2) * real(3.14159265358979) * x / L0); }
    real C(real x, real y) const { return real(0.5) + real(0.5) * sycl::tanh(real(2) * (y - yi(x)) / W); }
    real P(real x, real y, real) const {
        const real s = real(-4) * (y - yi(x)) / W;
        return -dRhoG * real(0.25) * W * (sycl::max(s, real(0)) + sycl::log(real(1) + sycl::exp(-sycl::fabs(s))));
    }
};

constexpr real c13 = real(1) / 3, c112 = real(1) / 12, c23 = real(2) / 3, c16 = real(1) / 6;
constexpr real w0 = real(4) / 9, ws = real(1) / 9, wd = real(1) / 36;
constexpr real grad2_min = real(1e-32);

//== chemical potential at cell (i,j) (Phase_Field_Mod: Chemical_Potential_AMR)
inline real chem_pot(const real *C, int n, int i, int j, real inv_dx, const Params &p) {
    auto Cv = [&](int a, int b) { return C[cix(a, b, n)]; };
    const real C0 = Cv(i, j);
    const real S = (Cv(i - 1, j - 1) + Cv(i + 1, j - 1) + Cv(i - 1, j + 1) + Cv(i + 1, j + 1))
                 + real(4) * (Cv(i, j - 1) + Cv(i - 1, j) + Cv(i + 1, j) + Cv(i, j + 1)) - real(20) * C0;
    return real(4) * p.Beta * C0 * (C0 - real(1)) * (C0 - real(0.5)) - p.kappa * inv_dx * inv_dx / real(6) * S;
}

//== pressure and velocity of cell (i,j) from g (Model_JCP2016_Mod: Macroscopic_Values_g_AMR)
inline void moments_cell(const real *C, const real *g, int n, int i, int j, real inv_dx, real mu,
                         const Params &p, real &Pc, real &vx, real &vy) {
    auto Cv = [&](int a, int b) { return C[cix(a, b, n)]; };
    auto G  = [&](int q) { return g[iix(i, j, q, n)]; };
    const real c1 = inv_dx / real(3), c2 = inv_dx / real(12);
    const real gx = c1 * (Cv(i + 1, j) - Cv(i - 1, j))
                  + c2 * (Cv(i + 1, j + 1) + Cv(i + 1, j - 1) - Cv(i - 1, j + 1) - Cv(i - 1, j - 1));
    const real gy = c1 * (Cv(i, j + 1) - Cv(i, j - 1))
                  + c2 * (Cv(i + 1, j + 1) + Cv(i - 1, j + 1) - Cv(i + 1, j - 1) - Cv(i - 1, j - 1));
    const real S0 = G(1) + ((G(2) + G(3)) + (G(4) + G(5))) + ((G(6) + G(7)) + (G(8) + G(9)));
    const real Sx = (G(2) - G(4)) + ((G(6) - G(7)) + (G(9) - G(8)));
    const real Sy = (G(3) - G(5)) + ((G(6) + G(7)) - (G(8) + G(9)));
    const real dRho = p.Rhoh - p.Rhol;
    const real Rho = p.Rhol + Cv(i, j) * dRho, rRho = real(1) / Rho;
    const real Fx = mu * gx, Fy = mu * gy + (p.Rhoh - Rho) * p.Gy;
    vx = (real(3) * Sx + real(0.5) * Fx) * rRho;
    vy = (real(3) * Sy + real(0.5) * Fy) * rRho;
    Pc = S0 + dRho * c16 * (vx * gx + vy * gy);
}

//------------------------------------------------------------------------------
//  refinement criterion of each leaf (Phase_Field_Mod: Max_Gradient). The
//  flagging depends on it only through err > refine_tol and err > deref_tol
//  (store_criterion)
//------------------------------------------------------------------------------
inline void criterion(sycl::queue &q, Fields F, Mesh M, Params p, int nLeaf, real refine_tol, real deref_tol,
                      real *err, std::uint8_t *code, int *changed) {
    constexpr int WG = NX * NY;
    const real *C = F.C[0];
    q.parallel_for(sycl::nd_range<1>(size_t(nLeaf) * WG, WG), [=](sycl::nd_item<1> it) {
        const int m = int(it.get_group(0)), n = M.leaves[m];
        const int t = int(it.get_local_id(0)), i = t % NX + 1, j = t / NX + 1;
        const real dx = C[cix(i + 1, j, n)] - C[cix(i - 1, j, n)];
        const real dy = C[cix(i, j + 1, n)] - C[cix(i, j - 1, n)];
        const real e2 = sycl::reduce_over_group(it.get_group(), dx * dx + dy * dy, sycl::maximum<real>());
        store_criterion(t, m, sycl::sqrt(e2) * (real(0.5) * p.inv_dx[M.level[n - 1]]), refine_tol, deref_tol, err, code,
                        changed);
    });
}

//------------------------------------------------------------------------------
//  moments of g (optional) and collision of h and g, one cell per work-item
//  (Model_JCP2016_Mod: Moments_and_Collide_Block / Collide_Block)
//------------------------------------------------------------------------------
inline void collide(sycl::queue &q, Fields F, Mesh M, Params p, int nLeaf, bool moments) {
    constexpr int WG = NX * NY;
    q.parallel_for(sycl::nd_range<1>(size_t(nLeaf) * WG, WG), [=](sycl::nd_item<1> it) {
        const int n = M.leaves[it.get_group(0)];
        const int t = int(it.get_local_id(0)), i = t % NX + 1, j = t / NX + 1;
        const real inv_dx = p.inv_dx[M.level[n - 1]];
        const real *C = F.C[0];
        auto Cv = [&](int a, int b) { return C[cix(a, b, n)]; };
        auto G  = [&](int qq) -> real & { return F.f[1].in[iix(i, j, qq, n)]; };
        auto H  = [&](int qq) -> real & { return F.f[0].in[iix(i, j, qq, n)]; };

        const real C0 = Cv(i, j);
        const real mu = chem_pot(C, n, i, j, inv_dx, p);
        const real dRho = p.Rhoh - p.Rhol, dRho3 = dRho * c13;

        real vx, vy, Pc;
        const size_t pc = pix(i, j, n);
        if (moments) {
            moments_cell(C, F.f[1].in, n, i, j, inv_dx, mu, p, Pc, vx, vy);
            F.P[pc] = Pc; F.Ux[pc] = vx; F.Uy[pc] = vy;
        } else {
            vx = F.Ux[pc]; vy = F.Uy[pc]; Pc = F.P[pc];
        }

        const real cC = real(0.5) * inv_dx, cM = real(0.25) * inv_dx;
        const real invTaul = real(1) / p.taul, dInvTau = real(1) / p.tauh - invTaul;
        const real cT = real(12) * p.M0 / p.W, omwc = real(1) - p.wc;
        const real w0d = w0 * dRho3, wsd = ws * dRho3, wdd = wd * dRho3;

        // central differences along e2, e3, e6, e7; mixed differences along all e_k
        const real eC2 = (Cv(i + 1, j) - Cv(i - 1, j)) * cC;
        const real eC3 = (Cv(i, j + 1) - Cv(i, j - 1)) * cC;
        const real eC6 = (Cv(i + 1, j + 1) - Cv(i - 1, j - 1)) * cC;
        const real eC7 = (Cv(i - 1, j + 1) - Cv(i + 1, j - 1)) * cC;
        const real C3 = real(3) * C0;
        const real eM2 = (real(5) * Cv(i + 1, j) - Cv(i - 1, j) - C3 - Cv(i + 2, j)) * cM;
        const real eM4 = (real(5) * Cv(i - 1, j) - Cv(i + 1, j) - C3 - Cv(i - 2, j)) * cM;
        const real eM3 = (real(5) * Cv(i, j + 1) - Cv(i, j - 1) - C3 - Cv(i, j + 2)) * cM;
        const real eM5 = (real(5) * Cv(i, j - 1) - Cv(i, j + 1) - C3 - Cv(i, j - 2)) * cM;
        const real eM6 = (real(5) * Cv(i + 1, j + 1) - Cv(i - 1, j - 1) - C3 - Cv(i + 2, j + 2)) * cM;
        const real eM8 = (real(5) * Cv(i - 1, j - 1) - Cv(i + 1, j + 1) - C3 - Cv(i - 2, j - 2)) * cM;
        const real eM7 = (real(5) * Cv(i - 1, j + 1) - Cv(i + 1, j - 1) - C3 - Cv(i - 2, j + 2)) * cM;
        const real eM9 = (real(5) * Cv(i + 1, j - 1) - Cv(i - 1, j + 1) - C3 - Cv(i + 2, j - 2)) * cM;

        const real DcDx = c23 * eC2 + c16 * (eC6 - eC7);
        const real DcDy = c23 * eC3 + c16 * (eC6 + eC7);
        const real DmDx = c13 * (eM2 - eM4) + c112 * ((eM6 - eM8) - (eM7 - eM9));
        const real DmDy = c13 * (eM3 - eM5) + c112 * ((eM6 - eM8) + (eM7 - eM9));

        const real g2 = DcDx * DcDx + DcDy * DcDy;
        const real rn = (g2 > grad2_min ? real(1) : real(0)) / sycl::sqrt(sycl::max(g2, grad2_min));
        const real ni = DcDx * rn, nj = DcDy * rn;

        const real base0 = real(1) - real(1.5) * (vx * vx + vy * vy);
        const real uDc = vx * DcDx + vy * DcDy, uDm = vx * DmDx + vy * DmDy;
        const real Rho = p.Rhol + C0 * dRho, rho3 = Rho * c13;
        const real Fb = (p.Rhoh - Rho) * p.Gy;
        const real amu = dRho3 + mu;
        const real cH = p.wc * C0, tCH = cT * (real(1) - C0);
        const real tt = sycl::min(sycl::max(C0, real(0)), real(1)) * dInvTau + invTaul;
        const real s9 = tt / (real(1) + real(0.5) * tt);

        {   // rest direction
            const real Gp = w0 * base0;
            H(1) = omwc * H(1) + cH * Gp;
            const real Apl = Gp * amu - w0d, fp = -vy * Gp * Fb;
            const real geqp = Pc * w0 + (Gp - w0) * rho3 + real(0.5) * (Apl * uDc - fp);
            G(1) = geqp - Apl * uDm + fp;
        }
        // a pair (qp, qm) of opposite directions: h relaxed; g - geq of both
        // and g before the relaxation of the stress moments
        auto pair = [&](int qp, int qm, real v, real w, real wdr, real sqn, real fyp, real fym,
                        real eC, real eMp, real eMm, real &dgp, real &dgm, real &gbp, real &gbm) {
            const real base = base0 + real(4.5) * v * v;
            const real Gp = w * (base + real(3) * v), Gm = w * (base - real(3) * v);
            const real sq = tCH * w * sqn;
            H(qp) = omwc * H(qp) + cH * (Gp + sq);
            H(qm) = omwc * H(qm) + cH * (Gm - sq);
            const real Apl = Gp * amu - wdr, Ami = Gm * amu - wdr;
            const real fp = fyp * Gp * Fb, fm = fym * Gm * Fb;
            const real geqp = Pc * w + (Gp - w) * rho3 - real(0.5) * (Apl * (eC - uDc) + fp);
            const real geqm = Pc * w + (Gm - w) * rho3 - real(0.5) * (Ami * (-eC - uDc) + fm);
            dgp = G(qp) - geqp;
            dgm = G(qm) - geqm;
            gbp = geqp + Apl * (eMp - uDm) + fp;
            gbm = geqm + Ami * (eMm - uDm) + fm;
        };
        real dg2, dg4, gb2, gb4, dg3, dg5, gb3, gb5, dg6, dg8, gb6, gb8, dg7, dg9, gb7, gb9;
        pair(2, 4, vx,      ws, wsd, ni,      -vy,          -vy,           eC2, eM2, eM4, dg2, dg4, gb2, gb4);
        pair(3, 5, vy,      ws, wsd, nj,      real(1) - vy, real(-1) - vy, eC3, eM3, eM5, dg3, dg5, gb3, gb5);
        pair(6, 8, vx + vy, wd, wdd, ni + nj, real(1) - vy, real(-1) - vy, eC6, eM6, eM8, dg6, dg8, gb6, gb8);
        pair(7, 9, vy - vx, wd, wdd, nj - ni, real(1) - vy, real(-1) - vy, eC7, eM7, eM9, dg7, dg9, gb7, gb9);

        const real m8 = ((dg2 + dg4) - (dg3 + dg5)) * (real(1) - s9) * real(0.25);
        const real m9 = ((dg6 + dg8) - (dg7 + dg9)) * (real(1) - s9) * real(0.25);
        G(2) = gb2 + m8; G(4) = gb4 + m8; G(3) = gb3 - m8; G(5) = gb5 - m8;
        G(6) = gb6 + m9; G(8) = gb8 + m9; G(7) = gb7 - m9; G(9) = gb9 - m9;
    });
}

//------------------------------------------------------------------------------
//  moments of g of the listed blocks (Model_JCP2016_Mod: Moments_g_Block)
//------------------------------------------------------------------------------
inline void moments(sycl::queue &q, Fields F, Params p, const BlockRef *list, int count) {
    constexpr int WG = NX * NY;
    q.parallel_for(sycl::nd_range<1>(size_t(count) * WG, WG), [=](sycl::nd_item<1> it) {
        const BlockRef b = list[it.get_group(0)];
        const int t = int(it.get_local_id(0)), i = t % NX + 1, j = t / NX + 1;
        const real inv_dx = p.inv_dx[b.L];
        const real mu = chem_pot(F.C[0], b.n, i, j, inv_dx, p);
        real Pc, vx, vy;
        moments_cell(F.C[0], F.f[1].in, b.n, i, j, inv_dx, mu, p, Pc, vx, vy);
        const size_t pc = pix(i, j, b.n);
        F.P[pc] = Pc; F.Ux[pc] = vx; F.Uy[pc] = vy;
    });
}

//------------------------------------------------------------------------------
//  initial state of the listed blocks: the case's fluids at rest (C with
//  its ghost layers, P; Init: BubbleInit, RTInit), then h and g at
//  equilibrium (Model_JCP2016_Mod: Initialize_Block). L0 / NX >> L: the cell
//  size of level L
//------------------------------------------------------------------------------
template <class Init>
void initial_state(sycl::queue &q, Fields F, Params p, Init init, int L0, const BlockRef *list, int count) {
    real *C = F.C[0];
    real *h = F.f[0].in, *g = F.f[1].in;
    constexpr int WG = NX * NY;
    q.parallel_for(sycl::nd_range<1>(size_t(count) * WG, WG), [=](sycl::nd_item<1> it) {
        const BlockRef b = list[it.get_group(0)];
        const int t = int(it.get_local_id(0));
        const int dx = (L0 / NX) >> b.L;
        auto xc = [&](int i) { return real(b.I * dx * NX) + real(dx) * (real(i) - real(0.5)); };
        auto yc = [&](int j) { return real(b.J * dx * NY) + real(dx) * (real(j) - real(0.5)); };
        for (int c = t; c < CX * CY; c += WG) {             // C with its ghost layers
            const int i = c % CX - 1, j = c / CX - 1;
            C[cix(i, j, b.n)] = init.C(xc(i), yc(j));
        }
        sycl::group_barrier(it.get_group());
        const int i = t % NX + 1, j = t / NX + 1;
        const size_t pc = pix(i, j, b.n);
        const real C0 = C[cix(i, j, b.n)];
        F.Ux[pc] = 0;
        F.Uy[pc] = 0;
        F.P[pc]  = init.P(xc(i), yc(j), C0);

        const real inv_dx = p.inv_dx[b.L];
        const real mu = chem_pot(C, b.n, i, j, inv_dx, p);
        const real dRho3 = (p.Rhoh - p.Rhol) * c13, cT = real(12) * p.M0 / p.W;
        real eDc[10];
        for (int k = 1; k <= 9; ++k)
            eDc[k] = (C[cix(i + ex[k], j + ey[k], b.n)] - C[cix(i - ex[k], j - ey[k], b.n)]) * (real(0.5) * inv_dx);
        const real DcDx = c13 * (eDc[2] - eDc[4]) + c112 * (eDc[6] - eDc[7] - eDc[8] + eDc[9]);
        const real DcDy = c13 * (eDc[3] - eDc[5]) + c112 * (eDc[6] + eDc[7] - eDc[8] - eDc[9]);
        const real g2 = DcDx * DcDx + DcDy * DcDy;
        const real rn = (g2 > grad2_min ? real(1) : real(0)) / sycl::sqrt(sycl::max(g2, grad2_min));
        const real ni = DcDx * rn, nj = DcDy * rn;
        const real U = 0, V = 0, Pc = F.P[pc];
        const real Rho = p.Rhol + C0 * (p.Rhoh - p.Rhol);
        for (int k = 1; k <= 9; ++k) {
            const real eU = real(ex[k]) * U + real(ey[k]) * V;
            const real eN = real(ex[k]) * ni + real(ey[k]) * nj;
            const real Gam = Wa[k] * (real(1) + eU * (real(3) + real(4.5) * eU) - real(1.5) * (U * U + V * V));
            h[iix(i, j, k, b.n)] = C0 * (Gam + cT * (real(1) - C0) * Wa[k] * eN);
            const real hlp = ((Gam - Wa[k]) * dRho3 + Gam * mu) * (eDc[k] - (U * DcDx + V * DcDy));
            g[iix(i, j, k, b.n)] = Pc * Wa[k] + (Gam - Wa[k]) * (Rho * c13) - real(0.5) * hlp;
        }
    });
}

} // namespace jcp

//== the model as the solver sees it (amr_lbm_sycl.cpp: Solver)
struct JCP2016 {
    using Params = JCPParams;
    using Fields = jcp::Fields;
    static constexpr int NG = 2;            // mixed differences reach two cells
    static constexpr real refine_tol = real(0.002), deref_tol = real(0.001);
    static constexpr int max_island = 0;
    static const char *name() { return "mass-conserving phase-field LBM (Fakhari, Geier & Lee, JCP 2016)"; }

    static void collide(sycl::queue &q, Fields F, Mesh M, const Params &p, int nLeaf, bool moments) {
        jcp::collide(q, F, M, p, nLeaf, moments);
    }
    static void moments(sycl::queue &q, Fields F, const Params &p, const BlockRef *list, int count) {
        jcp::moments(q, F, p, list, count);
    }
    static void criterion(sycl::queue &q, Fields F, Mesh M, const Params &p, int nLeaf, real refine_tol,
                          real deref_tol, real *err, std::uint8_t *code, int *changed) {
        jcp::criterion(q, F, M, p, nLeaf, refine_tol, deref_tol, err, code, changed);
    }
};

} // namespace amr
