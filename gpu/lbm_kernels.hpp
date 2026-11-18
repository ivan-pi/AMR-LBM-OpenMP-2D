//==============================================================================
//  Device kernels of the AMR-LBM solver (JCP 2016 model, SYCL)
//
//  The fields live on the device in the layout of the Fortran block pool
//  (src/Blocks_Mod.f90): column-major, last index = block, block indices from 1.
//      h, g      (0:nx+1, 0:ny+1, 9, cap)
//      C         (-1:nx+2, -1:ny+2, cap)
//      P, Ux, Uy (1:nx, 1:ny, cap)
//  Each kernel follows a Fortran routine (named in its comment) with the
//  same 1-based index arithmetic, so the two can be compared line by line.
//  Block-wise kernels run one work-group per block.
//  Single precision: the Arc A-series GPUs have no FP64 hardware.
//==============================================================================
#pragma once
#include <cstdint>
#include <sycl/sycl.hpp>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifndef AMR_BLOCK_SIZE
#error "AMR_BLOCK_SIZE must be defined"
#endif

namespace amr {

using real = float;

constexpr int NX = AMR_BLOCK_SIZE, NY = AMR_BLOCK_SIZE;
constexpr int HX = NX + 2, HY = NY + 2;     // h, g: 0..nx+1
constexpr int CX = NX + 4, CY = NY + 4;     // C: -1..nx+2
constexpr size_t HB = size_t(HX) * HY * 9;  // values per block
constexpr size_t CB = size_t(CX) * CY;
constexpr size_t PB = size_t(NX) * NY;
constexpr int MAXLEV = 32;

// D2Q9 (D2Q9_Mod), directions 1..9 at index 1..9
constexpr int ex[10]  = {0, 0, 1, 0, -1, 0, 1, -1, -1, 1};
constexpr int ey[10]  = {0, 0, 0, 1, 0, -1, 1, 1, -1, -1};
constexpr int opp[10] = {0, 1, 4, 5, 2, 3, 8, 9, 6, 7};
// neighbour directions k = 1..8 (Blocks_Mod): offsets and the first facing child
constexpr int di[9]      = {0, -1, 1, 0, 0, -1, -1, 1, 1};
constexpr int dj[9]      = {0, 0, 0, -1, 1, -1, 1, -1, 1};
constexpr int facing1[9] = {0, 1, 0, 2, 0, 3, 1, 2, 0};

// Fortran-order offsets (1-based block n and direction q)
inline size_t hix(int i, int j, int q, int n) {
    return size_t(i) + HX * (size_t(j) + HY * (size_t(q - 1) + 9 * size_t(n - 1)));
}
inline size_t cix(int i, int j, int n) {
    return size_t(i + 1) + CX * (size_t(j + 1) + CY * size_t(n - 1));
}
inline size_t pix(int i, int j, int n) {
    return size_t(i - 1) + NX * (size_t(j - 1) + NY * size_t(n - 1));
}

struct Params {
    real W, M0, wc, Gy, taul, tauh, Rhol, Rhoh, Beta, kappa;
    int maxLevels;
    real inv_dx[MAXLEV];
};

// block metadata (index n-1 for block n)
struct Mesh {
    int *level, *bI, *bJ, *child, *nbr, *leaf, *coarser;   // nbr: 8 per block
    int *leaves;                                            // leaf blocks, nLeaf
};

inline int nbr_of(const Mesh &m, int n, int k) { return m.nbr[8 * (n - 1) + (k - 1)]; }

//------------------------------------------------------------------------------
//  Biquadratic.inc: strip interpolations, input I2 x J2, output X2 x Y2
//  (C arrays a[i][j] = Fortran a(i+1,j+1)); first along the second index
//  with the integer weights, then along the first with the scaled ones
//------------------------------------------------------------------------------
template <int I2, int J2, int X2, int Y2>
inline void from_fine(const real (&in)[I2][J2], real (&out)[X2][Y2]) {
    real tmp[I2][Y2];
    for (int Y = 1; Y <= Y2; ++Y) {
        const int Jp = (Y + 1) % 2, j0 = 2 * Y - 1 - Jp;
        const real a = real(3 - 4 * Jp), c = real(2 - (3 - 4 * Jp));
        for (int i = 0; i < I2; ++i)
            tmp[i][Y - 1] = a * in[i][j0 - 1] + real(6) * in[i][j0] + c * in[i][j0 + 1];
    }
    const real wb = real(6) / real(64);
    for (int X = 1; X <= X2; ++X) {
        const int Jp = (X + 1) % 2, i0 = 2 * X - 1 - Jp;
        const real wa = real(3 - 4 * Jp) / real(64), wc = real(2 - (3 - 4 * Jp)) / real(64);
        for (int Y = 0; Y < Y2; ++Y)
            out[X - 1][Y] = wa * tmp[i0 - 1][Y] + wb * tmp[i0][Y] + wc * tmp[i0 + 1][Y];
    }
}

template <int I2, int J2, int X2, int Y2>
inline void from_coarse(const real (&in)[I2][J2], real (&out)[X2][Y2]) {
    real tmp[I2][Y2];
    for (int y = 1; y <= Y2; ++y) {
        const int I = (y + 1) / 2, a = 5 - 8 * ((y + 1) % 2), c = 2 - a;
        for (int i = 0; i < I2; ++i)
            tmp[i][y - 1] = real(a) * in[i][I - 1] + real(30) * in[i][I] + real(c) * in[i][I + 1];
    }
    const real wb = real(30) / real(1024);
    for (int x = 1; x <= X2; ++x) {
        const int I = (x + 1) / 2, a = 5 - 8 * ((x + 1) % 2), c = 2 - a;
        const real wa = real(a) / real(1024), wc = real(c) / real(1024);
        for (int y = 0; y < Y2; ++y)
            out[x - 1][y] = wa * tmp[I - 1][y] + wb * tmp[I][y] + wc * tmp[I + 1][y];
    }
}

//------------------------------------------------------------------------------
//  Ghost cells of C (Ghost_C_Mod): side or corner k of leaf n, NG layers.
//  second = false: copies and restrictions, true: interpolations
//------------------------------------------------------------------------------
template <int NG>
inline void fill_C(real *C, const Mesh &M, int n, int k, bool second) {
    const int nb = nbr_of(M, n, k);
    if (nb <= 0) return;
    const bool coarse = M.level[nb - 1] < M.level[n - 1];
    if (coarse != second) return;

    if (k <= 2) {                                           // Fill_Side_X
        const int s = di[k], ig = s < 0 ? 0 : NX + 1;
        real out[NG][NY];
        if (coarse) {
            const int i0 = s < 0 ? NX + 1 : 0, Y = (NY / 2) * (M.bJ[n - 1] % 2);
            real in[3][NY / 2 + 2];
            for (int l = 0; l < 3; ++l)
                for (int j = 0; j < NY / 2 + 2; ++j) in[l][j] = C[cix(i0 + s * l, Y + j, nb)];
            from_coarse<3, NY / 2 + 2, NG, NY>(in, out);
        } else {
            const int i0 = s < 0 ? NX : 1;
            if (M.leaf[nb - 1]) {
                for (int j = 1; j <= NY; ++j)
                    for (int l = 0; l < NG; ++l) C[cix(ig + s * l, j, n)] = C[cix(i0 + s * l, j, nb)];
                return;
            }
            const int c = M.child[nb - 1] + (s < 0 ? 1 : 0);     // children c, c+2
            real in[2 + NG][2 * NY];
            for (int l = 0; l < 2 + NG; ++l)
                for (int j = 1; j <= NY; ++j) {
                    in[l][j - 1]      = C[cix(i0 + s * l, j, c)];
                    in[l][NY + j - 1] = C[cix(i0 + s * l, j, c + 2)];
                }
            from_fine<2 + NG, 2 * NY, NG, NY>(in, out);
        }
        for (int j = 1; j <= NY; ++j)
            for (int l = 0; l < NG; ++l) C[cix(ig + s * l, j, n)] = out[l][j - 1];

    } else if (k <= 4) {                                    // Fill_Side_Y
        const int s = dj[k], jg = s < 0 ? 0 : NY + 1;
        real out[NX][NG];
        if (coarse) {
            const int j0 = s < 0 ? NY + 1 : 0, X = (NX / 2) * (M.bI[n - 1] % 2);
            real in[NX / 2 + 2][3];
            for (int i = 0; i < NX / 2 + 2; ++i)
                for (int l = 0; l < 3; ++l) in[i][l] = C[cix(X + i, j0 + s * l, nb)];
            from_coarse<NX / 2 + 2, 3, NX, NG>(in, out);
        } else {
            const int j0 = s < 0 ? NY : 1;
            if (M.leaf[nb - 1]) {
                for (int l = 0; l < NG; ++l)
                    for (int i = 1; i <= NX; ++i) C[cix(i, jg + s * l, n)] = C[cix(i, j0 + s * l, nb)];
                return;
            }
            const int c = M.child[nb - 1] + (s < 0 ? 2 : 0);     // children c, c+1
            real in[2 * NX][2 + NG];
            for (int l = 0; l < 2 + NG; ++l)
                for (int i = 1; i <= NX; ++i) {
                    in[i - 1][l]      = C[cix(i, j0 + s * l, c)];
                    in[NX + i - 1][l] = C[cix(i, j0 + s * l, c + 1)];
                }
            from_fine<2 * NX, 2 + NG, NX, NG>(in, out);
        }
        for (int l = 0; l < NG; ++l)
            for (int i = 1; i <= NX; ++i) C[cix(i, jg + s * l, n)] = out[i - 1][l];

    } else {                                                // Fill_Corner
        const int sx = di[k], sy = dj[k];
        const int ig = sx < 0 ? 0 : NX + 1, jg = sy < 0 ? 0 : NY + 1;
        real out[NG][NG];
        if (coarse) {
            const int px = M.bI[n - 1] % 2, py = M.bJ[n - 1] % 2;
            const int i0 = sx < 0 ? NX + 1 - (NX / 2) * px : (NX / 2) * (1 - px);
            const int j0 = sy < 0 ? NY + 1 - (NY / 2) * py : (NY / 2) * (1 - py);
            real in[3][3];
            for (int m = 0; m < 3; ++m)
                for (int l = 0; l < 3; ++l) in[l][m] = C[cix(i0 + sx * l, j0 + sy * m, nb)];
            from_coarse<3, 3, NG, NG>(in, out);
        } else {
            const int i0 = sx < 0 ? NX : 1, j0 = sy < 0 ? NY : 1;
            if (M.leaf[nb - 1]) {
                for (int m = 0; m < NG; ++m)
                    for (int l = 0; l < NG; ++l) out[l][m] = C[cix(i0 + sx * l, j0 + sy * m, nb)];
            } else {
                const int c = M.child[nb - 1] + facing1[k];
                real in[2 + NG][2 + NG];
                for (int m = 0; m < 2 + NG; ++m)
                    for (int l = 0; l < 2 + NG; ++l) in[l][m] = C[cix(i0 + sx * l, j0 + sy * m, c)];
                from_fine<2 + NG, 2 + NG, NG, NG>(in, out);
            }
        }
        for (int m = 0; m < NG; ++m)
            for (int l = 0; l < NG; ++l) C[cix(ig + sx * l, jg + sy * m, n)] = out[l][m];
    }
}

//------------------------------------------------------------------------------
//  One output value of the strip interpolations (as from_fine / from_coarse,
//  same operation order), in(i,j) the 0-based input: for the ghost fills
//  of h and g one work-item per ghost cell
//------------------------------------------------------------------------------
template <class In>
inline real fine_at(int X, int Y, In in) {
    const int JpY = (Y + 1) % 2, j0 = 2 * Y - 1 - JpY, aY = 3 - 4 * JpY, cY = 2 - aY;
    const int JpX = (X + 1) % 2, i0 = 2 * X - 1 - JpX, aX = 3 - 4 * JpX, cX = 2 - aX;
    auto tmp = [&](int i) { return real(aY) * in(i, j0 - 1) + real(6) * in(i, j0) + real(cY) * in(i, j0 + 1); };
    return real(aX) / real(64) * tmp(i0 - 1) + real(6) / real(64) * tmp(i0) + real(cX) / real(64) * tmp(i0 + 1);
}

template <class In>
inline real coarse_at(int x, int y, In in) {
    const int Iy = (y + 1) / 2, ay = 5 - 8 * ((y + 1) % 2), cy = 2 - ay;
    const int Ix = (x + 1) / 2, ax = 5 - 8 * ((x + 1) % 2), cx = 2 - ax;
    auto tmp = [&](int i) { return real(ay) * in(i, Iy - 1) + real(30) * in(i, Iy) + real(cy) * in(i, Iy + 1); };
    return real(ax) / real(1024) * tmp(Ix - 1) + real(30) / real(1024) * tmp(Ix) + real(cx) / real(1024) * tmp(Ix + 1);
}

//------------------------------------------------------------------------------
//  Ghost cells of h or g (Ghost_hg_Mod): the ghost cell at position pos
//  (1..n along a side, ignored at a corner) of side or corner k of leaf n,
//  direction q. All 8 moving directions are filled (the CPU copies only
//  those the streaming reads; the others are never read)
//------------------------------------------------------------------------------
inline void fill_f_cell(real *f, const Mesh &M, int n, int k, int pos, int q, bool second) {
    const int nb = nbr_of(M, n, k);
    if (nb <= 0) return;
    const bool coarse = M.level[nb - 1] < M.level[n - 1];
    if (coarse != second) return;
    const bool leaf = M.leaf[nb - 1];

    if (k <= 2) {                                                   // Fill_Side_X
        const int s = di[k], ig = s < 0 ? 0 : NX + 1, j = pos;
        real v;
        if (coarse) {
            const int i0 = s < 0 ? NX + 1 : 0, Y = (NY / 2) * (M.bJ[n - 1] % 2);
            v = coarse_at(1, j, [&](int l, int jj) { return f[hix(i0 + s * l, Y + jj, q, nb)]; });
        } else if (leaf) {
            v = f[hix(s < 0 ? NX : 1, j, q, nb)];
        } else {
            const int i0 = s < 0 ? NX : 1, c = M.child[nb - 1] + (s < 0 ? 1 : 0);     // children c, c+2
            v = fine_at(1, j, [&](int l, int jj) {
                return jj < NY ? f[hix(i0 + s * l, jj + 1, q, c)] : f[hix(i0 + s * l, jj - NY + 1, q, c + 2)];
            });
        }
        f[hix(ig, j, q, n)] = v;

    } else if (k <= 4) {                                            // Fill_Side_Y
        const int s = dj[k], jg = s < 0 ? 0 : NY + 1, i = pos;
        real v;
        if (coarse) {
            const int j0 = s < 0 ? NY + 1 : 0, X = (NX / 2) * (M.bI[n - 1] % 2);
            v = coarse_at(i, 1, [&](int ii, int l) { return f[hix(X + ii, j0 + s * l, q, nb)]; });
        } else if (leaf) {
            v = f[hix(i, s < 0 ? NY : 1, q, nb)];
        } else {
            const int j0 = s < 0 ? NY : 1, c = M.child[nb - 1] + (s < 0 ? 2 : 0);     // children c, c+1
            v = fine_at(i, 1, [&](int ii, int l) {
                return ii < NX ? f[hix(ii + 1, j0 + s * l, q, c)] : f[hix(ii - NX + 1, j0 + s * l, q, c + 1)];
            });
        }
        f[hix(i, jg, q, n)] = v;

    } else {                                                        // Fill_Corner
        const int sx = di[k], sy = dj[k];
        const int ig = sx < 0 ? 0 : NX + 1, jg = sy < 0 ? 0 : NY + 1;
        real v;
        if (coarse) {
            const int px = M.bI[n - 1] % 2, py = M.bJ[n - 1] % 2;
            const int i0 = sx < 0 ? NX + 1 - (NX / 2) * px : (NX / 2) * (1 - px);
            const int j0 = sy < 0 ? NY + 1 - (NY / 2) * py : (NY / 2) * (1 - py);
            v = coarse_at(1, 1, [&](int l, int m) { return f[hix(i0 + sx * l, j0 + sy * m, q, nb)]; });
        } else {
            const int i0 = sx < 0 ? NX : 1, j0 = sy < 0 ? NY : 1;
            if (leaf) {
                v = f[hix(i0, j0, q, nb)];
            } else {
                const int c = M.child[nb - 1] + facing1[k];
                v = fine_at(1, 1, [&](int l, int m) { return f[hix(i0 + sx * l, j0 + sy * m, q, c)]; });
            }
        }
        f[hix(ig, jg, q, n)] = v;
    }
}

//==============================================================================
//  Kernels
//==============================================================================
struct Fields { real *h, *g, *C, *P, *Ux, *Uy; };

// a block with its level and position (lists of blocks that are not, or not
// yet, in the uploaded leaf list)
struct BlockRef { int n, L, I, J; };

// physical parameters of the rising bubble (Rising_Bubble_Mod: Setup)
struct Bubble { real R, X0, Y0, Sigma, W; int L0; };

constexpr real c13 = real(1) / 3, c112 = real(1) / 12, c23 = real(2) / 3, c16 = real(1) / 6;
constexpr real w0 = real(4) / 9, ws = real(1) / 9, wd = real(1) / 36;
constexpr real grad2_min = real(1e-32);
constexpr real Wa[10] = {0, real(16) / 36, real(4) / 36, real(4) / 36, real(4) / 36, real(4) / 36,
                         real(1) / 36, real(1) / 36, real(1) / 36, real(1) / 36};

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
    auto G  = [&](int q) { return g[hix(i, j, q, n)]; };
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
//  Ghost cells of C (Ghost_C_Mod: Fill_Ghost_Cells_C), one pass over all
//  leaves: work-group of 16 per leaf, items 0..7 the neighbours k = 1..8,
//  then the wall rows (BC_Mod: Symmetric_C_3/4) after a group barrier
//------------------------------------------------------------------------------
template <int NG>
void ghost_C_pass(sycl::queue &q, Fields F, Mesh M, int nLeaf, bool second) {
    constexpr int WG = 16;
    real *C = F.C;
    q.parallel_for(sycl::nd_range<1>(size_t(nLeaf) * WG, WG), [=](sycl::nd_item<1> it) {
        const int n = M.leaves[it.get_group(0)];
        const int t = int(it.get_local_id(0));
        if (second && !M.coarser[n - 1]) return;            // uniform over the group
        if (t < 8) fill_C<NG>(C, M, n, t + 1, second);
        sycl::group_barrier(it.get_group());
        const bool bottom = nbr_of(M, n, 3) < 0, top = nbr_of(M, n, 4) < 0;
        if (!bottom && !top) return;
        // the whole row (first pass), or the end columns next to a side
        // interpolated from a coarser block (second pass)
        for (int X = 1 - NG + t; X <= NX + NG; X += WG) {
            if (second) {
                const int side = X < 1 ? 1 : (X > NX ? 2 : 0);
                if (side == 0) continue;
                const int nb = nbr_of(M, n, side);
                if (nb <= 0 || M.level[nb - 1] == M.level[n - 1]) continue;
            }
            for (int l = 1; l <= NG; ++l) {
                if (bottom) C[cix(X, 1 - l, n)]  = C[cix(X, l, n)];
                if (top)    C[cix(X, NY + l, n)] = C[cix(X, NY + 1 - l, n)];
            }
        }
    });
}

//------------------------------------------------------------------------------
//  Ghost cells of h and g (Ghost_hg_Mod: Fill_Ghost_Cells_hg), one pass:
//  work-group of 4*8*NX per leaf, one item per side cell and direction
//  (position fastest: contiguous along the bottom and top sides), the
//  first 32 also one corner and direction each; then the bounce-back rows
//  (BC_Mod: Bounce_Back_3/4, first pass only) after a group barrier
//------------------------------------------------------------------------------
inline void ghost_hg_pass(sycl::queue &q, Fields F, Mesh M, int nLeaf, bool second) {
    constexpr int WG = 4 * 8 * NX;
    static_assert(NX == NY && WG >= 32);
    real *h = F.h, *g = F.g;
    q.parallel_for(sycl::nd_range<1>(size_t(nLeaf) * WG, WG), [=](sycl::nd_item<1> it) {
        const int n = M.leaves[it.get_group(0)];
        const int t = int(it.get_local_id(0));
        if (second && !M.coarser[n - 1]) return;          // uniform over the group
        {
            const int pos = t % NX + 1, qd = (t / NX) % 8 + 2, k = t / (8 * NX) + 1;
            fill_f_cell(h, M, n, k, pos, qd, second);
            fill_f_cell(g, M, n, k, pos, qd, second);
        }
        if (t < 32) {
            const int qd = t % 8 + 2, k = t / 8 + 5;
            fill_f_cell(h, M, n, k, 0, qd, second);
            fill_f_cell(g, M, n, k, 0, qd, second);
        }
        if (second) return;
        sycl::group_barrier(it.get_group());
        const bool bottom = nbr_of(M, n, 3) < 0, top = nbr_of(M, n, 4) < 0;
        if (!bottom && !top) return;
        for (int w = t; w < (NX + 2) * 9; w += WG) {
            const int X = w % (NX + 2), i = w / (NX + 2) + 1;
            for (real *f : {h, g}) {
                if (bottom) {
                    const int Xn = sycl::min(sycl::max(0, X + ex[i] * ey[i]), NX + 1);
                    f[hix(X, 0, i, n)] = (i == 1) ? f[hix(X, 1, 1, n)] : f[hix(Xn, 1, opp[i], n)];
                }
                if (top) {
                    const int Xn = sycl::min(sycl::max(0, X - ex[i] * ey[i]), NX + 1);
                    f[hix(X, NY + 1, i, n)] = (i == 1) ? f[hix(X, NY, 1, n)] : f[hix(Xn, NY, opp[i], n)];
                }
            }
        }
    });
}

//------------------------------------------------------------------------------
//  refinement criterion of each leaf (Phase_Field_Mod: Max_Gradient). The
//  flagging depends on it only through err > refine_tol and err > deref_tol:
//  these two bits are kept per leaf (code) and *changed is set when any of
//  them differs from the previous call
//------------------------------------------------------------------------------
inline void criterion(sycl::queue &q, Fields F, Mesh M, Params p, int nLeaf, real refine_tol, real deref_tol,
                      real *err, std::uint8_t *code, int *changed) {
    constexpr int WG = NX * NY;
    const real *C = F.C;
    q.parallel_for(sycl::nd_range<1>(size_t(nLeaf) * WG, WG), [=](sycl::nd_item<1> it) {
        const int m = int(it.get_group(0)), n = M.leaves[m];
        const int t = int(it.get_local_id(0)), i = t % NX + 1, j = t / NX + 1;
        const real dx = C[cix(i + 1, j, n)] - C[cix(i - 1, j, n)];
        const real dy = C[cix(i, j + 1, n)] - C[cix(i, j - 1, n)];
        const real e2 = sycl::reduce_over_group(it.get_group(), dx * dx + dy * dy, sycl::maximum<real>());
        if (t == 0) {
            const real e = sycl::sqrt(e2) * (real(0.5) * p.inv_dx[M.level[n - 1]]);
            err[m] = e;
            const std::uint8_t c = std::uint8_t((e > refine_tol ? 1 : 0) | (e > deref_tol ? 2 : 0));
            if (c != code[m]) {
                code[m] = c;
                sycl::atomic_ref<int, sycl::memory_order::relaxed, sycl::memory_scope::device> a(*changed);
                a.store(1);
            }
        }
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
        const real *C = F.C;
        auto Cv = [&](int a, int b) { return C[cix(a, b, n)]; };
        auto G  = [&](int qq) -> real & { return F.g[hix(i, j, qq, n)]; };
        auto H  = [&](int qq) -> real & { return F.h[hix(i, j, qq, n)]; };

        const real C0 = Cv(i, j);
        const real mu = chem_pot(C, n, i, j, inv_dx, p);
        const real dRho = p.Rhoh - p.Rhol, dRho3 = dRho * c13;

        real vx, vy, Pc;
        const size_t pc = pix(i, j, n);
        if (moments) {
            moments_cell(C, F.g, n, i, j, inv_dx, mu, p, Pc, vx, vy);
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
//  streaming of h and g, then C = sum(h) (LBM_Common_Mod): exact shift on
//  the finest level, Lax-Wendroff with CFL number 1/dx on the coarser ones.
//  In place: each work-item reads its upwind values, then a group barrier
//------------------------------------------------------------------------------
inline void stream(sycl::queue &q, Fields F, Mesh M, Params p, int nLeaf) {
    constexpr int WG = NX * NY;
    q.parallel_for(sycl::nd_range<1>(size_t(nLeaf) * WG, WG), [=](sycl::nd_item<1> it) {
        const int n = M.leaves[it.get_group(0)];
        const int t = int(it.get_local_id(0)), X = t % NX + 1, Y = t / NX + 1;
        const int L = M.level[n - 1];
        real *h = F.h, *g = F.g;
        real hn[10], gn[10];
        if (L == p.maxLevels) {
            for (int qd = 2; qd <= 9; ++qd) {
                hn[qd] = h[hix(X - ex[qd], Y - ey[qd], qd, n)];
                gn[qd] = g[hix(X - ex[qd], Y - ey[qd], qd, n)];
            }
        } else {
            const real cfl = p.inv_dx[L];
            const real cU = real(0.5) * cfl * (real(1) + cfl), c0 = real(1) - cfl * cfl;
            const real cD = real(0.5) * cfl * (real(1) - cfl);
            for (int qd = 2; qd <= 9; ++qd) {
                hn[qd] = cU * h[hix(X - ex[qd], Y - ey[qd], qd, n)] + c0 * h[hix(X, Y, qd, n)]
                       - cD * h[hix(X + ex[qd], Y + ey[qd], qd, n)];
                gn[qd] = cU * g[hix(X - ex[qd], Y - ey[qd], qd, n)] + c0 * g[hix(X, Y, qd, n)]
                       - cD * g[hix(X + ex[qd], Y + ey[qd], qd, n)];
            }
        }
        sycl::group_barrier(it.get_group());
        real Cs = h[hix(X, Y, 1, n)];
        for (int qd = 2; qd <= 9; ++qd) {
            h[hix(X, Y, qd, n)] = hn[qd];
            g[hix(X, Y, qd, n)] = gn[qd];
            Cs += hn[qd];
        }
        F.C[cix(X, Y, n)] = Cs;
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
        const real mu = chem_pot(F.C, b.n, i, j, inv_dx, p);
        real Pc, vx, vy;
        moments_cell(F.C, F.g, b.n, i, j, inv_dx, mu, p, Pc, vx, vy);
        const size_t pc = pix(i, j, b.n);
        F.P[pc] = Pc; F.Ux[pc] = vx; F.Uy[pc] = vy;
    });
}

//------------------------------------------------------------------------------
//  initial state of the listed blocks: the bubble at rest
//  (Rising_Bubble_Mod: Initial_State), then h and g at equilibrium
//  (Model_JCP2016_Mod: Initialize_Block)
//------------------------------------------------------------------------------
inline void initial_state(sycl::queue &q, Fields F, Params p, Bubble bb, const BlockRef *list, int count) {
    constexpr int WG = NX * NY;
    q.parallel_for(sycl::nd_range<1>(size_t(count) * WG, WG), [=](sycl::nd_item<1> it) {
        const BlockRef b = list[it.get_group(0)];
        const int t = int(it.get_local_id(0));
        const int dx = (bb.L0 / NX) >> b.L;
        for (int c = t; c < CX * CY; c += WG) {             // C with its ghost layers
            const int i = c % CX - 1, j = c / CX - 1;
            const real Xc = real(b.I * dx * NX) + real(dx) * (real(i) - real(0.5));
            const real Yc = real(b.J * dx * NY) + real(dx) * (real(j) - real(0.5));
            const real Ri = sycl::sqrt((Xc - bb.X0) * (Xc - bb.X0) + (Yc - bb.Y0) * (Yc - bb.Y0));
            F.C[cix(i, j, b.n)] = real(0.5) + real(0.5) * sycl::tanh(real(2) * (Ri - bb.R) / bb.W);
        }
        sycl::group_barrier(it.get_group());
        const int i = t % NX + 1, j = t / NX + 1;
        const size_t pc = pix(i, j, b.n);
        const real C0 = F.C[cix(i, j, b.n)];
        F.Ux[pc] = 0;
        F.Uy[pc] = 0;
        F.P[pc]  = -bb.Sigma / bb.R * C0;

        const real inv_dx = p.inv_dx[b.L];
        const real mu = chem_pot(F.C, b.n, i, j, inv_dx, p);
        const real dRho3 = (p.Rhoh - p.Rhol) * c13, cT = real(12) * p.M0 / p.W;
        real eDc[10];
        for (int k = 1; k <= 9; ++k)
            eDc[k] = (F.C[cix(i + ex[k], j + ey[k], b.n)] - F.C[cix(i - ex[k], j - ey[k], b.n)]) * (real(0.5) * inv_dx);
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
            F.h[hix(i, j, k, b.n)] = C0 * (Gam + cT * (real(1) - C0) * Wa[k] * eN);
            const real hlp = ((Gam - Wa[k]) * dRho3 + Gam * mu) * (eDc[k] - (U * DcDx + V * DcDy));
            F.g[hix(i, j, k, b.n)] = Pc * Wa[k] + (Gam - Wa[k]) * (Rho * c13) - real(0.5) * hlp;
        }
    });
}

//------------------------------------------------------------------------------
//  zero all fields of the listed slots (Blocks_Mod: Clear_Slot)
//------------------------------------------------------------------------------
inline void clear_slots(sycl::queue &q, Fields F, const int *slots, int count) {
    constexpr int WG = 128;
    q.parallel_for(sycl::nd_range<1>(size_t(count) * WG, WG), [=](sycl::nd_item<1> it) {
        const size_t n = size_t(slots[it.get_group(0)]);
        const int t = int(it.get_local_id(0));
        for (size_t c = t; c < HB; c += WG) { F.h[HB * (n - 1) + c] = 0; F.g[HB * (n - 1) + c] = 0; }
        for (size_t c = t; c < CB; c += WG) F.C[CB * (n - 1) + c] = 0;
        for (size_t c = t; c < PB; c += WG) {
            F.P[PB * (n - 1) + c] = 0; F.Ux[PB * (n - 1) + c] = 0; F.Uy[PB * (n - 1) + c] = 0;
        }
    });
}

//------------------------------------------------------------------------------
//  parent -> its 4 new children c1..c1+3 (Prolong_Restrict_Mod:
//  Prolong_Bilinear): h, g bilinear from the interior, extrapolated at the
//  edges (Bilinear_from_Coarse_f); C with one ghost layer (Bilinear_from_Coarse)
//------------------------------------------------------------------------------
inline void prolong(sycl::queue &q, Fields F, const int *pairs /* parent, c1 */, int count) {
    constexpr int WG = 4 * NX * NY;
    q.parallel_for(sycl::nd_range<1>(size_t(count) * WG, WG), [=](sycl::nd_item<1> it) {
        const int m = int(it.get_group(0)), n = pairs[2 * m], c1 = pairs[2 * m + 1];
        const int t = int(it.get_local_id(0));
        {   // fine cell (fi, fj) of the 2nx x 2ny children
            const int fi = t % (2 * NX) + 1, fj = t / (2 * NX) + 1;
            const int X = sycl::max(1, sycl::min(NX - 1, fi / 2)), a = 2 * fi - 4 * X + 1, b = 4 - a;
            const int Y = sycl::max(1, sycl::min(NY - 1, fj / 2)), c = 2 * fj - 4 * Y + 1, d = 4 - c;
            const int qc = (fi - 1) / NX + 2 * ((fj - 1) / NY);       // child c1+qc
            const int ci = fi - NX * ((fi - 1) / NX), cj = fj - NY * ((fj - 1) / NY);
            for (real *f : {F.h, F.g})
                for (int k = 1; k <= 9; ++k) {
                    const real t0 = real(b) * f[hix(X, Y, k, n)]     + real(a) * f[hix(X + 1, Y, k, n)];
                    const real t1 = real(b) * f[hix(X, Y + 1, k, n)] + real(a) * f[hix(X + 1, Y + 1, k, n)];
                    f[hix(ci, cj, k, c1 + qc)] = (real(d) * t0 + real(c) * t1) / real(16);
                }
        }
        // C: fine cells 0..2nx+1 from coarse 0..nx+1, each child gets nx+2 x ny+2
        for (int w = t; w < (2 * NX + 2) * (2 * NY + 2); w += WG) {
            const int fi = w % (2 * NX + 2), fj = w / (2 * NX + 2);
            const int X = fi / 2, a = 2 * fi - 4 * X + 1, b = 4 - a;
            const int Y = fj / 2, c = 2 * fj - 4 * Y + 1, d = 4 - c;
            const real t0 = real(b) * F.C[cix(X, Y, n)]     + real(a) * F.C[cix(X + 1, Y, n)];
            const real t1 = real(b) * F.C[cix(X, Y + 1, n)] + real(a) * F.C[cix(X + 1, Y + 1, n)];
            const real v = (real(d) * t0 + real(c) * t1) / real(16);
            // the children whose 0..nx+1 range contains (fi, fj): offsets 0 and nx
            for (int qy = 0; qy < 2; ++qy)
                for (int qx = 0; qx < 2; ++qx) {
                    const int ci = fi - NX * qx, cj = fj - NY * qy;
                    if (ci < 0 || ci > NX + 1 || cj < 0 || cj > NY + 1) continue;
                    F.C[cix(ci, cj, c1 + qx + 2 * qy)] = v;
                }
        }
    });
}

//------------------------------------------------------------------------------
//  children -> parent (Prolong_Restrict_Mod: Restrict_Bilinear): the slot
//  of the parent is cleared, then each child's 2x2 means go to its quadrant
//------------------------------------------------------------------------------
inline void restrict_children(sycl::queue &q, Fields F, const int *pairs /* parent, c */, int count) {
    constexpr int WG = NX * NY;
    q.parallel_for(sycl::nd_range<1>(size_t(count) * WG, WG), [=](sycl::nd_item<1> it) {
        const int m = int(it.get_group(0)), n = pairs[2 * m], c = pairs[2 * m + 1];
        const int t = int(it.get_local_id(0));
        const size_t nn = size_t(n);
        for (size_t w = t; w < HB; w += WG) { F.h[HB * (nn - 1) + w] = 0; F.g[HB * (nn - 1) + w] = 0; }
        for (size_t w = t; w < CB; w += WG) F.C[CB * (nn - 1) + w] = 0;
        sycl::group_barrier(it.get_group());
        // parent cell (i, j) <- child q, cells (2ii+1..2ii+2, 2jj+1..2jj+2)
        const int i = t % NX + 1, j = t / NX + 1;
        const int qx = (i - 1) / (NX / 2), qy = (j - 1) / (NY / 2), ch = c + qx + 2 * qy;
        const int a = 2 * (i - 1 - qx * NX / 2) + 1, b = 2 * (j - 1 - qy * NY / 2) + 1;
        auto avg = [](real f11, real f12, real f21, real f22) {
            return real(0.5) * (real(0.5) * (f11 + f12) + real(0.5) * (f21 + f22));
        };
        for (real *f : {F.h, F.g})
            for (int k = 1; k <= 9; ++k)
                f[hix(i, j, k, n)] = avg(f[hix(a, b, k, ch)], f[hix(a, b + 1, k, ch)],
                                         f[hix(a + 1, b, k, ch)], f[hix(a + 1, b + 1, k, ch)]);
        F.C[cix(i, j, n)] = avg(F.C[cix(a, b, ch)], F.C[cix(a, b + 1, ch)],
                                F.C[cix(a + 1, b, ch)], F.C[cix(a + 1, b + 1, ch)]);
        for (real *f : {F.P, F.Ux, F.Uy})
            f[pix(i, j, n)] = avg(f[pix(a, b, ch)], f[pix(a, b + 1, ch)], f[pix(a + 1, b, ch)], f[pix(a + 1, b + 1, ch)]);
    });
}

} // namespace amr
