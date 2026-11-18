//==============================================================================
//  Device kernels of the AMR-LBM solver shared by the models (SYCL)
//
//  The fields live on the device in block pools like the Fortran ones
//  (src/Blocks_Mod.f90): column-major, last index = block, block indices from 1.
//      f[0..NF-1] distributions (h, g, ...): interior (1:nx, 1:ny, 9, cap),
//                dense planes, and apart a ghost ring per moving direction
//                (Pool, iix, gix, at)
//      C[0..NC-1] scalars with ghost layers (-1:nx+2, -1:ny+2, cap);
//                C[s] = sum of f[s] after the streaming
//      P, Ux, Uy (1:nx, 1:ny, cap)
//  The JCP 2016 model (model_jcp2016.hpp) has NF = 2 (h, g) and NC = 1 (the
//  phase field C); the multicomponent model of Parker (model_parker.hpp) one
//  distribution and one concentration per species. The kernels here (ghost
//  cells, streaming, prolongation, restriction) work on any FieldsT<NF, NC>.
//  Unlike the Fortran (0:nx+1, 0:ny+1) planes of h and g, the dense planes
//  fill whole cache lines: the halo padding doubled the memory traffic of
//  the streaming. Each kernel follows a Fortran routine (named in its
//  comment) with the same 1-based indices, so the two can be compared.
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
constexpr int CX = NX + 4, CY = NY + 4;     // C: -1..nx+2
constexpr size_t HB = size_t(NX) * NY * 9;  // a distribution: interior values per block
constexpr int NR = 2 * (NX + 2) + 2 * NY;   // ghost ring cells per direction
constexpr size_t GB = size_t(NR) * 8;       // a distribution: ghost values per block (directions 2..9)
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
constexpr real Wa[10] = {0, real(16) / 36, real(4) / 36, real(4) / 36, real(4) / 36, real(4) / 36,
                         real(1) / 36, real(1) / 36, real(1) / 36, real(1) / 36};

// a distribution: per block 9 dense nx x ny planes (an 8x8 plane of floats is
// four 64-byte cache lines; no padding) and, apart, a ring of ghost cells for
// the moving directions 2..9: rows 0 and ny+1 (x = 0..nx+1), then columns 0
// and nx+1 (y = 1..ny). Indices 1-based as in the Fortran code: interior
// 1..nx x 1..ny, ghost layer 0 and nx+1 / ny+1.
struct Pool { real *in, *gh; };
inline size_t iix(int i, int j, int q, int n) {
    return size_t(i - 1) + NX * (size_t(j - 1) + NY * (size_t(q - 1) + 9 * size_t(n - 1)));
}
inline size_t gix(int i, int j, int q, int n) {
    const int r = j == 0 ? i : j == NY + 1 ? (NX + 2) + i : i == 0 ? 2 * (NX + 2) + (j - 1) : 2 * (NX + 2) + NY + (j - 1);
    return size_t(r) + NR * (size_t(q - 2) + 8 * size_t(n - 1));
}
inline bool inside(int i, int j) { return i >= 1 && i <= NX && j >= 1 && j <= NY; }
//== a cell of block n or of its ghost layer
inline real &at(Pool f, int i, int j, int q, int n) { return inside(i, j) ? f.in[iix(i, j, q, n)] : f.gh[gix(i, j, q, n)]; }
inline size_t cix(int i, int j, int n) {
    return size_t(i + 1) + CX * (size_t(j + 1) + CY * size_t(n - 1));
}
inline size_t pix(int i, int j, int n) {
    return size_t(i - 1) + NX * (size_t(j - 1) + NY * size_t(n - 1));
}

//== the levels: the finest one and the inverse cell sizes (finest-level
//   units); the models' parameter structs derive from it
struct Levels {
    int maxLevels;
    real inv_dx[MAXLEV];
};

//== the field pool of a model with NF distributions and NC scalars
template <int NF_, int NC_>
struct FieldsT {
    static constexpr int NF = NF_, NC = NC_;
    static_assert(NC <= NF, "a scalar is the sum of the distribution of the same index");
    Pool f[NF];
    real *C[NC];
    real *P, *Ux, *Uy;
};

// a block with its level and position (lists of blocks that are not, or not
// yet, in the uploaded leaf list)
struct BlockRef { int n, L, I, J; };

// block metadata (index n-1 for block n)
struct Mesh {
    int *level, *bI, *bJ, *child, *nbr, *leaf;             // nbr: 8 per block
    int *coarser;       // a neighbour is coarser: ghost cells interpolated from it
    int *finer;         // a same-level neighbour is refined: its children read our ghost cells
    int *pullnb;        // 8 per block: the same-level leaf neighbour in direction k, else 0
    int *leaves;        // leaf blocks, nLeaf
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
//  Ghost cells of a scalar C (Ghost_C_Mod): side or corner k of leaf n, NG
//  layers. second = false: copies and restrictions, true: interpolations
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
//  of the distributions one work-item per ghost cell
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
//  Ghost cells of a distribution (Ghost_hg_Mod): the ghost cell at position
//  pos (1..n along a side, ignored at a corner) of side or corner k of leaf
//  n, direction q. All 8 moving directions are filled (the CPU copies only
//  those the streaming reads; the others are never read)
//------------------------------------------------------------------------------
inline void fill_f_cell(Pool f, const Mesh &M, int n, int k, int pos, int q, bool second) {
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
            v = coarse_at(1, j, [&](int l, int jj) { return at(f, i0 + s * l, Y + jj, q, nb); });
        } else if (leaf) {
            v = at(f, s < 0 ? NX : 1, j, q, nb);
        } else {
            const int i0 = s < 0 ? NX : 1, c = M.child[nb - 1] + (s < 0 ? 1 : 0);     // children c, c+2
            v = fine_at(1, j, [&](int l, int jj) {
                return jj < NY ? at(f, i0 + s * l, jj + 1, q, c) : at(f, i0 + s * l, jj - NY + 1, q, c + 2);
            });
        }
        at(f, ig, j, q, n) = v;

    } else if (k <= 4) {                                            // Fill_Side_Y
        const int s = dj[k], jg = s < 0 ? 0 : NY + 1, i = pos;
        real v;
        if (coarse) {
            const int j0 = s < 0 ? NY + 1 : 0, X = (NX / 2) * (M.bI[n - 1] % 2);
            v = coarse_at(i, 1, [&](int ii, int l) { return at(f, X + ii, j0 + s * l, q, nb); });
        } else if (leaf) {
            v = at(f, i, s < 0 ? NY : 1, q, nb);
        } else {
            const int j0 = s < 0 ? NY : 1, c = M.child[nb - 1] + (s < 0 ? 2 : 0);     // children c, c+1
            v = fine_at(i, 1, [&](int ii, int l) {
                return ii < NX ? at(f, ii + 1, j0 + s * l, q, c) : at(f, ii - NX + 1, j0 + s * l, q, c + 1);
            });
        }
        at(f, i, jg, q, n) = v;

    } else {                                                        // Fill_Corner
        const int sx = di[k], sy = dj[k];
        const int ig = sx < 0 ? 0 : NX + 1, jg = sy < 0 ? 0 : NY + 1;
        real v;
        if (coarse) {
            const int px = M.bI[n - 1] % 2, py = M.bJ[n - 1] % 2;
            const int i0 = sx < 0 ? NX + 1 - (NX / 2) * px : (NX / 2) * (1 - px);
            const int j0 = sy < 0 ? NY + 1 - (NY / 2) * py : (NY / 2) * (1 - py);
            v = coarse_at(1, 1, [&](int l, int m) { return at(f, i0 + sx * l, j0 + sy * m, q, nb); });
        } else {
            const int i0 = sx < 0 ? NX : 1, j0 = sy < 0 ? NY : 1;
            if (leaf) {
                v = at(f, i0, j0, q, nb);
            } else {
                const int c = M.child[nb - 1] + facing1[k];
                v = fine_at(1, 1, [&](int l, int m) { return at(f, i0 + sx * l, j0 + sy * m, q, c); });
            }
        }
        at(f, ig, jg, q, n) = v;
    }
}

//==============================================================================
//  Kernels
//==============================================================================

//------------------------------------------------------------------------------
//  Ghost cells of the scalars (Ghost_C_Mod: Fill_Ghost_Cells_C), NG layers,
//  one pass over all leaves: work-group of 16 per leaf, items 0..7 the
//  neighbours k = 1..8 of every scalar, then the wall rows (BC_Mod:
//  Symmetric_C_3/4) after a group barrier
//------------------------------------------------------------------------------
template <int NG, class F>
void ghost_C_pass(sycl::queue &q, F fields, Mesh M, const int *list, int count, bool second) {
    constexpr int WG = 16;
    q.parallel_for(sycl::nd_range<1>(size_t(count) * WG, WG), [=](sycl::nd_item<1> it) {
        const int n = list[it.get_group(0)];
        const int t = int(it.get_local_id(0));
        if (second && !M.coarser[n - 1]) return;            // uniform over the group
        if (t < 8)
            for (int s = 0; s < F::NC; ++s) fill_C<NG>(fields.C[s], M, n, t + 1, second);
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
            for (int s = 0; s < F::NC; ++s) {
                real *C = fields.C[s];
                for (int l = 1; l <= NG; ++l) {
                    if (bottom) C[cix(X, 1 - l, n)]  = C[cix(X, l, n)];
                    if (top)    C[cix(X, NY + l, n)] = C[cix(X, NY + 1 - l, n)];
                }
            }
        }
    });
}

//------------------------------------------------------------------------------
//  Ghost cells of the distributions (Ghost_hg_Mod: Fill_Ghost_Cells_hg), one
//  pass over the listed leaves: work-group of 4*8*NX per leaf, one item per
//  side cell and direction (position fastest: contiguous along the bottom
//  and top sides), the first 32 also one corner and direction each; then the
//  bounce-back rows (BC_Mod: Bounce_Back_3/4, first pass) after a group
//  barrier. The streaming reads same-level leaf neighbours directly
//  (pull()), so the first pass copies only for blocks whose ghost cells
//  finer blocks interpolate from (Mesh::finer); the list holds those and
//  the wall blocks. The second pass lists the blocks with a coarser neighbour.
//------------------------------------------------------------------------------
template <class F>
void ghost_hg_pass(sycl::queue &q, F fields, Mesh M, const int *list, int count, bool second) {
    constexpr int WG = 4 * 8 * NX;
    static_assert(NX == NY && WG >= 32);
    q.parallel_for(sycl::nd_range<1>(size_t(count) * WG, WG), [=](sycl::nd_item<1> it) {
        const int n = list[it.get_group(0)];
        const int t = int(it.get_local_id(0));
        if (second || M.finer[n - 1]) {                   // uniform over the group
            const int pos = t % NX + 1, qd = (t / NX) % 8 + 2, k = t / (8 * NX) + 1;
            for (int s = 0; s < F::NF; ++s) fill_f_cell(fields.f[s], M, n, k, pos, qd, second);
        }
        if (t < 32 && (second || M.finer[n - 1])) {
            const int qd = t % 8 + 2, k = t / 8 + 5;
            for (int s = 0; s < F::NF; ++s) fill_f_cell(fields.f[s], M, n, k, 0, qd, second);
        }
        if (second) return;
        sycl::group_barrier(it.get_group());
        const bool bottom = nbr_of(M, n, 3) < 0, top = nbr_of(M, n, 4) < 0;
        if (!bottom && !top) return;
        for (int w = t; w < (NX + 2) * 8; w += WG) {      // the rest direction is never read there
            const int X = w % (NX + 2), i = w / (NX + 2) + 2;
            for (int s = 0; s < F::NF; ++s) {
                const Pool f = fields.f[s];
                if (bottom) {
                    const int Xn = sycl::min(sycl::max(0, X + ex[i] * ey[i]), NX + 1);
                    at(f, X, 0, i, n) = at(f, Xn, 1, opp[i], n);
                }
                if (top) {
                    const int Xn = sycl::min(sycl::max(0, X - ex[i] * ey[i]), NX + 1);
                    at(f, X, NY + 1, i, n) = at(f, Xn, NY, opp[i], n);
                }
            }
        }
    });
}

//------------------------------------------------------------------------------
//  The end of a refinement criterion kernel (one work-group per leaf m): the
//  block error e of the group's leader is stored, and so are the two bits
//  the flagging depends on (e > refine_tol, e > deref_tol); *changed is set
//  when any of them differs from the previous call
//------------------------------------------------------------------------------
inline void store_criterion(int t, int m, real e, real refine_tol, real deref_tol, real *err, std::uint8_t *code,
                            int *changed) {
    if (t != 0) return;
    err[m] = e;
    const std::uint8_t c = std::uint8_t((e > refine_tol ? 1 : 0) | (e > deref_tol ? 2 : 0));
    if (c != code[m]) {
        code[m] = c;
        sycl::atomic_ref<int, sycl::memory_order::relaxed, sycl::memory_scope::device> a(*changed);
        a.store(1);
    }
}

//------------------------------------------------------------------------------
//  streaming of the distributions, then C[s] = sum(f[s]) (LBM_Common_Mod):
//  exact shift on the finest level, Lax-Wendroff with CFL number 1/dx on the
//  coarser ones. From F into Out (rest direction copied). Each work-group
//  first stages its block and a border of one cell in shared local memory,
//  every value loaded once: the border from same-level leaf neighbours
//  directly (pull), else from the ghost ring
//------------------------------------------------------------------------------
//== f at (x,y), a cell of block n or of its ghost layer: a ghost cell towards a
//   same-level leaf is read from that leaf (Ghost_hg_Mod's copy), any other
//   from the ghost layer (interpolated, restricted or a wall value)
inline real pull(Pool f, const Mesh &M, int n, int x, int y, int q) {
    const int sx = x < 1 ? -1 : (x > NX ? 1 : 0), sy = y < 1 ? -1 : (y > NY ? 1 : 0);
    if (sx != 0 || sy != 0) {
        constexpr int k_of[3][3] = {{5, 1, 6}, {3, 0, 4}, {7, 2, 8}};
        const int nb = M.pullnb[8 * (n - 1) + k_of[sx + 1][sy + 1] - 1];
        if (nb > 0) return f.in[iix(x - sx * NX, y - sy * NY, q, nb)];
    }
    return at(f, x, y, q, n);
}

//== the streaming stages its block with a border of one cell in shared local
//   memory: one (nx+2) x (ny+2) plane per distribution and moving direction
constexpr int SX = NX + 2, SP = SX * (NY + 2);

template <class F>
void stream(sycl::queue &q, F fields, F out, Mesh M, Levels p, int nLeaf) {
    constexpr int WG = NX * NY, NF = F::NF;
    q.submit([&](sycl::handler &cgh) {
        sycl::local_accessor<real, 1> sl(sycl::range<1>(NF * 8 * SP), cgh);
        cgh.parallel_for(sycl::nd_range<1>(size_t(nLeaf) * WG, WG), [=](sycl::nd_item<1> it) {
            const int n = M.leaves[it.get_group(0)];
            const int t = int(it.get_local_id(0)), X = t % NX + 1, Y = t / NX + 1;
            const int L = M.level[n - 1];
            auto S = [&](int fi, int qd, int x, int y) -> real & { return sl[(8 * fi + qd - 2) * SP + x + SX * y]; };

            //== stage: the interior (one cell per item, coalesced), then the
            //   border ring, each value once (pull: a same-level leaf
            //   neighbour, else the ghost ring)
            for (int fi = 0; fi < NF; ++fi)
                for (int qd = 2; qd <= 9; ++qd) S(fi, qd, X, Y) = fields.f[fi].in[iix(X, Y, qd, n)];
            for (int w = t; w < NF * 8 * NR; w += WG) {
                const int fi = w / (8 * NR), qd = (w / NR) % 8 + 2, r = w % NR;
                const int x = r < SX ? r : r < 2 * SX ? r - SX : r < 2 * SX + NY ? 0 : NX + 1;
                const int y = r < SX ? 0 : r < 2 * SX ? NY + 1 : r < 2 * SX + NY ? r - 2 * SX + 1 : r - 2 * SX - NY + 1;
                S(fi, qd, x, y) = pull(fields.f[fi], M, n, x, y, qd);
            }
            sycl::group_barrier(it.get_group());

            const bool finest = L == p.maxLevels;
            const real cfl = p.inv_dx[L];
            const real cU = real(0.5) * cfl * (real(1) + cfl), c0 = real(1) - cfl * cfl;
            const real cD = real(0.5) * cfl * (real(1) - cfl);
            for (int fi = 0; fi < NF; ++fi) {
                real *fo = out.f[fi].in;
                const real f1 = fields.f[fi].in[iix(X, Y, 1, n)];
                fo[iix(X, Y, 1, n)] = f1;
                real Cs = f1;
                for (int qd = 2; qd <= 9; ++qd) {
                    const real fn = finest ? S(fi, qd, X - ex[qd], Y - ey[qd])
                                           : cU * S(fi, qd, X - ex[qd], Y - ey[qd]) + c0 * S(fi, qd, X, Y)
                                                 - cD * S(fi, qd, X + ex[qd], Y + ey[qd]);
                    fo[iix(X, Y, qd, n)] = fn;
                    Cs += fn;
                }
                if (fi < F::NC) fields.C[fi][cix(X, Y, n)] = Cs;
            }
        });
    });
}

//------------------------------------------------------------------------------
//  zero all fields of the listed slots (Blocks_Mod: Clear_Slot)
//------------------------------------------------------------------------------
template <class F>
void clear_slots(sycl::queue &q, F fields, const int *slots, int count) {
    constexpr int WG = 128;
    q.parallel_for(sycl::nd_range<1>(size_t(count) * WG, WG), [=](sycl::nd_item<1> it) {
        const size_t n = size_t(slots[it.get_group(0)]);
        const int t = int(it.get_local_id(0));
        for (int s = 0; s < F::NF; ++s) {
            for (size_t c = t; c < HB; c += WG) fields.f[s].in[HB * (n - 1) + c] = 0;
            for (size_t c = t; c < GB; c += WG) fields.f[s].gh[GB * (n - 1) + c] = 0;
        }
        for (int s = 0; s < F::NC; ++s)
            for (size_t c = t; c < CB; c += WG) fields.C[s][CB * (n - 1) + c] = 0;
        for (size_t c = t; c < PB; c += WG) {
            fields.P[PB * (n - 1) + c] = 0; fields.Ux[PB * (n - 1) + c] = 0; fields.Uy[PB * (n - 1) + c] = 0;
        }
    });
}

//------------------------------------------------------------------------------
//  parent -> its 4 new children c1..c1+3 (Prolong_Restrict_Mod:
//  Prolong_Bilinear): the distributions bilinear from the interior,
//  extrapolated at the edges (Bilinear_from_Coarse_f); the scalars with one
//  ghost layer (Bilinear_from_Coarse)
//------------------------------------------------------------------------------
template <class F>
void prolong(sycl::queue &q, F fields, const int *pairs /* parent, c1 */, int count) {
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
            for (int s = 0; s < F::NF; ++s) {
                real *f = fields.f[s].in;
                for (int k = 1; k <= 9; ++k) {
                    const real t0 = real(b) * f[iix(X, Y, k, n)]     + real(a) * f[iix(X + 1, Y, k, n)];
                    const real t1 = real(b) * f[iix(X, Y + 1, k, n)] + real(a) * f[iix(X + 1, Y + 1, k, n)];
                    f[iix(ci, cj, k, c1 + qc)] = (real(d) * t0 + real(c) * t1) / real(16);
                }
            }
        }
        // scalars: fine cells 0..2nx+1 from coarse 0..nx+1, each child gets nx+2 x ny+2
        for (int w = t; w < (2 * NX + 2) * (2 * NY + 2); w += WG) {
            const int fi = w % (2 * NX + 2), fj = w / (2 * NX + 2);
            const int X = fi / 2, a = 2 * fi - 4 * X + 1, b = 4 - a;
            const int Y = fj / 2, c = 2 * fj - 4 * Y + 1, d = 4 - c;
            for (int s = 0; s < F::NC; ++s) {
                real *C = fields.C[s];
                const real t0 = real(b) * C[cix(X, Y, n)]     + real(a) * C[cix(X + 1, Y, n)];
                const real t1 = real(b) * C[cix(X, Y + 1, n)] + real(a) * C[cix(X + 1, Y + 1, n)];
                const real v = (real(d) * t0 + real(c) * t1) / real(16);
                // the children whose 0..nx+1 range contains (fi, fj): offsets 0 and nx
                for (int qy = 0; qy < 2; ++qy)
                    for (int qx = 0; qx < 2; ++qx) {
                        const int ci = fi - NX * qx, cj = fj - NY * qy;
                        if (ci < 0 || ci > NX + 1 || cj < 0 || cj > NY + 1) continue;
                        C[cix(ci, cj, c1 + qx + 2 * qy)] = v;
                    }
            }
        }
    });
}

//------------------------------------------------------------------------------
//  children -> parent (Prolong_Restrict_Mod: Restrict_Bilinear): the slot
//  of the parent is cleared, then each child's 2x2 means go to its quadrant
//------------------------------------------------------------------------------
template <class F>
void restrict_children(sycl::queue &q, F fields, const int *pairs /* parent, c */, int count) {
    constexpr int WG = NX * NY;
    q.parallel_for(sycl::nd_range<1>(size_t(count) * WG, WG), [=](sycl::nd_item<1> it) {
        const int m = int(it.get_group(0)), n = pairs[2 * m], c = pairs[2 * m + 1];
        const int t = int(it.get_local_id(0));
        const size_t nn = size_t(n);
        for (int s = 0; s < F::NF; ++s) {
            for (size_t w = t; w < HB; w += WG) fields.f[s].in[HB * (nn - 1) + w] = 0;
            for (size_t w = t; w < GB; w += WG) fields.f[s].gh[GB * (nn - 1) + w] = 0;
        }
        for (int s = 0; s < F::NC; ++s)
            for (size_t w = t; w < CB; w += WG) fields.C[s][CB * (nn - 1) + w] = 0;
        sycl::group_barrier(it.get_group());
        // parent cell (i, j) <- child q, cells (2ii+1..2ii+2, 2jj+1..2jj+2)
        const int i = t % NX + 1, j = t / NX + 1;
        const int qx = (i - 1) / (NX / 2), qy = (j - 1) / (NY / 2), ch = c + qx + 2 * qy;
        const int a = 2 * (i - 1 - qx * NX / 2) + 1, b = 2 * (j - 1 - qy * NY / 2) + 1;
        auto avg = [](real f11, real f12, real f21, real f22) {
            return real(0.5) * (real(0.5) * (f11 + f12) + real(0.5) * (f21 + f22));
        };
        for (int s = 0; s < F::NF; ++s) {
            real *f = fields.f[s].in;
            for (int k = 1; k <= 9; ++k)
                f[iix(i, j, k, n)] = avg(f[iix(a, b, k, ch)], f[iix(a, b + 1, k, ch)],
                                         f[iix(a + 1, b, k, ch)], f[iix(a + 1, b + 1, k, ch)]);
        }
        for (int s = 0; s < F::NC; ++s) {
            real *C = fields.C[s];
            C[cix(i, j, n)] = avg(C[cix(a, b, ch)], C[cix(a, b + 1, ch)], C[cix(a + 1, b, ch)], C[cix(a + 1, b + 1, ch)]);
        }
        for (real *f : {fields.P, fields.Ux, fields.Uy})
            f[pix(i, j, n)] = avg(f[pix(a, b, ch)], f[pix(a, b + 1, ch)], f[pix(a + 1, b, ch)], f[pix(a + 1, b + 1, ch)]);
    });
}

} // namespace amr
