//==============================================================================
//  AMR-LBM rising bubble on a GPU (SYCL): the JCP 2016 model of the Fortran
//  solver (src/), single precision
//
//  The host keeps the quadtree of blocks (metadata only) and runs the mesh
//  logic of Mesh_Mod / Refine_Mod / Blocks_Mod step for step: the same slot
//  allocation, Z-order leaf lists, neighbour table and flagging, so a run
//  can be compared with the Fortran program built with AMR_PRECISION=4.
//  All field work runs on the device (lbm_kernels.hpp): the time step, the
//  initial state, prolongation and restriction. Per step the host fetches
//  the refinement criterion of the leaves; after a mesh change it uploads
//  the block metadata, at diagnostic steps it fetches C and Uy.
//
//  The pieces:
//      DeviceArray<T>   a USM device allocation (RAII)
//      DeviceFields     the field pool h, g, C, P, Ux, Uy
//      DeviceMesh       the block metadata and leaf list on the device
//      Quadtree         the mesh on the host; adapt() returns a MeshChange,
//                       the field work it implies
//      BubbleCase       parameters (namelist), initial state, diagnostics
//      HostFields       a copy of C, P, Ux, Uy on the host (output steps)
//      VtkSeries        VTK XML output, the files of the Fortran VTK writer
//      StageTimers      wall time per stage of the time step
//      Solver           the time loop
//
//  Usage:  amr_lbm_sycl [number_of_time_steps [output_interval_in_steps]]
//          AMR_CASE=case.nml   namelist &rising_bubble (as the Fortran code)
//          AMR_GPU_SYNC=1      wait after every stage (stage timings)
//          AMR_OUTPUT=vtk|none field output (default vtk)
//  The host thread waits for the GPU many times per step: on a loaded
//  machine pin it to an idle CPU (taskset -c <cpu>), or every wait costs a
//  scheduler time slice. Writes bubble_rise.dat and, at the output steps,
//  vtk/Rising_Bubble_2D_NNNN.vtu, vtk/AMR_Blocks_NNNN.vtu with the series
//  Rising_Bubble_2D.pvd, AMR_Blocks.pvd (as VTK_Output_Mod).
//==============================================================================
#include "lbm_kernels.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#if !defined(AMR_L0) || !defined(AMR_MIN_LEVEL)
#error "AMR_L0 and AMR_MIN_LEVEL must be defined"
#endif

using namespace amr;

//==============================================================================
//  Configuration (Config_Mod)
//==============================================================================
constexpr int L0 = AMR_L0, minLevels = AMR_MIN_LEVEL;
constexpr int nBaseX = 1, nBaseY = 3;
constexpr int ilog2(int v) { return v <= 1 ? 0 : 1 + ilog2(v / 2); }
constexpr int maxLevels = ilog2(L0 / NX);
constexpr int MaxBlks = nBaseX * nBaseY * ((1 << (2 * (maxLevels + 1))) - 1) / 3;
static_assert((NX << maxLevels) == L0, "L0 must be a power-of-two multiple of the block size");
static_assert(maxLevels < MAXLEV, "too many levels");

inline int dx_level(int L) { return (L0 / NX) >> L; }

//==============================================================================
//  Device memory
//==============================================================================
template <class T>
class DeviceArray {
public:
    DeviceArray(sycl::queue &q, size_t n) : q_(&q), n_(n), p_(sycl::malloc_device<T>(n, q)) {
        if (!p_) throw std::runtime_error("device allocation of " + std::to_string(n * sizeof(T)) + " bytes failed");
        q.memset(p_, 0, n * sizeof(T));
    }
    ~DeviceArray() { if (p_) sycl::free(p_, *q_); }
    DeviceArray(const DeviceArray &) = delete;
    DeviceArray &operator=(const DeviceArray &) = delete;
    DeviceArray(DeviceArray &&o) noexcept : q_(o.q_), n_(o.n_), p_(std::exchange(o.p_, nullptr)) {}

    T *get() const { return p_; }
    size_t size() const { return n_; }

    //== host -> device (blocking: the source may be a temporary)
    template <class V> T *put(const V &v) {
        static_assert(sizeof(v[0]) % sizeof(T) == 0);
        const size_t bytes = v.size() * sizeof(v[0]);
        if (bytes > n_ * sizeof(T)) throw std::runtime_error("DeviceArray::put: too large");
        q_->memcpy(p_, v.data(), bytes).wait();
        return p_;
    }
    //== device -> host, the first count elements (blocking)
    void get(T *dst, size_t count) const { q_->memcpy(dst, p_, count * sizeof(T)).wait(); }

private:
    sycl::queue *q_;
    size_t n_;
    T *p_;
};

//== the field pool, slots 1..cap. h and g are double-buffered: the
//   streaming reads neighbouring blocks directly, so it writes the other
//   buffer (next()) and the two are swapped afterwards
struct DeviceFields {
    DeviceArray<real> h0, g0, h1, g1, C, P, Ux, Uy;
    int cur = 0;

    DeviceFields(sycl::queue &q, size_t cap)
        : h0(q, HB * cap), g0(q, HB * cap), h1(q, HB * cap), g1(q, HB * cap),
          C(q, CB * cap), P(q, PB * cap), Ux(q, PB * cap), Uy(q, PB * cap) {}

    Fields view() const {
        return {cur ? h1.get() : h0.get(), cur ? g1.get() : g0.get(), C.get(), P.get(), Ux.get(), Uy.get()};
    }
    Fields next() const {
        return {cur ? h0.get() : h1.get(), cur ? g0.get() : g1.get(), C.get(), P.get(), Ux.get(), Uy.get()};
    }
    void swap() { cur = 1 - cur; }
    static double mib_per_block() { return double(4 * HB + CB + 3 * PB) * sizeof(real) / (1 << 20); }
};

//==============================================================================
//  The mesh on the host (Blocks_Mod, Mesh_Mod, Refine_Mod)
//==============================================================================
struct Block {
    int level = -1, I = 0, J = 0, parent = -1, child = -1;
    int nbr[9] = {};                // nbr[k], k = 1..8
    bool leaf = false, coarser = false, finer = false, refin = false, deref = false;
};

// facing(1:2,k): children of a refined neighbour in direction k touching the block
constexpr int facing[9][2] = {{0, 0}, {1, 3}, {0, 2}, {2, 3}, {0, 1}, {3, 3}, {1, 1}, {2, 2}, {0, 0}};

//== the field work a mesh adaptation implies, in this order: restrictions
//   (parent, first child), clearing of the new slots, prolongations
//   (parent, first new child), moments of g of the new blocks
struct MeshChange {
    std::vector<int> restrict_pairs, cleared, prolong_pairs;
    std::vector<BlockRef> new_blocks;
};

class Quadtree {
public:
    Quadtree() : B_(MaxBlks + 1), at_(MaxBlks + 1, 0) {}

    const Block &operator[](int n) const { return B_[n]; }
    int nblks() const { return nblks_; }
    int unused() const { return 4 * int(free_first_.size()); }
    int levelMin() const { return levelMin_; }
    int levelMax() const { return levelMax_; }
    const std::vector<int> &leaves() const { return leaves_; }
    int nLeaf() const { return int(leaves_.size()); }

    std::vector<BlockRef> leaf_refs() const {
        std::vector<BlockRef> r;
        r.reserve(leaves_.size());
        for (int n : leaves_) r.push_back(ref(n));
        return r;
    }
    BlockRef ref(int n) const { return {n, B_[n].level, B_[n].I, B_[n].J}; }

    //== Mesh_Mod: Create_Domain (fresh slots: their fields are zero already)
    void create_domain() {
        for (int J = 0; J < nBaseY; ++J)
            for (int I = 0; I < nBaseX; ++I) {
                Block &b = B_[++nblks_];
                b = Block{};
                b.level = 0; b.I = I; b.J = J; b.leaf = true;
            }
        for (int L = 1; L <= minLevels; ++L) {
            const int nb = nblks_;
            for (int n = 1; n <= nb; ++n)
                if (B_[n].level == L - 1) add_children(n, nblks_ + 1);
        }
        leaf_lists();
    }

    //== Refine_Mod: Refine_Derefine_Mesh with the criterion err (per slot);
    //   returns the field work, empty if nothing changed
    std::optional<MeshChange> adapt(const std::vector<real> &err, real refine_tol, real deref_tol) {
        if (!flag(err, refine_tol, deref_tol)) return std::nullopt;
        MeshChange mc;

        //== Derefine: restrict the flagged groups, queue their slots
        for (int n = 1; n <= nblks_; ++n) {
            if (!B_[n].deref) continue;
            const int c = B_[n].child;
            mc.restrict_pairs.insert(mc.restrict_pairs.end(), {n, c});
            std::fill(&B_[c], &B_[c] + 4, Block{});
            B_[n].child = -1;
            B_[n].leaf = true;
            free_first_.push_back(c);
        }
        //== Refine: children of the flagged leaves, coarse levels first, in free slots
        for (int L = levelMin_; L <= levelMax_; ++L) {
            const int nb = nblks_;
            for (int n = 1; n <= nb; ++n) {
                if (B_[n].level != L || !B_[n].refin) continue;
                int c1 = nblks_ + 1;
                if (!free_first_.empty()) { c1 = free_first_.front(); free_first_.pop_front(); }
                add_children(n, c1);
                mc.prolong_pairs.insert(mc.prolong_pairs.end(), {n, c1});
                for (int c = c1; c < c1 + 4; ++c) {
                    mc.cleared.push_back(c);
                    mc.new_blocks.push_back(ref(c));
                }
            }
        }
        leaf_lists();
        return mc;
    }

private:
    std::vector<Block> B_;          // slots 1..MaxBlks
    int nblks_ = 0;
    std::deque<int> free_first_;    // first slots of free sibling groups (FIFO)
    std::vector<int> leaves_;       // by level, Z-order within a level
    int levelMin_ = 0, levelMax_ = 0;
    std::vector<int> at_;           // block at Slot(L,I,J) while rebuilding the neighbours

    static int slot(int L, int I, int J) {
        return 1 + nBaseX * nBaseY * ((1 << (2 * L)) - 1) / 3 + I + J * nBaseX * (1 << L);
    }

    //== Mesh_Mod: Add_Child_Blocks
    void add_children(int par, int c1) {
        if (c1 + 3 > MaxBlks) throw std::runtime_error("more blocks than in the fully refined quadtree");
        B_[par].child = c1;
        B_[par].leaf = false;
        for (int q = 0; q < 4; ++q) {
            Block &b = B_[c1 + q];
            b = Block{};
            b.level = B_[par].level + 1;
            b.I = 2 * B_[par].I + q % 2;
            b.J = 2 * B_[par].J + q / 2;
            b.parent = par;
            b.leaf = true;
        }
        nblks_ = std::max(nblks_, c1 + 3);
    }

    //== Mesh_Mod: Array_of_Leaf_Blocks, then the neighbour table
    void leaf_lists() {
        std::vector<int> z;
        z.reserve(nblks_);
        auto visit = [&](auto &&self, int b) -> void {        // Z-order: depth first
            if (B_[b].leaf) { z.push_back(b); return; }
            for (int c = 0; c < 4; ++c) self(self, B_[b].child + c);
        };
        for (int n = 1; n <= nBaseX * nBaseY; ++n) visit(visit, n);

        std::vector<int> start(maxLevels + 2, 0);              // stable counting sort by level
        for (int b : z) ++start[B_[b].level + 1];
        for (int L = 0; L <= maxLevels; ++L) start[L + 1] += start[L];
        levelMin_ = B_[z.front()].level;
        levelMax_ = B_[z.front()].level;
        for (int b : z) {
            levelMin_ = std::min(levelMin_, B_[b].level);
            levelMax_ = std::max(levelMax_, B_[b].level);
        }
        leaves_.assign(z.size(), 0);
        for (int b : z) leaves_[start[B_[b].level]++] = b;
        neighbours();
    }

    //== Mesh_Mod: Rebuild_Neighbor_Table (periodic in x, walls -k at the bottom
    //   and top): the same level, else the coarser block covering the position
    void neighbours() {
        auto live = [&](auto &&f) { for (int n = 1; n <= nblks_; ++n) if (B_[n].level >= 0) f(n, B_[n]); };
        live([&](int n, Block &b) { at_[slot(b.level, b.I, b.J)] = n; });
        live([&](int, Block &b) {
            const int L = b.level, w = nBaseX * (1 << L);
            b.coarser = false;
            b.finer = false;
            for (int k = 1; k <= 8; ++k) {
                const int J = b.J + dj[k];
                if (J < 0 || J >= nBaseY * (1 << L)) { b.nbr[k] = -k; continue; }
                const int I = ((b.I + di[k]) % w + w) % w;
                int nb = at_[slot(L, I, J)];
                if (nb == 0 && L > 0) nb = at_[slot(L - 1, I / 2, J / 2)];
                if (nb == 0) throw std::runtime_error("no neighbour (2:1 balance violated?)");
                b.nbr[k] = nb;
            }
        });
        live([&](int, Block &b) {
            for (int k = 1; k <= 8; ++k) {
                if (b.nbr[k] <= 0) continue;
                const Block &nb = B_[b.nbr[k]];
                if (nb.level < b.level) b.coarser = true;
                if (nb.level == b.level && !nb.leaf) b.finer = true;
            }
        });
        live([&](int, Block &b) { at_[slot(b.level, b.I, b.J)] = 0; });
    }

    //== Refine_Mod: Flag_for_refinement, Flag_for_derefinement
    bool flag(const std::vector<real> &err, real refine_tol, real deref_tol) {
        for (int n = 1; n <= nblks_; ++n)
            B_[n].refin = B_[n].leaf && B_[n].level < maxLevels && err[n] > refine_tol;
        for (int n = 1; n <= nblks_; ++n)
            if (B_[n].refin) enforce_jump(n);

        for (int n = 1; n <= nblks_; ++n) {
            Block &b = B_[n];
            b.deref = false;
            const int c = b.child;
            if (c < 0 || b.level < minLevels) continue;
            const Block *k = &B_[c];
            if (!std::all_of(k, k + 4, [](const Block &x) { return x.leaf && !x.refin; })) continue;
            b.deref = std::none_of(&err[c], &err[c] + 4, [&](real e) { return e > deref_tol; });
        }
        for (int L = levelMax_ - 1; L >= levelMin_ - 1; --L)       // finest first
            for (int n = 1; n <= nblks_; ++n)
                if (B_[n].level == L && B_[n].deref) keep_jump(n);

        return std::any_of(&B_[1], &B_[1] + nblks_, [](const Block &b) { return b.refin || b.deref; });
    }

    //== a block to be refined takes its coarser neighbours along
    void enforce_jump(int n) {
        for (int k = 1; k <= 8; ++k) {
            const int b = B_[n].nbr[k];
            if (b > 0 && !B_[b].refin && B_[b].level == B_[n].level - 1) {
                B_[b].refin = true;
                enforce_jump(b);
            }
        }
    }

    //== n stays refined if a neighbour's child touching it will be refined or stays refined
    void keep_jump(int n) {
        for (int k = 1; k <= 8; ++k) {
            const int b = B_[n].nbr[k];
            if (b <= 0 || B_[b].leaf) continue;
            for (int i = 0; i < 2; ++i) {
                const Block &c = B_[B_[b].child + facing[k][i]];
                if (c.refin || (!c.leaf && !c.deref)) { B_[n].deref = false; return; }
            }
        }
    }
};

//== the block metadata on the device, uploaded after every mesh change, and
//   the leaf lists of the ghost passes: hg1 the leaves whose h, g ghost cells
//   are copied (finer neighbour) or set by the walls, coarse those with a
//   coarser neighbour (the interpolating second passes)
class DeviceMesh {
public:
    DeviceMesh(sycl::queue &q, size_t cap)
        : level_(q, cap), I_(q, cap), J_(q, cap), child_(q, cap), nbr_(q, 8 * cap), leaf_(q, cap),
          coarser_(q, cap), finer_(q, cap), pullnb_(q, 8 * cap), leaves_(q, cap), hg1_(q, cap), coarse_(q, cap) {}

    void upload(const Quadtree &T) {
        const int nb = T.nblks();
        std::vector<int> lev(nb), I(nb), J(nb), ch(nb), nbr(8 * size_t(nb)), lf(nb), co(nb), fi(nb);
        std::vector<int> pull(8 * size_t(nb), 0);
        for (int n = 1; n <= nb; ++n) {
            const Block &b = T[n];
            lev[n - 1] = b.level; I[n - 1] = b.I; J[n - 1] = b.J; ch[n - 1] = b.child;
            lf[n - 1] = b.leaf; co[n - 1] = b.coarser; fi[n - 1] = b.finer;
            std::copy(b.nbr + 1, b.nbr + 9, &nbr[8 * size_t(n - 1)]);
            for (int k = 1; k <= 8; ++k) {
                const int x = b.nbr[k];
                if (x > 0 && T[x].leaf && T[x].level == b.level) pull[8 * size_t(n - 1) + k - 1] = x;
            }
        }
        std::vector<int> hg1, coarse;
        for (int n : T.leaves()) {
            const Block &b = T[n];
            if (b.finer || b.nbr[3] < 0 || b.nbr[4] < 0) hg1.push_back(n);
            if (b.coarser) coarse.push_back(n);
        }
        level_.put(lev); I_.put(I); J_.put(J); child_.put(ch); nbr_.put(nbr);
        leaf_.put(lf); coarser_.put(co); finer_.put(fi); pullnb_.put(pull); leaves_.put(T.leaves());
        hg1_.put(hg1); coarse_.put(coarse);
        n_hg1_ = int(hg1.size());
        n_coarse_ = int(coarse.size());
    }

    Mesh view() const {
        return {level_.get(), I_.get(), J_.get(), child_.get(), nbr_.get(), leaf_.get(),
                coarser_.get(), finer_.get(), pullnb_.get(), leaves_.get()};
    }
    const int *leaves() const { return leaves_.get(); }
    const int *hg1() const { return hg1_.get(); }
    const int *coarse() const { return coarse_.get(); }
    int n_hg1() const { return n_hg1_; }
    int n_coarse() const { return n_coarse_; }

private:
    DeviceArray<int> level_, I_, J_, child_, nbr_, leaf_, coarser_, finer_, pullnb_, leaves_, hg1_, coarse_;
    int n_hg1_ = 0, n_coarse_ = 0;
};

//==============================================================================
//  Rising bubble case (Rising_Bubble_Mod), in single precision as the
//  Fortran program built with AMR_PRECISION=4
//==============================================================================
class BubbleCase {
public:
    real t_end = 8, time_per_step = 0;

    //== Rising_Bubble_Mod: Setup; the physical parameters for the kernels
    Params setup() {
        R_ = real(L0) / real(10);
        D_ = real(2) * R_;
        read_namelist();

        Params p{};
        p.Rhoh = 1;
        p.Rhol = p.Rhoh / rho_ratio_;
        p.Gy = Gy_lu_;
        p.W = Cn_ > 0 ? Cn_ * real(L0) : Width_;
        p.M0 = Pe_ > 0 ? std::sqrt(p.Gy * D_) * p.W / Pe_ : Mobility_;
        p.wc = real(1) / (real(0.5) + p.M0 * real(3));
        Sigma_ = p.Gy * (p.Rhoh - p.Rhol) * D_ * D_ / Eo_;
        mu_h_ = std::sqrt(p.Gy * p.Rhoh * (p.Rhoh - p.Rhol) * D_ * D_ * D_) / Re_;
        const real mul = mu_h_ / mu_ratio_;
        p.tauh = real(3) * (mu_h_ / p.Rhoh);
        p.taul = real(3) * (mul / p.Rhol);
        p.Beta = real(12) * Sigma_ / p.W;
        p.kappa = real(1.5) * Sigma_ * p.W;
        p.maxLevels = maxLevels;
        for (int L = 0; L <= maxLevels; ++L) p.inv_dx[L] = real(NX) / real(L0) * std::ldexp(real(1), L);
        time_per_step = std::sqrt(p.Gy / D_);
        std::printf("\nEo, Re_h, muh, mul, Pe =%12.3f%12.3f%12.3f%12.3f%12.3f\n", Eo_, Re_, mu_h_, mul,
                    std::sqrt(p.Gy * D_) * p.W / p.M0);
        params_ = p;
        return p;
    }

    real sigma() const { return Sigma_; }

    //== the bubble for the initial-state kernel
    Bubble bubble() const { return {R_, real(L0 / 2), real(L0 / 2), Sigma_, params_.W, L0}; }

    //== one line of diagnostics (Solver_Mod: Print_Diagnostics) from C and Uy
    void print_line(real time, const Quadtree &T, const std::vector<real> &C, const std::vector<real> &Uy) const {
        const Params &p = params_;
        real rmin = 10, rmax = -10, umax = 0;
        double mass = 0;
        for_leaf_cells(T, [&](int n, int i, int j, int dx) {
            const real rho = p.Rhol + C[cix(i, j, n)] * (p.Rhoh - p.Rhol);
            rmin = std::min(rmin, rho);
            rmax = std::max(rmax, rho);
            mass += double(rho) * dx * dx;
            umax = std::max(umax, std::fabs(Uy[pix(i, j, n)]));
        });
        std::printf("%7.2f%7d%7d%6d%11.5f%11.5f%9.5f%12.1f%4d%4d\n", time, T.nblks(), T.nLeaf(), T.unused(),
                    rmin, rmax, umax, mass, T.levelMin(), T.levelMax());
    }

    //== Rising_Bubble_Mod: Bubble_Diagnostics, appended to bubble_rise.dat
    void rise(real time, const Quadtree &T, const std::vector<real> &C, const std::vector<real> &Uy) {
        if (!rise_) {
            rise_.reset(std::fopen("bubble_rise.dat", "w"));
            std::fprintf(rise_.get(), "#      t*           y_c          u_r          Re_r\n");
        }
        double A = 0, Ay = 0, Av = 0;
        for_leaf_cells(T, [&](int n, int i, int j, int dx) {
            const real c = C[cix(i, j, n)];
            if (c >= real(0.9)) return;
            const double w = double(real(1) - std::max(c, real(0))) * dx * dx;
            A += w;
            Ay += w * (double(T[n].J) * dx * NY + dx * (j - 0.5));
            Av += w * Uy[pix(i, j, n)];
        });
        std::fprintf(rise_.get(), "%13.5E%13.5E%13.5E%13.5E\n", double(time), Ay / A, Av / A,
                     double(params_.Rhoh) * (Av / A) * D_ / mu_h_);
        std::fflush(rise_.get());
    }

private:
    real Eo_ = 40, Re_ = 40, rho_ratio_ = 1000, mu_ratio_ = 100, Gy_lu_ = real(5e-6);
    real Mobility_ = real(0.03), Pe_ = 0, Width_ = 5, Cn_ = 0;
    real R_ = 0, D_ = 0, mu_h_ = 0, Sigma_ = 0;
    Params params_{};
    std::unique_ptr<FILE, int (*)(FILE *)> rise_{nullptr, &std::fclose};

    template <class F> static void for_leaf_cells(const Quadtree &T, F &&f) {
        for (int n : T.leaves()) {
            const int dx = dx_level(T[n].level);
            for (int j = 1; j <= NY; ++j)
                for (int i = 1; i <= NX; ++i) f(n, i, j, dx);
        }
    }

    //== &rising_bubble ... / from the file named by AMR_CASE: "name = value"
    //   items, comments (!) and Fortran d exponents allowed
    void read_namelist() {
        const char *file = std::getenv("AMR_CASE");
        if (!file || !*file) return;
        std::ifstream in(file);
        if (!in) throw std::runtime_error("AMR_CASE: cannot open the case file");
        const std::map<std::string, real *> vars = {
            {"eo", &Eo_}, {"re", &Re_}, {"rho_ratio", &rho_ratio_}, {"mu_ratio", &mu_ratio_},
            {"gy_lu", &Gy_lu_}, {"mobility", &Mobility_}, {"pe", &Pe_}, {"width", &Width_},
            {"cn", &Cn_}, {"t_end", &t_end}};
        auto trim = [](std::string s) {
            s.erase(0, s.find_first_not_of(" \t\r"));
            s.erase(s.find_last_not_of(" \t\r") + 1);
            return s;
        };
        std::string body;           // the group, comments removed, lower case, items separated by ','
        bool inside = false;
        for (std::string line; std::getline(in, line);) {
            if (auto p = line.find('!'); p != std::string::npos) line.erase(p);
            std::transform(line.begin(), line.end(), line.begin(), [](unsigned char c) { return char(std::tolower(c)); });
            if (!inside) {
                const auto p = line.find("&rising_bubble");
                if (p == std::string::npos) continue;
                inside = true;
                line.erase(0, p + 14);
            }
            if (const auto p = line.find('/'); p != std::string::npos) { body += line.substr(0, p); break; }
            body += line + ",";
        }
        std::istringstream items(body);
        for (std::string item; std::getline(items, item, ',');) {
            const auto eq = item.find('=');
            if (eq == std::string::npos) continue;
            const std::string key = trim(item.substr(0, eq));
            std::string val = trim(item.substr(eq + 1));
            std::replace(val.begin(), val.end(), 'd', 'e');
            const auto v = vars.find(key);
            if (v == vars.end()) throw std::runtime_error("AMR_CASE: unknown name " + key);
            *v->second = real(std::strtod(val.c_str(), nullptr));
        }
        std::printf("\nCase file: %s\n", file);
    }
};

//==============================================================================
//  Host copy of the fields of blocks 1..nblks needed for output: C and Uy at
//  every diagnostic step, P and Ux as well at the output steps
//==============================================================================
struct HostFields {
    std::vector<real> C, P, Ux, Uy;

    void fetch(const DeviceFields &d, int nblks, bool all) {
        C.resize(CB * nblks);
        Uy.resize(PB * nblks);
        d.C.get(C.data(), C.size());
        d.Uy.get(Uy.data(), Uy.size());
        if (!all) return;
        P.resize(PB * nblks);
        Ux.resize(PB * nblks);
        d.P.get(P.data(), P.size());
        d.Ux.get(Ux.data(), Ux.size());
    }
};

//==============================================================================
//  VTK XML output (VTK_Output_Mod: Write_VTU_Step): per output step the leaf
//  cells and the leaf blocks as UnstructuredGrids of quads in vtk/, raw
//  binary appended data with UInt64 size headers, and the two .pvd time
//  series rewritten every step (valid if the run stops). Coordinates /L0.
//==============================================================================
class VtkSeries {
public:
    void step(real time, const Quadtree &T, const HostFields &f) {
        if (times_.empty()) std::filesystem::create_directories(dir_);
        times_.push_back(time);
        char id[8];
        std::snprintf(id, sizeof id, "%04d", int(times_.size()) - 1);
        cells(T, f).write(dir_ + "/Rising_Bubble_2D_" + id + ".vtu");
        blocks(T).write(dir_ + "/AMR_Blocks_" + id + ".vtu");
        pvd("Rising_Bubble_2D");
        pvd("AMR_Blocks");
    }

private:
    const std::string dir_ = "vtk";
    std::vector<real> times_;

    //== quads with scalar real fields, an optional 3-vector and integer fields
    struct UGrid {
        std::vector<real> pts;              // 3 per point
        std::vector<std::int32_t> conn, offs;
        std::vector<std::uint8_t> types;
        std::vector<std::pair<std::string, std::vector<real>>> rfields;
        std::string vname;
        std::vector<real> vdata;            // 3 per cell
        std::vector<std::pair<std::string, std::vector<std::int32_t>>> ifields;

        void add_quad(std::int32_t p0, std::int32_t p1, std::int32_t p2, std::int32_t p3) {
            conn.insert(conn.end(), {p0, p1, p2, p3});
            offs.push_back(std::int32_t(conn.size()));
            types.push_back(9);             // VTK_QUAD
        }

        void write(const std::string &file) const {
            static_assert(sizeof(real) == 4, "Float32 output");
            std::ofstream out(file, std::ios::binary);
            if (!out) throw std::runtime_error("cannot write " + file);
            const size_t nCells = types.size(), nPts = pts.size() / 3;
            std::uint64_t off = 0;
            auto header = [&](const std::string &name, const char *type, int ncomp, std::uint64_t bytes) {
                out << "  <DataArray type=\"" << type << "\" Name=\"" << name << "\" NumberOfComponents=\""
                    << ncomp << "\" format=\"appended\" offset=\"" << off << "\"/>\n";
                off += 8 + bytes;
            };
            out << "<?xml version=\"1.0\"?>\n"
                << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
                << "<UnstructuredGrid>\n"
                << "<Piece NumberOfPoints=\"" << nPts << "\" NumberOfCells=\"" << nCells << "\">\n"
                << "<CellData>\n";
            for (const auto &[name, v] : rfields) header(name, "Float32", 1, v.size() * 4);
            if (!vname.empty()) header(vname, "Float32", 3, vdata.size() * 4);
            for (const auto &[name, v] : ifields) header(name, "Int32", 1, v.size() * 4);
            out << "</CellData>\n<Points>\n";
            header("Points", "Float32", 3, pts.size() * 4);
            out << "</Points>\n<Cells>\n";
            header("connectivity", "Int32", 1, conn.size() * 4);
            header("offsets", "Int32", 1, offs.size() * 4);
            header("types", "UInt8", 1, types.size());
            out << "</Cells>\n</Piece>\n</UnstructuredGrid>\n<AppendedData encoding=\"raw\">\n_";
            auto blob = [&](const auto &v) {
                const std::uint64_t n = v.size() * sizeof(v[0]);
                out.write(reinterpret_cast<const char *>(&n), 8);
                out.write(reinterpret_cast<const char *>(v.data()), std::streamsize(n));
            };
            for (const auto &[name, v] : rfields) blob(v);
            if (!vname.empty()) blob(vdata);
            for (const auto &[name, v] : ifields) blob(v);
            blob(pts); blob(conn); blob(offs); blob(types);
            out << "\n</AppendedData>\n</VTKFile>\n";
        }
    };

    //== all leaf cells, with C, P, Velocity, Level and Block
    static UGrid cells(const Quadtree &T, const HostFields &f) {
        UGrid g;
        const size_t nCells = size_t(T.nLeaf()) * NX * NY;
        std::vector<real> C, P;
        std::vector<std::int32_t> level, block;
        C.reserve(nCells); P.reserve(nCells); level.reserve(nCells); block.reserve(nCells);
        g.vname = "Velocity";
        std::int32_t p0 = 0;
        for (int b : T.leaves()) {
            const Block &blk = T[b];
            const int dx = dx_level(blk.level);
            for (int j = 1; j <= NY + 1; ++j)
                for (int i = 1; i <= NX + 1; ++i)
                    g.pts.insert(g.pts.end(), {real(blk.I * dx * NX + dx * (i - 1)) / real(L0),
                                               real(blk.J * dx * NY + dx * (j - 1)) / real(L0), real(0)});
            for (int j = 1; j <= NY; ++j)
                for (int i = 1; i <= NX; ++i) {
                    const std::int32_t k = p0 + (j - 1) * (NX + 1) + i - 1;    // lower-left point
                    g.add_quad(k, k + 1, k + 1 + NX + 1, k + NX + 1);
                    C.push_back(f.C[cix(i, j, b)]);
                    P.push_back(f.P[pix(i, j, b)]);
                    g.vdata.insert(g.vdata.end(), {f.Ux[pix(i, j, b)], f.Uy[pix(i, j, b)], real(0)});
                    level.push_back(blk.level);
                    block.push_back(b);
                }
            p0 += (NX + 1) * (NY + 1);
        }
        g.rfields = {{"C", std::move(C)}, {"P", std::move(P)}};
        g.ifields = {{"Level", std::move(level)}, {"Block", std::move(block)}};
        return g;
    }

    //== one quad per leaf block, with dx, Level, Block and Parent
    static UGrid blocks(const Quadtree &T) {
        UGrid g;
        std::vector<real> dxs;
        std::vector<std::int32_t> level, block, parent;
        std::int32_t k = 0;
        for (int b : T.leaves()) {
            const Block &blk = T[b];
            const int dx = dx_level(blk.level);
            const real x0 = real(blk.I * dx * NX) / real(L0), x1 = real((blk.I + 1) * dx * NX) / real(L0);
            const real y0 = real(blk.J * dx * NY) / real(L0), y1 = real((blk.J + 1) * dx * NY) / real(L0);
            g.pts.insert(g.pts.end(), {x0, y0, 0, x1, y0, 0, x1, y1, 0, x0, y1, 0});
            g.add_quad(k, k + 1, k + 2, k + 3);
            k += 4;
            dxs.push_back(real(dx) / real(L0));
            level.push_back(blk.level);
            block.push_back(b);
            parent.push_back(blk.parent);
        }
        g.rfields = {{"dx", std::move(dxs)}};
        g.ifields = {{"Level", std::move(level)}, {"Block", std::move(block)}, {"Parent", std::move(parent)}};
        return g;
    }

    //== the time series of the files written so far (timestep as Fortran ES16.8E3)
    void pvd(const std::string &base) const {
        std::ofstream out(base + ".pvd");
        out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\">\n<Collection>\n";
        for (size_t k = 0; k < times_.size(); ++k) {
            char t[32], id[8];
            std::snprintf(t, sizeof t, "%.8E", double(times_[k]));
            std::string ts(t);                                  // exponent with 3 digits
            const auto e = ts.find('E');
            const std::string mant = ts.substr(0, e + 2), ex = ts.substr(e + 2);
            ts = mant + std::string(3 - std::min<size_t>(3, ex.size()), '0') + ex;
            std::snprintf(id, sizeof id, "%04zu", k);
            out << "  <DataSet timestep=\"" << ts << "\" part=\"0\" file=\"" << dir_ << "/" << base << "_" << id
                << ".vtu\"/>\n";
        }
        out << "</Collection>\n</VTKFile>\n";
    }
};

//==============================================================================
//  Wall time per stage; with sync the queue is drained before each reading
//==============================================================================
class StageTimers {
public:
    using clock = std::chrono::steady_clock;
    StageTimers(sycl::queue &q, bool sync) : q_(q), sync_(sync), last_(clock::now()) {}

    template <class F> void operator()(const std::string &stage, F &&work) {
        work();
        if (sync_) q_.wait();
        const auto now = clock::now();
        acc_[stage] += std::chrono::duration<double>(now - last_).count();
        last_ = now;
    }

    void report() const {
        std::printf("\n------------------------------\n");
        for (const auto &[stage, t] : acc_) std::printf(" t_%-15s%12.3f\n", (stage + ":").c_str(), t);
        if (!sync_) std::printf(" (the stages overlap on the device: set AMR_GPU_SYNC=1 for stage times)\n");
    }

private:
    sycl::queue &q_;
    bool sync_;
    clock::time_point last_;
    std::map<std::string, double> acc_;
};

//==============================================================================
//  The solver
//==============================================================================
class Solver {
public:
    Solver()
        : q_(sycl::gpu_selector_v, sycl::property::queue::in_order{}),
          fields_(q_, MaxBlks), mesh_(q_, MaxBlks),
          err_(q_, MaxBlks), code_(q_, MaxBlks), changed_(q_, 1),
          list_(q_, 4 * size_t(MaxBlks)), list2_(q_, 4 * size_t(MaxBlks)),
          err_slot_(MaxBlks + 1, 0) {}

    void run(int argc, char **argv);

private:
    sycl::queue q_;
    DeviceFields fields_;
    DeviceMesh mesh_;
    DeviceArray<real> err_;
    DeviceArray<std::uint8_t> code_;    // criterion bits per leaf (see criterion())
    DeviceArray<int> changed_;
    bool reflag_ = true;                // the bits are stale: flag on the host anyway
    DeviceArray<int> list_, list2_;     // block lists for the kernels
    Quadtree tree_;
    BubbleCase case_;
    Params p_{};
    std::vector<real> err_leaf_, err_slot_;
    HostFields host_;
    std::optional<VtkSeries> vtk_;     // empty: AMR_OUTPUT=none

    static constexpr int ghost_C = 2;   // JCP 2016: mixed differences reach two cells
    static constexpr real refine_tol = real(0.002), deref_tol = real(0.001);

    Fields F() const { return fields_.view(); }
    Mesh M() const { return mesh_.view(); }
    int nLeaf() const { return tree_.nLeaf(); }
    const BlockRef *put_refs(const std::vector<BlockRef> &r) {
        return reinterpret_cast<const BlockRef *>(list_.put(r));
    }

    //== Solver_Mod: Initialize_State
    void initialize_state() {
        const auto r = tree_.leaf_refs();
        initial_state(q_, F(), p_, case_.bubble(), put_refs(r), int(r.size()));
    }

    void ghost_cells_C() {
        ghost_C_pass<ghost_C>(q_, F(), M(), mesh_.leaves(), nLeaf(), false);
        if (mesh_.n_coarse() > 0) ghost_C_pass<ghost_C>(q_, F(), M(), mesh_.coarse(), mesh_.n_coarse(), true);
    }

    void ghost_cells_hg() {
        if (mesh_.n_hg1() > 0) ghost_hg_pass(q_, F(), M(), mesh_.hg1(), mesh_.n_hg1(), false);
        if (mesh_.n_coarse() > 0) ghost_hg_pass(q_, F(), M(), mesh_.coarse(), mesh_.n_coarse(), true);
    }

    //== criterion on the device, flagging on the host, field work on the device.
    //   The flagging is a function of the criterion bits of the leaves: if
    //   none changed since a flagging without mesh change, it is skipped
    bool adapt_mesh() {
        q_.fill(changed_.get(), 0, 1);
        criterion(q_, F(), M(), p_, nLeaf(), refine_tol, deref_tol, err_.get(), code_.get(), changed_.get());
        int changed = 0;
        changed_.get(&changed, 1);
        if (!changed && !reflag_) return false;

        err_leaf_.resize(nLeaf());
        err_.get(err_leaf_.data(), err_leaf_.size());
        for (int m = 0; m < nLeaf(); ++m) err_slot_[tree_.leaves()[m]] = err_leaf_[m];

        const auto change = tree_.adapt(err_slot_, refine_tol, deref_tol);
        reflag_ = change.has_value();       // a new mesh: the bits belong to the old leaves
        if (!change) return false;
        if (!change->restrict_pairs.empty())
            restrict_children(q_, F(), list_.put(change->restrict_pairs), int(change->restrict_pairs.size() / 2));
        if (!change->prolong_pairs.empty()) {
            clear_slots(q_, F(), list_.put(change->cleared), int(change->cleared.size()));
            prolong(q_, F(), list2_.put(change->prolong_pairs), int(change->prolong_pairs.size() / 2));
            moments(q_, F(), p_, put_refs(change->new_blocks), int(change->new_blocks.size()));
        }
        mesh_.upload(tree_);
        return true;
    }

};

void Solver::run(int argc, char **argv) {
    using clock = std::chrono::steady_clock;
    const auto t_start = clock::now();
    const char *s = std::getenv("AMR_GPU_SYNC");
    const bool sync = s && *s && *s != '0';
    std::printf("\nGPU: %s\n", q_.get_device().get_info<sycl::info::device::name>().c_str());
    std::printf("\nCase: rising bubble (Fakhari, Geier & Lee, JCP 2016, Sec. 5.1)\n");
    p_ = case_.setup();
    const char *o = std::getenv("AMR_OUTPUT");
    const std::string writer = o && *o ? o : "vtk";
    if (writer == "vtk") vtk_.emplace();
    else if (writer != "none") throw std::runtime_error("AMR_OUTPUT: unknown writer " + writer + " (vtk, none)");
    std::printf("\nOutput writer: %s\n", writer.c_str());
    std::printf("\nBlock pool: %d blocks (%.0f MiB on the device)\n", MaxBlks,
                MaxBlks * DeviceFields::mib_per_block());

    //== domain, initial state, refinement level by level from the initial state
    tree_.create_domain();
    mesh_.upload(tree_);
    initialize_state();
    for (int L = minLevels; L < maxLevels; ++L) {
        adapt_mesh();
        initialize_state();
    }
    q_.wait();
    const auto t_init = clock::now();

    //== steps and intervals (Solver_Mod: Main_Loop)
    const int Tf = argc >= 2 ? std::atoi(argv[1]) : int(std::lround(case_.t_end / case_.time_per_step / 20)) * 20;
    const int nOutStep = argc >= 3 ? std::atoi(argv[2]) : std::max(Tf / 20, 1);
    if (Tf < 1 || nOutStep < 1) throw std::runtime_error("usage: amr_lbm_sycl [steps [output_interval]]");
    const int nDiagStep = std::max(Tf / 120, 1), nProgStep = std::max(Tf / 20, 1);

    std::printf("\nModel: mass-conserving phase-field LBM (Fakhari, Geier & Lee, JCP 2016)\n");
    std::printf("\n    Tf     Sigma    W      M0    tauc   taul   tauh     Ref_Tols\n\n");
    std::printf("%8d%8.4f%6.2f%7.3f%7.3f%7.3f%7.3f%8.4f%8.4f\n", Tf, case_.sigma(), p_.W, p_.M0, 1 / p_.wc,
                p_.taul, p_.tauh, refine_tol, deref_tol);
    std::printf("\n     T   Nblks  nLeaf  Unused  Rho_min    Rho_max   U_max      Mass     Levels\n\n");

    StageTimers timed(q_, sync);
    auto tp_last = clock::now();
    long long cell_updates = 0;

    //== progress after Tdone steps: MLUPS over the last interval, time left
    auto progress = [&](int Tdone) {
        q_.wait();
        const auto now = clock::now();
        const double el = std::chrono::duration<double>(now - t_init).count();
        const double mlups = cell_updates / std::max(std::chrono::duration<double>(now - tp_last).count(), 1e-9) / 1e6;
        std::printf(" step%9d /%9d (%5.1f%%)  t* =%8.3f  leaf cells%9d  MLUPS%8.2f  elapsed %7.1f s  left ~ %7.1f s\n",
                    Tdone, Tf, 100.0 * Tdone / Tf, double(Tdone * case_.time_per_step), nLeaf() * NX * NY,
                    mlups, el, el / Tdone * (Tf - Tdone));
        std::fflush(stdout);
        tp_last = now;
        cell_updates = 0;
    };

    for (int Tn = 0; Tn <= Tf; ++Tn) {
        timed("Ghost_C", [&] { ghost_cells_C(); });
        if (Tn > 0) timed("refine", [&] { if (adapt_mesh()) ghost_cells_C(); });

        bool moments_done = (Tn == 0);
        const bool is_output = (Tn % nOutStep == 0) || Tn == Tf;
        if (is_output || Tn % nDiagStep == 0) {
            timed("Output", [&] {
                if (!moments_done) {
                    const auto r = tree_.leaf_refs();
                    moments(q_, F(), p_, put_refs(r), int(r.size()));
                    moments_done = true;
                }
                host_.fetch(fields_, tree_.nblks(), is_output && vtk_);
                const real time = real(Tn) * case_.time_per_step;
                if (is_output) {
                    if (vtk_) vtk_->step(time, tree_, host_);
                    case_.print_line(time, tree_, host_.C, host_.Uy);
                }
                case_.rise(time, tree_, host_.C, host_.Uy);
            });
            if (Tn == Tf) break;
        }

        timed("collision", [&] { collide(q_, F(), M(), p_, nLeaf(), !moments_done); });
        timed("Ghost_hg", [&] { ghost_cells_hg(); });
        timed("advection", [&] {
            stream(q_, F(), fields_.next(), M(), p_, nLeaf());
            fields_.swap();
        });

        cell_updates += (long long)nLeaf() * NX * NY;
        if ((Tn + 1) % nProgStep == 0 && Tn + 1 < Tf) progress(Tn + 1);
    }
    q_.wait();

    timed.report();
    std::printf(" t_initialization:%12.3f\n", std::chrono::duration<double>(t_init - t_start).count());
    std::printf("\n *****************************************\n");
    std::printf(" Time Elapsed:%12.1f\n", std::chrono::duration<double>(clock::now() - t_start).count());
    std::printf(" *****************************************\n");
}

int main(int argc, char **argv) {
    try {
        Solver solver;
        solver.run(argc, argv);
    } catch (const sycl::exception &e) {
        std::fprintf(stderr, "SYCL error: %s\n", e.what());
        return 1;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
