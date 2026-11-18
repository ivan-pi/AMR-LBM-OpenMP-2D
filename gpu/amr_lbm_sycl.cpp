//==============================================================================
//  AMR-LBM on a GPU (SYCL): the models of the Fortran solver (src/), single
//  precision
//
//      JCP 2016   mass-conserving phase-field model (model_jcp2016.hpp):
//                 rising bubble, Rayleigh-Taylor instability
//      Parker     multicomponent model with chemical-potential driven
//                 diffusion (model_parker.hpp): two species (droplets,
//                 spinodal decomposition) and three (an ouzo-like mixture of
//                 an oil, ethanol and water, diluted with water)
//
//  The host keeps the quadtree of blocks (metadata only) and runs the mesh
//  logic of Mesh_Mod / Refine_Mod / Blocks_Mod step for step: the same slot
//  allocation, Z-order leaf lists, neighbour table and flagging (with the
//  refinement of coarse islands), so a run can be compared with the Fortran
//  program built with AMR_PRECISION=4. All field work runs on the device:
//  the time step, the initial state, prolongation and restriction. Per step
//  the host fetches the refinement criterion of the leaves; after a mesh
//  change it uploads the block metadata, at diagnostic steps it fetches the
//  scalars and the velocity.
//
//  The pieces:
//      DeviceArray<T>       a USM device allocation (RAII)
//      DeviceFields<F>      the field pool of a model (FieldsT<NF, NC>)
//      DeviceMesh           the block metadata and leaf list on the device
//      Quadtree             the mesh on the host; adapt() returns a MeshChange,
//                           the field work it implies
//      Case<Model>          a test case: parameters (namelist), initial state,
//                           diagnostics, output fields; RisingBubble,
//                           RayleighTaylor (JCP 2016), ParkerBinary,
//                           OuzoCase (Parker, 2 and 3 species)
//      HostFields<NC>       a copy of the scalars, P, Ux, Uy on the host
//      VtkSeries            VTK XML output, the files of the Fortran VTK writer
//      FrameStream          in-situ rendering: frames to a renderer process
//      StageTimers          wall time per stage of the time step
//      Solver<Model>        the time loop
//
//  Usage:  amr_lbm_sycl [number_of_time_steps [output_interval_in_steps]]
//          AMR_CASE=case.nml   namelist &rising_bubble, &rayleigh_taylor, &parker
//                              (as the Fortran code) or &ouzo; its group selects
//                              the case and the model
//          AMR_GPU_SYNC=1      wait after every stage (stage timings)
//          AMR_OUTPUT=vtk|none field output (default vtk)
//          AMR_RENDER="cmd"    stream frames to a renderer at the output steps,
//                              e.g. "python3 tools/render_stream.py -o run.mp4"
//  Without a GPU the program runs on a CPU device (e.g. Intel's OpenCL CPU
//  runtime); ONEAPI_DEVICE_SELECTOR picks the device. The host thread waits
//  for the device many times per step: on a loaded machine pin it to an idle
//  CPU (taskset -c <cpu>), or every wait costs a scheduler time slice. At the
//  output steps it writes vtk/Rising_Bubble_2D_NNNN.vtu, vtk/AMR_Blocks_NNNN.vtu
//  with the series Rising_Bubble_2D.pvd, AMR_Blocks.pvd (as VTK_Output_Mod),
//  and each case its own record (bubble_rise.dat, rt_fronts.dat, droplet.dat,
//  spinodal.dat, ouzo.dat).
//==============================================================================
#include "model_jcp2016.hpp"
#include "model_parker.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <condition_variable>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <variant>
#include <vector>

#if !defined(AMR_L0) || !defined(AMR_MIN_LEVEL)
#error "AMR_L0 and AMR_MIN_LEVEL must be defined"
#endif

using namespace amr;

//==============================================================================
//  Configuration (Config_Mod)
//==============================================================================
constexpr int L0 = AMR_L0, minLevels = AMR_MIN_LEVEL;
#ifndef AMR_DOMAIN_HEIGHT
#define AMR_DOMAIN_HEIGHT 3
#endif
constexpr int nBaseX = 1, nBaseY = AMR_DOMAIN_HEIGHT;     // domain L0 x nBaseY*L0
#ifdef AMR_PERIODIC_Y
constexpr bool periodicY = true;                           // periodic in y as well
#else
constexpr bool periodicY = false;                          // walls at the bottom and top
#endif
constexpr int ilog2(int v) { return v <= 1 ? 0 : 1 + ilog2(v / 2); }
constexpr int maxLevels = ilog2(L0 / NX);
constexpr int MaxBlks = nBaseX * nBaseY * ((1 << (2 * (maxLevels + 1))) - 1) / 3;
static_assert((NX << maxLevels) == L0, "L0 must be a power-of-two multiple of the block size");
static_assert(maxLevels < MAXLEV, "too many levels");

inline int dx_level(int L) { return (L0 / NX) >> L; }

//== the inverse cell sizes of the levels (Config_Mod: inv_dx_level)
inline void set_levels(Levels &p) {
    p.maxLevels = maxLevels;
    for (int L = 0; L <= maxLevels; ++L) p.inv_dx[L] = real(NX) / real(L0) * std::ldexp(real(1), L);
}

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

//== the field pool of a model, slots 1..cap. The distributions are
//   double-buffered: the streaming reads neighbouring blocks directly, so it
//   writes the other buffer (next()) and the two are swapped afterwards.
//   Their ghost rings (lbm_kernels.hpp: Pool) are single: the streaming
//   writes interiors only
template <class F>
struct DeviceFields {
    static constexpr int NF = F::NF, NC = F::NC;
    std::vector<DeviceArray<real>> f0, f1, gh, C;
    DeviceArray<real> P, Ux, Uy;
    int cur = 0;

    DeviceFields(sycl::queue &q, size_t cap) : P(q, PB * cap), Ux(q, PB * cap), Uy(q, PB * cap) {
        f0.reserve(NF); f1.reserve(NF); gh.reserve(NF); C.reserve(NC);
        for (int s = 0; s < NF; ++s) {
            f0.emplace_back(q, HB * cap);
            f1.emplace_back(q, HB * cap);
            gh.emplace_back(q, GB * cap);
        }
        for (int s = 0; s < NC; ++s) C.emplace_back(q, CB * cap);
    }

    F view() const { return fields(cur); }
    F next() const { return fields(1 - cur); }
    void swap() { cur = 1 - cur; }
    static double mib_per_block() { return double(2 * NF * HB + NF * GB + NC * CB + 3 * PB) * sizeof(real) / (1 << 20); }

private:
    F fields(int b) const {
        F x{};
        for (int s = 0; s < NF; ++s) x.f[s] = {b ? f1[s].get() : f0[s].get(), gh[s].get()};
        for (int s = 0; s < NC; ++s) x.C[s] = C[s].get();
        x.P = P.get(); x.Ux = Ux.get(); x.Uy = Uy.get();
        return x;
    }
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
//   (parent, first new child), moments of the new blocks
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

    //== Refine_Mod: Refine_Derefine_Mesh with the criterion err (per slot) and
    //   the refinement of coarse islands of at most max_island leaves (0: off);
    //   returns the field work, empty if nothing changed
    std::optional<MeshChange> adapt(const std::vector<real> &err, real refine_tol, real deref_tol, int max_island) {
        if (!flag(err, refine_tol, deref_tol, max_island)) return std::nullopt;
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

    //== Mesh_Mod: Rebuild_Neighbor_Table (periodic in x, and in y with
    //   periodicY, else walls -k at the bottom and top): the same level, else
    //   the coarser block covering the position
    void neighbours() {
        auto live = [&](auto &&f) { for (int n = 1; n <= nblks_; ++n) if (B_[n].level >= 0) f(n, B_[n]); };
        live([&](int n, Block &b) { at_[slot(b.level, b.I, b.J)] = n; });
        live([&](int, Block &b) {
            const int L = b.level, w = nBaseX * (1 << L), h = nBaseY * (1 << L);
            b.coarser = false;
            b.finer = false;
            for (int k = 1; k <= 8; ++k) {
                int J = b.J + dj[k];
                if (periodicY) J = (J % h + h) % h;
                if (J < 0 || J >= h) { b.nbr[k] = -k; continue; }
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
    bool flag(const std::vector<real> &err, real refine_tol, real deref_tol, int max_island) {
        for (int n = 1; n <= nblks_; ++n)
            B_[n].refin = B_[n].leaf && B_[n].level < maxLevels && err[n] > refine_tol;
        if (max_island > 0) coarse_islands(max_island, false);
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
        //== parents that would form a small island of leaves enclosed by finer
        //   blocks keep their children (they would be refined again at the next step)
        if (max_island > 0) coarse_islands(max_island, true);
        for (int L = levelMax_ - 1; L >= levelMin_ - 1; --L)       // finest first
            for (int n = 1; n <= nblks_; ++n)
                if (B_[n].level == L && B_[n].deref) keep_jump(n);

        return std::any_of(&B_[1], &B_[1] + nblks_, [](const Block &b) { return b.refin || b.deref; });
    }

    //== Refine_Mod: Refine_Coarse_Islands. Connected groups (through the four
    //   sides) of at most max_island leaves on one level below the finest
    //   that touch only finer blocks, e.g. the interior of a large droplet
    //   whose gradient-based criterion is zero (a group that reaches a wall or
    //   a coarser block is not an island), are flagged for refinement.
    //   after_deref: the mesh after the planned derefinement, where parents
    //   flagged deref are leaves; an island made of them keeps its children
    void coarse_islands(int max_island, bool after_deref) {
        auto leaf_after = [&](int n) { return (B_[n].leaf && !B_[n].refin) || (after_deref && B_[n].deref); };
        auto coarse_leaf = [&](int n) { return B_[n].level >= 0 && B_[n].level < maxLevels && leaf_after(n); };
        std::vector<int> comp(nblks_ + 1, 0), stack, members;
        int ic = 0;
        for (int n = 1; n <= nblks_; ++n) {
            if (comp[n] != 0 || !coarse_leaf(n)) continue;
            comp[n] = ++ic;
            stack.assign(1, n);
            members.clear();
            bool enclosed = true;
            while (!stack.empty()) {
                const int m = stack.back();
                stack.pop_back();
                members.push_back(m);
                for (int k = 1; k <= 4; ++k) {
                    const int b = B_[m].nbr[k];
                    if (b <= 0) enclosed = false;                               // a wall: the far field
                    else if (B_[b].level < B_[m].level) enclosed = false;       // a coarser neighbour
                    else if (leaf_after(b) && comp[b] == 0) {                   // same-level leaf: same island
                        comp[b] = ic;
                        stack.push_back(b);
                    }                                                           // else finer
                }
            }
            if (!enclosed || int(members.size()) > max_island) continue;
            for (int m : members) {
                if (B_[m].deref) B_[m].deref = false;
                else B_[m].refin = true;
            }
        }
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
//   the leaf lists of the ghost passes: hg1 the leaves whose distribution
//   ghost cells are copied (finer neighbour) or set by the walls, coarse
//   those with a coarser neighbour (the interpolating second passes)
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
//  Host copy of the fields of blocks 1..nblks needed for output: the scalars
//  and Uy at every diagnostic step, P and Ux as well at the output steps (or
//  always, for cases that need them)
//==============================================================================
template <int NC>
struct HostFields {
    std::array<std::vector<real>, NC> C;
    std::vector<real> P, Ux, Uy;

    template <class D> void fetch(const D &d, int nblks, bool all) {
        for (int s = 0; s < NC; ++s) {
            C[s].resize(CB * nblks);
            d.C[s].get(C[s].data(), C[s].size());
        }
        Uy.resize(PB * nblks);
        d.Uy.get(Uy.data(), Uy.size());
        if (!all) return;
        P.resize(PB * nblks);
        Ux.resize(PB * nblks);
        d.P.get(P.data(), P.size());
        d.Ux.get(Ux.data(), Ux.size());
    }
};

//==============================================================================
//  Namelists: the file named by AMR_CASE holds one group &name ... /, with
//  "name = value" items (reals, integers, quoted strings), comments (!) and
//  Fortran d exponents
//==============================================================================
namespace namelist {

using Var = std::variant<real *, int *, std::string *>;

inline std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return char(std::tolower(c)); });
    return s;
}

//== the name of the first group in the AMR_CASE file, "" without a file
inline std::string group() {
    const char *file = std::getenv("AMR_CASE");
    if (!file || !*file) return "";
    std::ifstream in(file);
    if (!in) throw std::runtime_error("AMR_CASE: cannot open the case file");
    for (std::string line; std::getline(in, line);) {
        if (auto p = line.find('!'); p != std::string::npos) line.erase(p);
        if (auto p = line.find('&'); p != std::string::npos) {
            std::string g = lower(line.substr(p + 1));
            return g.substr(0, g.find_first_of(" \t\r/"));
        }
    }
    return "";
}

//== the values of group grp into vars (names in lower case)
inline void read(const std::string &grp, const std::map<std::string, Var> &vars) {
    const char *file = std::getenv("AMR_CASE");
    if (!file || !*file) return;
    std::ifstream in(file);
    auto trim = [](std::string s) {
        s.erase(0, s.find_first_not_of(" \t\r"));
        s.erase(s.find_last_not_of(" \t\r") + 1);
        return s;
    };
    std::string body;               // the group's items, separated by ','
    bool inside = false;
    for (std::string line; std::getline(in, line);) {
        if (auto p = line.find('!'); p != std::string::npos) line.erase(p);
        if (!inside) {
            const auto p = lower(line).find("&" + grp);
            if (p == std::string::npos) continue;
            inside = true;
            line.erase(0, p + 1 + grp.size());
        }
        if (const auto p = line.find('/'); p != std::string::npos) { body += line.substr(0, p); break; }
        body += line + ",";
    }
    if (!inside) throw std::runtime_error("AMR_CASE: no namelist &" + grp);
    std::istringstream items(body);
    for (std::string item; std::getline(items, item, ',');) {
        const auto eq = item.find('=');
        if (eq == std::string::npos) continue;
        const std::string key = lower(trim(item.substr(0, eq)));
        std::string val = trim(item.substr(eq + 1));
        const auto v = vars.find(key);
        if (v == vars.end()) throw std::runtime_error("AMR_CASE: unknown name " + key + " in &" + grp);
        if (auto *s = std::get_if<std::string *>(&v->second)) {
            if (!val.empty() && (val.front() == '\'' || val.front() == '"')) val = val.substr(1, val.size() - 2);
            **s = lower(trim(val));
            continue;
        }
        std::replace(val.begin(), val.end(), 'd', 'e');
        std::replace(val.begin(), val.end(), 'D', 'e');
        const double x = std::strtod(val.c_str(), nullptr);
        if (auto *r = std::get_if<real *>(&v->second)) **r = real(x);
        else **std::get_if<int *>(&v->second) = int(std::lround(x));
    }
    std::printf("\nCase file: %s\n", file);
}

} // namespace namelist

//==============================================================================
//  Test cases (Case_Mod), in single precision as the Fortran programs
//  built with AMR_PRECISION=4: parameters, initial state, diagnostics and
//  the fields of the output
//==============================================================================
using File = std::unique_ptr<FILE, int (*)(FILE *)>;
inline File open_record(const char *name, const char *header) {
    File f(std::fopen(name, "w"), &std::fclose);
    if (!f) throw std::runtime_error(std::string("cannot write ") + name);
    std::fprintf(f.get(), "%s\n", header);
    return f;
}

template <class F> void for_leaf_cells(const Quadtree &T, F &&f) {
    for (int n : T.leaves()) {
        const int dx = dx_level(T[n].level);
        for (int j = 1; j <= NY; ++j)
            for (int i = 1; i <= NX; ++i) f(n, i, j, dx);
    }
}

template <class Model>
class Case {
public:
    using Params = typename Model::Params;
    using Fields = typename Model::Fields;
    using Host = HostFields<Fields::NC>;
    //== a cell field of the output: its name and value at cell (i, j) of block n
    struct CellField {
        std::string name;
        std::function<real(const Host &, int n, int i, int j)> value;
    };

    virtual ~Case() = default;
    virtual std::string name() const = 0;
    //== the physical parameters for the kernels (reads the case's namelist);
    //   sets t_end, time_per_step and possibly the mesh policy
    virtual Params setup() = 0;
    //== the initial state of the listed blocks, the distributions at equilibrium
    virtual void initial_state(sycl::queue &q, Fields F, const BlockRef *list, int count) const = 0;
    //== the parameter table above the progress lines
    virtual void print_parameters(int Tf) const = 0;
    //== one line of the progress table (Solver_Mod: Print_Diagnostics)
    virtual void print_line(real time, const Quadtree &T, const Host &f) const = 0;
    //== the case's own record of the run, at every diagnostic step
    virtual void diagnostics(real time, const Quadtree &T, const Host &f) = 0;
    //== the scalar cell fields of the VTK output (with Velocity, Level, Block)
    //   and of the in-situ rendering
    virtual std::vector<CellField> vtk_fields() const = 0;
    virtual std::vector<CellField> render_fields() const = 0;
    //== the diagnostics read P and Ux too
    virtual bool diagnostics_need_all() const { return false; }

    real t_end = 0, time_per_step = 0;
    //== mesh policy (Refine_Mod), the model's defaults unless the case changes them
    real refine_tol = Model::refine_tol, deref_tol = Model::deref_tol;
    int max_island = Model::max_island;

protected:
    Params p_{};

    static real speed(const Host &f, int n, int i, int j) {
        const real u = f.Ux[pix(i, j, n)], v = f.Uy[pix(i, j, n)];
        return std::sqrt(u * u + v * v);
    }
};

//------------------------------------------------------------------------------
//  The cases of the JCP 2016 model
//------------------------------------------------------------------------------
class JcpCase : public Case<JCP2016> {
public:
    real sigma = 0;

    void print_parameters(int Tf) const override {
        std::printf("\n    Tf     Sigma    W      M0    tauc   taul   tauh     Ref_Tols\n\n");
        std::printf("%8d%8.4f%6.2f%7.3f%7.3f%7.3f%7.3f%8.4f%8.4f\n", Tf, sigma, p_.W, p_.M0, 1 / p_.wc,
                    p_.taul, p_.tauh, refine_tol, deref_tol);
        std::printf("\n     T   Nblks  nLeaf  Unused  Rho_min    Rho_max   U_max      Mass     Levels\n\n");
    }

    //== from C and Uy (Solver_Mod: Print_Diagnostics)
    void print_line(real time, const Quadtree &T, const Host &f) const override {
        real rmin = 10, rmax = -10, umax = 0;
        double mass = 0;
        for_leaf_cells(T, [&](int n, int i, int j, int dx) {
            const real rho = p_.Rhol + f.C[0][cix(i, j, n)] * (p_.Rhoh - p_.Rhol);
            rmin = std::min(rmin, rho);
            rmax = std::max(rmax, rho);
            mass += double(rho) * dx * dx;
            umax = std::max(umax, std::fabs(f.Uy[pix(i, j, n)]));
        });
        std::printf("%7.2f%7d%7d%6d%11.5f%11.5f%9.5f%12.1f%4d%4d\n", time, T.nblks(), T.nLeaf(), T.unused(),
                    rmin, rmax, umax, mass, T.levelMin(), T.levelMax());
    }

    std::vector<CellField> vtk_fields() const override {
        return {{"C", [](const Host &f, int n, int i, int j) { return f.C[0][cix(i, j, n)]; }},
                {"P", [](const Host &f, int n, int i, int j) { return f.P[pix(i, j, n)]; }}};
    }
    std::vector<CellField> render_fields() const override {
        const real rl = p_.Rhol, rh = p_.Rhoh;
        return {{"C", [](const Host &f, int n, int i, int j) { return f.C[0][cix(i, j, n)]; }},
                {"rho", [rl, rh](const Host &f, int n, int i, int j) { return rl + f.C[0][cix(i, j, n)] * (rh - rl); }},
                {"speed", speed}};
    }

protected:
    //== the model parameters from the fluid properties (common part of the Setups)
    void fluid_parameters(real W, real M0, real mu_h, real mu_l) {
        p_.W = W;
        p_.M0 = M0;
        p_.wc = real(1) / (real(0.5) + M0 * real(3));
        p_.tauh = real(3) * (mu_h / p_.Rhoh);
        p_.taul = real(3) * (mu_l / p_.Rhol);
        p_.Beta = real(12) * sigma / W;
        p_.kappa = real(1.5) * sigma * W;
        set_levels(p_);
    }
};

//== rising bubble (Rising_Bubble_Mod): a bubble of diameter D = L0/5 at
//   (L0/2, L0/2); bubble_rise.dat
class RisingBubble : public JcpCase {
public:
    std::string name() const override { return "rising bubble (Fakhari, Geier & Lee, JCP 2016, Sec. 5.1)"; }

    Params setup() override {
        R_ = real(L0) / real(10);
        D_ = real(2) * R_;
        t_end = 8;
        namelist::read("rising_bubble", {{"eo", &Eo_}, {"re", &Re_}, {"rho_ratio", &rho_ratio_},
                                         {"mu_ratio", &mu_ratio_}, {"gy_lu", &Gy_lu_}, {"mobility", &Mobility_},
                                         {"pe", &Pe_}, {"width", &Width_}, {"cn", &Cn_}, {"t_end", &t_end}});
        p_.Rhoh = 1;
        p_.Rhol = p_.Rhoh / rho_ratio_;
        p_.Gy = Gy_lu_;
        const real W = Cn_ > 0 ? Cn_ * real(L0) : Width_;
        const real M0 = Pe_ > 0 ? std::sqrt(p_.Gy * D_) * W / Pe_ : Mobility_;
        sigma = p_.Gy * (p_.Rhoh - p_.Rhol) * D_ * D_ / Eo_;
        mu_h_ = std::sqrt(p_.Gy * p_.Rhoh * (p_.Rhoh - p_.Rhol) * D_ * D_ * D_) / Re_;
        const real mul = mu_h_ / mu_ratio_;
        fluid_parameters(W, M0, mu_h_, mul);
        time_per_step = std::sqrt(p_.Gy / D_);
        std::printf("\nEo, Re_h, muh, mul, Pe =%12.3f%12.3f%12.3f%12.3f%12.3f\n", Eo_, Re_, mu_h_, mul,
                    std::sqrt(p_.Gy * D_) * W / M0);
        return p_;
    }

    void initial_state(sycl::queue &q, Fields F, const BlockRef *list, int count) const override {
        jcp::initial_state(q, F, p_, jcp::BubbleInit{R_, real(L0 / 2), real(L0 / 2), sigma, p_.W}, L0, list, count);
    }

    //== Rising_Bubble_Mod: Bubble_Diagnostics (cells with C >= 0.9 left out)
    void diagnostics(real time, const Quadtree &T, const Host &f) override {
        if (!rise_) rise_ = open_record("bubble_rise.dat", "#      t*           y_c          u_r          Re_r");
        double A = 0, Ay = 0, Av = 0;
        for_leaf_cells(T, [&](int n, int i, int j, int dx) {
            const real c = f.C[0][cix(i, j, n)];
            if (c >= real(0.9)) return;
            const double w = double(real(1) - std::max(c, real(0))) * dx * dx;
            A += w;
            Ay += w * (double(T[n].J) * dx * NY + dx * (j - 0.5));
            Av += w * f.Uy[pix(i, j, n)];
        });
        std::fprintf(rise_.get(), "%13.5E%13.5E%13.5E%13.5E\n", double(time), Ay / A, Av / A,
                     double(p_.Rhoh) * (Av / A) * D_ / mu_h_);
        std::fflush(rise_.get());
    }

private:
    real Eo_ = 40, Re_ = 40, rho_ratio_ = 1000, mu_ratio_ = 100, Gy_lu_ = real(5e-6);
    real Mobility_ = real(0.03), Pe_ = 0, Width_ = 5, Cn_ = 0;
    real R_ = 0, D_ = 0, mu_h_ = 0;
    File rise_{nullptr, &std::fclose};
};

//== Rayleigh-Taylor instability (Rayleigh_Taylor_Mod; Fakhari et al., PRE
//   2017, Sec. IV B): heavy fluid over light, the interface at mid-height
//   perturbed by A L0 cos(2 pi x / L0). U0 = sqrt(g L0), t0 = sqrt(L0/(g At)):
//   At, Re = rho_h U0 L0/mu_h, Ca = mu_h U0/sigma, Pe = U0 L0/M, t* = t/t0.
//   rt_fronts.dat: bubble front (x = 0) and spike tip (x = L0/2) relative to
//   the unperturbed interface, in units of L0
class RayleighTaylor : public JcpCase {
public:
    std::string name() const override { return "Rayleigh-Taylor instability (Fakhari et al., PRE 2017, Sec. IV B)"; }

    Params setup() override {
        t_end = real(2.5);
        namelist::read("rayleigh_taylor", {{"at", &At_}, {"re", &Re_}, {"mu_ratio", &mu_ratio_}, {"ca", &Ca_},
                                           {"gy_lu", &Gy_lu_}, {"amplitude", &Amplitude_}, {"mobility", &Mobility_},
                                           {"pe", &Pe_}, {"width", &Width_}, {"cn", &Cn_}, {"t_end", &t_end}});
        if (At_ <= 0 || At_ >= 1) throw std::runtime_error("Rayleigh-Taylor: At must lie in (0,1)");
        const real D = real(L0);
        Y0_ = real(0.5) * real(nBaseY * L0);
        p_.Rhoh = 1;
        p_.Rhol = p_.Rhoh * (1 - At_) / (1 + At_);
        p_.Gy = Gy_lu_;
        const real U0 = std::sqrt(p_.Gy * D);
        const real W = Cn_ > 0 ? Cn_ * real(L0) : Width_;
        const real M0 = Pe_ > 0 ? U0 * D / Pe_ : Mobility_;
        const real mu_h = p_.Rhoh * U0 * D / Re_, mu_l = mu_h / mu_ratio_;
        sigma = mu_h * U0 / Ca_;
        fluid_parameters(W, M0, mu_h, mu_l);
        std::printf("\nAt, Re, Ca, Pe =%7.4f%9.1f%7.3f%9.1f   tau_h, tau_l =%9.5f%9.5f   t0 =%9.1f\n", At_, Re_, Ca_,
                    U0 * D / M0, p_.tauh + real(0.5), p_.taul + real(0.5), std::sqrt(D / (p_.Gy * At_)));
        time_per_step = std::sqrt(p_.Gy * At_ / D);
        column_.assign(size_t(nBaseY) * L0, 0);
        return p_;
    }

    void initial_state(sycl::queue &q, Fields F, const BlockRef *list, int count) const override {
        const jcp::RTInit init{Y0_, Amplitude_ * real(L0), real(L0), p_.W, (p_.Rhoh - p_.Rhol) * p_.Gy};
        jcp::initial_state(q, F, p_, init, L0, list, count);
    }

    void diagnostics(real time, const Quadtree &T, const Host &f) override {
        if (!fronts_) fronts_ = open_record("rt_fronts.dat", "#      t*       bubble/L0     spike/L0");
        std::fprintf(fronts_.get(), "%13.5E%13.5E%13.5E\n", double(time), double((lowest_crossing(T, f, 0) - Y0_) / real(L0)),
                     double((lowest_crossing(T, f, L0 / 2) - Y0_) / real(L0)));
        std::fflush(fronts_.get());
    }

private:
    real At_ = real(0.5), Re_ = 3000, mu_ratio_ = 1, Ca_ = real(0.26), Gy_lu_ = real(2e-6);
    real Amplitude_ = real(0.1), Mobility_ = real(0.03), Pe_ = 1000, Width_ = 5, Cn_ = 0;
    real Y0_ = 0;
    std::vector<real> column_;          // C along a vertical line, finest-level resolution
    File fronts_{nullptr, &std::fclose};

    //== the lowest height at which C crosses 1/2 on the vertical line through
    //   the cells containing the finest-level position x (Lowest_Crossing):
    //   the leaf cells on the line spread to finest-level resolution, the
    //   crossing interpolated linearly between cell centres
    real lowest_crossing(const Quadtree &T, const Host &f, int x) {
        for (int n : T.leaves()) {
            const int dx = dx_level(T[n].level), xb = T[n].I * dx * NX;
            if (x < xb || x >= xb + dx * NX) continue;
            const int i = (x - xb) / dx + 1, yb = T[n].J * dx * NY;
            for (int j = 1; j <= NY; ++j)
                std::fill_n(&column_[size_t(yb + (j - 1) * dx)], dx, f.C[0][cix(i, j, n)]);
        }
        for (size_t k = 0; k < column_.size(); ++k)
            if (column_[k] >= real(0.5))
                return k == 0 ? real(0.5)
                              : real(k) - real(0.5) + (real(0.5) - column_[k - 1]) / (column_[k] - column_[k - 1]);
        return real(column_.size()) - real(0.5);       // no heavy fluid on the line
    }
};

//------------------------------------------------------------------------------
//  The cases of the Parker model. Lattice units of the finest level; one
//  time unit is one time step
//------------------------------------------------------------------------------
template <int NS>
class ParkerCaseBase : public Case<Parker<NS>> {
public:
    using Base = Case<Parker<NS>>;
    using typename Base::Host;
    using typename Base::CellField;

    void print_parameters(int Tf) const override {
        std::printf("\n      Tf       D     kappa      nu   omega_E  Ref_Tols          max_island\n\n");
        std::printf("%8d%9.4f%9.4f%9.4f%9.4f%10.4f%8.4f%6d\n", Tf, p_.Diff, p_.kappa, p_.nu, p_.omega_E,
                    this->refine_tol, this->deref_tol, this->max_island);
        std::printf("\n        T   Nblks  nLeaf  Unused     x1_min     x1_max    U_max      Mass_1   Levels\n\n");
    }

    //== x1 = c1/c min/max, |Uy| max, the integral of c1 (the columns of
    //   Print_Diagnostics; the Fortran code prints c1 itself, with c = 1)
    void print_line(real time, const Quadtree &T, const Host &f) const override {
        real xmin = 10, xmax = -10, umax = 0;
        double mass = 0;
        for_leaf_cells(T, [&](int n, int i, int j, int dx) {
            const real x1 = mole_fraction(f, 0, n, i, j);
            xmin = std::min(xmin, x1);
            xmax = std::max(xmax, x1);
            mass += double(f.C[0][cix(i, j, n)]) * dx * dx;
            umax = std::max(umax, std::fabs(f.Uy[pix(i, j, n)]));
        });
        std::printf("%10.1f%7d%7d%6d%11.5f%11.5f%9.5f%12.1f%4d%4d\n", time, T.nblks(), T.nLeaf(), T.unused(),
                    xmin, xmax, umax, mass, T.levelMin(), T.levelMax());
    }

    bool diagnostics_need_all() const override { return true; }

    std::vector<CellField> vtk_fields() const override {
        std::vector<CellField> v{{"C", [](const Host &f, int n, int i, int j) { return f.C[0][cix(i, j, n)]; }}};
        for (int s = 0; s < NS; ++s)
            v.push_back({"x_" + species_[s], [s](const Host &f, int n, int i, int j) { return mole_fraction(f, s, n, i, j); }});
        v.push_back({"P", [](const Host &f, int n, int i, int j) { return f.P[pix(i, j, n)]; }});
        return v;
    }
    std::vector<CellField> render_fields() const override {
        std::vector<CellField> v;
        for (int s = 0; s < NS; ++s)
            v.push_back({species_[s], [s](const Host &f, int n, int i, int j) { return mole_fraction(f, s, n, i, j); }});
        v.push_back({"speed", Base::speed});
        return v;
    }

protected:
    using Base::p_;
    std::array<std::string, NS> species_;
    //== namelist &parker / &ouzo: parameters of the model and the mesh policy
    real D_lu_ = real(0.2), kappa_ = 1, nu_lu_ = real(0.02), omega_E_ = 1, c_total_ = 1;
    int n_steps_ = 12000, max_island_ = -1;
    real refine_tol_ = -1, deref_tol_ = -1;

    static real mole_fraction(const Host &f, int s, int n, int i, int j) {
        real c = 0;
        for (int r = 0; r < NS; ++r) c += f.C[r][cix(i, j, n)];
        return f.C[s][cix(i, j, n)] / c;
    }

    void common_items(std::map<std::string, namelist::Var> &v) {
        v.insert({{"d_lu", &D_lu_}, {"kappa", &kappa_}, {"nu_lu", &nu_lu_}, {"omega_e", &omega_E_},
                  {"c_total", &c_total_}, {"n_steps", &n_steps_}, {"max_island", &max_island_},
                  {"refine_tol", &refine_tol_}, {"deref_tol", &deref_tol_}});
    }

    //== the model parameters and the mesh policy from the namelist items
    void model_parameters() {
        p_.Diff = D_lu_;
        p_.kappa = kappa_;
        p_.nu = nu_lu_;
        p_.omega_E = omega_E_;
        p_.wD = real(1) / (p_.Diff / parker::cs2 + real(0.5));
        p_.wV = real(1) / (p_.nu / parker::cs2 + real(0.5));
        p_.x_eps = real(1e-9);
        set_levels(p_);
        if (max_island_ >= 0) this->max_island = max_island_;
        if (refine_tol_ > 0) this->refine_tol = refine_tol_;
        if (deref_tol_ > 0) this->deref_tol = deref_tol_;
        this->t_end = real(n_steps_);
        this->time_per_step = 1;
        std::printf("mesh: refine_tol, deref_tol, max_island_blocks:%10.2E%10.2E%4d\n", this->refine_tol,
                    this->deref_tol, this->max_island);
    }

    //== ln gamma of the mole fractions x (double precision, host): the formulas
    //   of model_parker.hpp
    std::array<double, NS> ln_gamma(const std::array<double, NS> &x) const {
        std::array<double, NS> lg{};
        if constexpr (NS == 2) {
            const double xa = x[0], xb = 1 - xa, A12 = p_.A[0][1], A21 = p_.A[1][0];
            lg[0] = (A12 + 2 * (A21 - A12) * xa) * xb * xb;
            lg[1] = (A21 + 2 * (A12 - A21) * xb) * xa * xa;
        } else {
            double g = 0;
            std::array<double, NS> dg{};
            for (int a = 0; a < NS; ++a)
                for (int b = a + 1; b < NS; ++b) {
                    const double Aba = p_.A[b][a], Aab = p_.A[a][b];
                    g += x[a] * x[b] * (Aba * x[a] + Aab * x[b]);
                    dg[a] += x[b] * (2 * Aba * x[a] + Aab * x[b]);
                    dg[b] += x[a] * (Aba * x[a] + 2 * Aab * x[b]);
                }
            g += p_.C123 * x[0] * x[1] * x[2];
            dg[0] += p_.C123 * x[1] * x[2];
            dg[1] += p_.C123 * x[0] * x[2];
            dg[2] += p_.C123 * x[0] * x[1];
            for (int s = 0; s < NS; ++s) lg[s] = dg[s] - 2 * g;
        }
        return lg;
    }

    //== coexisting mole fractions on the binary edge of species a, b (the
    //   others absent), as fractions of a: equal activities x gamma of both
    //   (Newton; Parker_Cases_Mod: Coexistence)
    std::pair<double, double> coexistence(int a, int b) const {
        auto lgx = [&](double x1) {
            std::array<double, NS> x{};
            x[a] = x1;
            x[b] = 1 - x1;
            return ln_gamma(x);
        };
        auto residual = [&](const std::array<double, 2> &v) {
            const auto ga = lgx(v[0]), gb = lgx(v[1]);
            return std::array<double, 2>{std::log(v[0]) + ga[a] - std::log(v[1]) - gb[a],
                                         std::log(1 - v[0]) + ga[b] - std::log(1 - v[1]) - gb[b]};
        };
        //== Newton, the step halved to stay inside (0, 1); several start values
        //   (the first is Parker_Cases_Mod's)
        for (const double x0 : {0.05, 0.01, 0.2, 1e-3, 1e-4}) {
            std::array<double, 2> v{x0, 1 - x0};
            bool ok = false;
            for (int it = 0; it < 200 && !ok; ++it) {
                const auto F = residual(v);
                double J[2][2];
                for (int k = 0; k < 2; ++k) {
                    auto vp = v;
                    vp[k] += 1e-7 * std::min(v[k], 1 - v[k]) / 0.05;
                    const auto Fp = residual(vp);
                    J[0][k] = (Fp[0] - F[0]) / (vp[k] - v[k]);
                    J[1][k] = (Fp[1] - F[1]) / (vp[k] - v[k]);
                }
                const double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
                double d0 = (-J[1][1] * F[0] + J[0][1] * F[1]) / det, d1 = (J[1][0] * F[0] - J[0][0] * F[1]) / det;
                if (!std::isfinite(d0) || !std::isfinite(d1)) break;
                while (!(v[0] + d0 > 0 && v[0] + d0 < 1 && v[1] + d1 > 0 && v[1] + d1 < 1)) {
                    d0 *= 0.5;
                    d1 *= 0.5;
                }
                v[0] += d0;
                v[1] += d1;
                ok = std::max(std::fabs(d0), std::fabs(d1)) < 1e-13;
            }
            if (ok && v[0] < v[1] - 1e-6) return {v[0], v[1]};
        }
        throw std::runtime_error("Coexistence: no two-phase solution (Margules A too small?)");
    }

    //== the largest |u| of the leaf cells
    static real max_speed(const Quadtree &T, const Host &f) {
        real u2 = 0;
        for_leaf_cells(T, [&](int n, int i, int j, int) {
            u2 = std::max(u2, f.Ux[pix(i, j, n)] * f.Ux[pix(i, j, n)] + f.Uy[pix(i, j, n)] * f.Uy[pix(i, j, n)]);
        });
        return std::sqrt(u2);
    }
};

//== two species (Parker_Cases_Mod), namelist &parker:
//   kind = 'droplet'  : a droplet of radius R, rich in species 1, in the
//                       species-2-rich phase, elongated by eps along x
//                       (eps = 0: static droplet, Laplace test; eps > 0:
//                       n = 2 shape oscillation). droplet.dat: t,
//                       deformation (Ixx-Iyy)/(Ixx+Iyy), R_eq, p_in, p_out,
//                       mass of species 1, largest |u|.
//   kind = 'box'      : the whole domain filled with the mixture x1 = x_mix +
//                       noise. spinodal.dat like 'spinodal'.
//   kind = 'spinodal' : the middle L0 x L0 square filled with the mixture
//                       x1 = x_mix + noise; spinodal decomposition.
//                       spinodal.dat: t, number of leaf cells, mass of
//                       species 1, largest |u|.
//   The coexisting mole fractions of the Margules model are computed at
//   start-up (Newton) and used for the initial phases.
class ParkerBinary : public ParkerCaseBase<2> {
public:
    std::string name() const override { return "two-species droplet / spinodal decomposition (namelist &parker)"; }

    Params setup() override {
        std::map<std::string, namelist::Var> items{{"kind", &kind_}, {"a_12", &A12_}, {"a_21", &A21_},
                                                   {"radius", &Radius_}, {"eps", &eps_}, {"x_mix", &x_mix_},
                                                   {"noise", &noise_}};
        common_items(items);
        namelist::read("parker", items);
        if (kind_ == "droplet") kind_id_ = 0;
        else if (kind_ == "box") kind_id_ = 1;
        else if (kind_ == "spinodal") kind_id_ = 2;
        else throw std::runtime_error("namelist &parker: kind must be droplet, spinodal or box");
        species_ = {"1", "2"};
        p_.A[0][1] = A12_;
        p_.A[1][0] = A21_;
        p_.A[0][0] = p_.A[1][1] = 0;
        p_.C123 = 0;
        const auto [xa, xb] = coexistence(0, 1);
        xa_ = real(xa);
        xb_ = real(xb);
        X0_ = real(0.5) * real(nBaseX * L0);
        Y0_ = real(0.5) * real(nBaseY * L0);
        std::printf("\nkind:                    %s\n", kind_.c_str());
        std::printf("D, kappa, nu:            %10.5f%10.5f%10.5f\n", D_lu_, kappa_, nu_lu_);
        std::printf("Margules A12, A21:       %10.5f%10.5f\n", A12_, A21_);
        std::printf("coexistence x1 (a, b):   %10.5f%10.5f\n", xa_, xb_);
        std::printf("kappa D^2:               %12.4E\n", kappa_ * D_lu_ * D_lu_);
        model_parameters();
        return p_;
    }

    void initial_state(sycl::queue &q, Fields F, const BlockRef *list, int count) const override {
        const parker::BinaryInit init{kind_id_, xa_, xb_, X0_, Y0_, Radius_, eps_, x_mix_, noise_, real(L0)};
        parker::initial_state<2>(q, F, init, c_total_, L0, list, count);
    }

    //== Parker_Cases_Mod: Diagnostics
    void diagnostics(real time, const Quadtree &T, const Host &f) override {
        if (!rec_) {
            rec_ = kind_id_ == 0
                       ? open_record("droplet.dat", "#      t        deformation      R_eq          p_in         p_out        mass_1       max|u|")
                       : open_record("spinodal.dat", "#      t        leaf_cells     mass_1       max|u|");
        }
        double A = 0, Ixx = 0, Iyy = 0, Pin = 0, Nin = 0, Pout = 0, Nout = 0, Mass = 0;
        const double R = Radius_;
        for_leaf_cells(T, [&](int n, int i, int j, int dx) {
            const double area = double(dx) * dx;
            const double x = T[n].I * dx * NX + dx * (i - 0.5) - X0_, y = T[n].J * dx * NY + dx * (j - 0.5) - Y0_;
            const double r2 = x * x + y * y;
            const double x1 = double(f.C[0][cix(i, j, n)]) / (double(f.C[0][cix(i, j, n)]) + f.C[1][cix(i, j, n)]);
            //== droplet moments from r < 2.5 R only: the far field relaxes to the
            //   discrete coexistence value, slightly off xa
            double w = std::min(std::max((x1 - xa_) / (xb_ - xa_), 0.0), 1.0) * area;
            if (r2 > (2.5 * R) * (2.5 * R)) w = 0;
            A += w;
            Ixx += w * x * x;
            Iyy += w * y * y;
            Mass += double(f.C[0][cix(i, j, n)]) * area;
            const double P = f.P[pix(i, j, n)];
            if (r2 < (0.5 * R) * (0.5 * R)) {
                Pin += P * area;
                Nin += area;
            } else if (r2 > (R + 12) * (R + 12) && r2 < (2.5 * R) * (2.5 * R)) {
                Pout += P * area;
                Nout += area;
            }
        });
        const double umax = max_speed(T, f);
        if (kind_id_ == 0)
            std::fprintf(rec_.get(), "%10.1f%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E\n", double(time), (Ixx - Iyy) / (Ixx + Iyy),
                         std::sqrt(A / (4 * std::atan(1.0))), Pin / std::max(Nin, 1.0), Pout / std::max(Nout, 1.0), Mass, umax);
        else
            std::fprintf(rec_.get(), "%10.1f%12d%14.6E%14.6E\n", double(time), T.nLeaf() * NX * NY, Mass, umax);
        std::fflush(rec_.get());
    }

private:
    std::string kind_ = "droplet";
    int kind_id_ = 0;
    real A12_ = real(2.7627), A21_ = real(2.4474);
    real Radius_ = 20, eps_ = 0, x_mix_ = real(0.5), noise_ = real(0.05);
    real xa_ = 0, xb_ = 0, X0_ = 0, Y0_ = 0;
    File rec_{nullptr, &std::fclose};
};

//== three species, namelist &ouzo: an ouzo-like mixture of an oil (species 1,
//   e.g. trans-anethole), ethanol (2) and water (3) meeting pure water. The
//   oil mixes with ethanol, ethanol with water, but oil and water hardly
//   (Margules A_13 = A_31 = 4: coexistence x_oil = 0.021 / 0.979). The stock
//   (x1_stock, x2_stock) is one phase; diluted with water its mixing line
//   crosses the spinodal (tools/ternary_phase.py), so where ethanol diffuses
//   out into the water the oil comes out of solution as droplets
//   (spontaneous emulsification, the "ouzo effect").
//   kind = 'drop'  : a drop of stock of radius Radius in the centre
//   kind = 'layer' : the stock below the middle height, its edge perturbed
//                    by amplitude*L0 cos(2 pi x/L0)
//   kind = 'mix'   : the whole domain diluted to the fraction 'dilution' of
//                    stock (plus noise): the emulsion forms everywhere at once
//   Initial edges of width Width (tanh), and a relative perturbation of
//   amplitude noise of the oil and ethanol fractions on the finest cells.
//   ouzo.dat: t, leaf cells, the masses of the three species, max |u|, the
//   number of oil droplets (connected regions with x_oil > 1/2 on the
//   finest-level grid, at least 4 cells), their total area, the fraction of
//   the oil in them and their mean equivalent radius.
class OuzoCase : public ParkerCaseBase<3> {
public:
    std::string name() const override { return "ternary mixture: ouzo meeting water (namelist &ouzo)"; }

    Params setup() override {
        Radius_ = real(L0) / 4;
        std::map<std::string, namelist::Var> items{
            {"kind", &kind_}, {"a_12", &A_[0][1]}, {"a_21", &A_[1][0]}, {"a_13", &A_[0][2]}, {"a_31", &A_[2][0]},
            {"a_23", &A_[1][2]}, {"a_32", &A_[2][1]}, {"c_123", &C123_}, {"x1_stock", &x1_stock_},
            {"x2_stock", &x2_stock_}, {"radius", &Radius_}, {"width", &Width_}, {"amplitude", &amplitude_},
            {"dilution", &dilution_}, {"noise", &noise_}, {"species", &names_}};
        n_steps_ = 40000;
        common_items(items);
        namelist::read("ouzo", items);
        if (kind_ == "drop") kind_id_ = 0;
        else if (kind_ == "layer") kind_id_ = 1;
        else if (kind_ == "mix") kind_id_ = 2;
        else throw std::runtime_error("namelist &ouzo: kind must be drop, layer or mix");
        std::istringstream sn(names_);
        for (auto &s : species_)
            if (!std::getline(sn, s, ':')) throw std::runtime_error("namelist &ouzo: species = 'oil:ethanol:water'");
        for (int r = 0; r < 3; ++r)
            for (int s = 0; s < 3; ++s) p_.A[r][s] = r == s ? 0 : A_[r][s];
        p_.C123 = C123_;
        const real x3 = 1 - x1_stock_ - x2_stock_;
        if (x1_stock_ <= 0 || x2_stock_ < 0 || x3 < 0) throw std::runtime_error("namelist &ouzo: bad stock composition");
        X0_ = real(0.5) * real(nBaseX * L0);
        Y0_ = real(0.5) * real(nBaseY * L0);
        const auto [xa, xb] = coexistence(0, 2);
        const auto lg = ln_gamma({x1_stock_, x2_stock_, x3});
        std::printf("\nkind:                          %s\n", kind_.c_str());
        std::printf("species 1, 2, 3:               %s, %s, %s\n", species_[0].c_str(), species_[1].c_str(), species_[2].c_str());
        std::printf("D, kappa, nu:                  %10.5f%10.5f%10.5f\n", D_lu_, kappa_, nu_lu_);
        std::printf("Margules A12 A21 A13 A31 A23 A32, C123:%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n", A_[0][1], A_[1][0],
                    A_[0][2], A_[2][0], A_[1][2], A_[2][1], C123_);
        std::printf("coexistence x_%s on the %s-%s edge: %10.5f%10.5f\n", species_[0].c_str(), species_[0].c_str(),
                    species_[2].c_str(), xa, xb);
        std::printf("stock x:                       %10.5f%10.5f%10.5f   ln gamma:%9.4f%9.4f%9.4f\n", x1_stock_, x2_stock_,
                    x3, lg[0], lg[1], lg[2]);
        if (kind_id_ == 2) std::printf("dilution (stock fraction):     %10.5f\n", dilution_);
        std::printf("kappa D^2:                     %12.4E\n", kappa_ * D_lu_ * D_lu_);
        model_parameters();
        raster_.assign(size_t(nBaseX) * L0 * nBaseY * L0, 0);
        return p_;
    }

    void initial_state(sycl::queue &q, Fields F, const BlockRef *list, int count) const override {
        const real x3 = 1 - x1_stock_ - x2_stock_;
        const parker::OuzoInit init{kind_id_,  {x1_stock_, x2_stock_, x3}, {0, 0, 1}, X0_, Y0_, Radius_, Width_,
                                    amplitude_ * real(L0), dilution_, noise_, real(L0)};
        parker::initial_state<3>(q, F, init, c_total_, L0, list, count);
    }

    void diagnostics(real time, const Quadtree &T, const Host &f) override {
        if (!rec_)
            rec_ = open_record("ouzo.dat", "#      t     leaf_cells     mass_1        mass_2        mass_3        max|u|"
                                           "     droplets   oil_area   oil_in_drops  mean_radius");
        double mass[3] = {0, 0, 0};
        const int W = nBaseX * L0, H = nBaseY * L0;
        for (int n : T.leaves()) {
            const int dx = dx_level(T[n].level), x0 = T[n].I * dx * NX, y0 = T[n].J * dx * NY;
            for (int j = 1; j <= NY; ++j)
                for (int i = 1; i <= NX; ++i) {
                    real c = 0;
                    for (int s = 0; s < 3; ++s) {
                        mass[s] += double(f.C[s][cix(i, j, n)]) * dx * dx;
                        c += f.C[s][cix(i, j, n)];
                    }
                    const real x1 = f.C[0][cix(i, j, n)] / c;
                    for (int b = 0; b < dx; ++b)
                        std::fill_n(&raster_[size_t(y0 + (j - 1) * dx + b) * W + x0 + (i - 1) * dx], dx, x1);
                }
        }
        //== oil droplets: 4-connected regions of x1 > 1/2, periodic in x (and y)
        std::vector<int> label(raster_.size(), 0), stack;
        int nd = 0;
        double area = 0, oil = 0;
        for (size_t p0 = 0; p0 < raster_.size(); ++p0) {
            if (label[p0] || raster_[p0] <= real(0.5)) continue;
            stack.assign(1, int(p0));
            label[p0] = -1;
            size_t cells = 0;
            double xsum = 0;
            while (!stack.empty()) {
                const int p = stack.back();
                stack.pop_back();
                ++cells;
                xsum += raster_[p];
                const int x = p % W, y = p / W;
                const int nb[4][2] = {{x - 1, y}, {x + 1, y}, {x, y - 1}, {x, y + 1}};
                for (const auto &q : nb) {
                    int xx = (q[0] + W) % W, yy = q[1];
                    if (periodicY) yy = (yy + H) % H;
                    if (yy < 0 || yy >= H) continue;
                    const size_t pp = size_t(yy) * W + xx;
                    if (!label[pp] && raster_[pp] > real(0.5)) {
                        label[pp] = -1;
                        stack.push_back(int(pp));
                    }
                }
            }
            if (cells >= 4) {
                ++nd;
                area += double(cells);
                oil += xsum * c_total_;                       // c x1 per cell, c ~ c_total
            }
        }
        const double r_mean = nd ? std::sqrt(area / nd / (4 * std::atan(1.0))) : 0;
        std::fprintf(rec_.get(), "%10.1f%12d%14.6E%14.6E%14.6E%14.6E%8d%12.1f%14.6E%12.4f\n", double(time),
                     T.nLeaf() * NX * NY, mass[0], mass[1], mass[2], double(max_speed(T, f)), nd, area,
                     oil / mass[0], r_mean);
        std::fflush(rec_.get());
    }

private:
    std::string kind_ = "drop", names_ = "anethole:ethanol:water";
    int kind_id_ = 0;
    //== Margules A_rs: oil-ethanol, oil-water, ethanol-water
    real A_[3][3] = {{0, real(0.3), 4}, {real(0.3), 0, real(0.6)}, {4, real(0.6), 0}};
    real C123_ = 0;
    real x1_stock_ = real(0.25), x2_stock_ = real(0.3);
    real Radius_ = 0, Width_ = 4, amplitude_ = real(0.02), dilution_ = real(0.5), noise_ = real(0.02);
    real X0_ = 0, Y0_ = 0;
    std::vector<real> raster_;
    File rec_{nullptr, &std::fclose};
};

//==============================================================================
//  VTK XML output (VTK_Output_Mod: Write_VTU_Step): per output step the leaf
//  cells and the leaf blocks as UnstructuredGrids of quads in vtk/, raw
//  binary appended data with UInt64 size headers, and the two .pvd time
//  series rewritten every step (valid if the run stops). Coordinates /L0.
//==============================================================================
class VtkSeries {
public:
    template <class Host, class Fields>
    void step(real time, const Quadtree &T, const Host &f, const Fields &fields) {
        if (times_.empty()) std::filesystem::create_directories(dir_);
        times_.push_back(time);
        char id[8];
        std::snprintf(id, sizeof id, "%04d", int(times_.size()) - 1);
        cells(T, f, fields).write(dir_ + "/Rising_Bubble_2D_" + id + ".vtu");
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

    //== all leaf cells, with the case's scalar fields, Velocity, Level and Block
    template <class Host, class Fields>
    static UGrid cells(const Quadtree &T, const Host &f, const Fields &fields) {
        UGrid g;
        const size_t nCells = size_t(T.nLeaf()) * NX * NY;
        std::vector<std::vector<real>> vals(fields.size());
        for (auto &v : vals) v.reserve(nCells);
        std::vector<std::int32_t> level, block;
        level.reserve(nCells); block.reserve(nCells);
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
                    for (size_t m = 0; m < fields.size(); ++m) vals[m].push_back(fields[m].value(f, b, i, j));
                    g.vdata.insert(g.vdata.end(), {f.Ux[pix(i, j, b)], f.Uy[pix(i, j, b)], real(0)});
                    level.push_back(blk.level);
                    block.push_back(b);
                }
            p0 += (NX + 1) * (NY + 1);
        }
        for (size_t m = 0; m < fields.size(); ++m) g.rfields.emplace_back(fields[m].name, std::move(vals[m]));
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
//  In-situ rendering: at each output step a frame (cell fields of the leaves
//  and the block layout) is handed to a background thread that writes it to
//  the standard input of a renderer process, so the solver goes on while the
//  frame is drawn. At most two frames wait; then push() blocks. Frame (little
//  endian): "AMRF", int32 version 2, float64 time, int32 L0, nBaseX, nBaseY,
//  block size, number of leaves m, number of fields; per field a name of 16
//  bytes (blank padded); per leaf int32 x0, y0, dx (finest-level lattice
//  units); per field m x n x n float32 (x fastest). The fields are the case's
//  (JCP 2016: C, the density rho and the speed |u|; Parker: the mole
//  fractions of the species and the speed). The renderer sees end of file at
//  the end.
//==============================================================================
class FrameStream {
public:
    explicit FrameStream(const std::string &command) : pipe_(popen(command.c_str(), "w")) {
        if (!pipe_) throw std::runtime_error("AMR_RENDER: cannot start " + command);
        worker_ = std::thread([this] { run(); });
        std::printf("\nRenderer: %s\n", command.c_str());
    }
    ~FrameStream() {
        {
            std::lock_guard<std::mutex> lk(m_);
            done_ = true;
        }
        cv_.notify_all();
        worker_.join();
        if (const int rc = pclose(pipe_); rc != 0) std::fprintf(stderr, "AMR_RENDER: renderer exit status %d\n", rc);
    }
    FrameStream(const FrameStream &) = delete;
    FrameStream &operator=(const FrameStream &) = delete;

    template <class Host, class Fields>
    void push(real time, const Quadtree &T, const Host &f, const Fields &fields) {
        Frame fr;
        fr.time = time;
        const size_t ncell = size_t(T.nLeaf()) * NX * NY;
        fr.blocks.reserve(3 * size_t(T.nLeaf()));
        fr.names.reserve(fields.size());
        fr.fields.resize(fields.size());
        for (const auto &fd : fields) fr.names.push_back(fd.name);
        for (auto &v : fr.fields) v.reserve(ncell);
        for (int n : T.leaves()) {
            const int dx = dx_level(T[n].level);
            fr.blocks.insert(fr.blocks.end(), {T[n].I * dx * NX, T[n].J * dx * NY, dx});
            for (size_t m = 0; m < fields.size(); ++m)
                for (int j = 1; j <= NY; ++j)
                    for (int i = 1; i <= NX; ++i) fr.fields[m].push_back(fields[m].value(f, n, i, j));
        }
        std::unique_lock<std::mutex> lk(m_);
        cv_.wait(lk, [this] { return queue_.size() < 2; });
        queue_.push_back(std::move(fr));
        cv_.notify_all();
    }

private:
    struct Frame {
        double time = 0;
        std::vector<std::int32_t> blocks;
        std::vector<std::string> names;
        std::vector<std::vector<float>> fields;
    };
    FILE *pipe_;
    std::thread worker_;
    std::mutex m_;
    std::condition_variable cv_;
    std::deque<Frame> queue_;
    bool done_ = false, broken_ = false;

    void run() {
        for (;;) {
            Frame fr;
            {
                std::unique_lock<std::mutex> lk(m_);
                cv_.wait(lk, [this] { return done_ || !queue_.empty(); });
                if (queue_.empty()) return;             // done and drained
                fr = std::move(queue_.front());
                queue_.pop_front();
            }
            cv_.notify_all();
            if (!broken_) write(fr);
        }
    }

    void write(const Frame &fr) {
        const std::int32_t head[7] = {2, L0, nBaseX, nBaseY, NX, std::int32_t(fr.blocks.size() / 3),
                                      std::int32_t(fr.fields.size())};
        bool ok = std::fwrite("AMRF", 1, 4, pipe_) == 4 && std::fwrite(&head[0], 4, 1, pipe_) == 1 &&
                  std::fwrite(&fr.time, 8, 1, pipe_) == 1 && std::fwrite(&head[1], 4, 6, pipe_) == 6;
        for (const std::string &name : fr.names) {
            char padded[16];
            std::memset(padded, ' ', 16);
            std::memcpy(padded, name.data(), std::min<size_t>(16, name.size()));
            ok = ok && std::fwrite(padded, 1, 16, pipe_) == 16;
        }
        ok = ok && std::fwrite(fr.blocks.data(), 4, fr.blocks.size(), pipe_) == fr.blocks.size();
        for (const auto &v : fr.fields) ok = ok && std::fwrite(v.data(), 4, v.size(), pipe_) == v.size();
        ok = ok && std::fflush(pipe_) == 0;
        if (!ok) {
            broken_ = true;
            std::fprintf(stderr, "AMR_RENDER: the renderer stopped reading; no more frames\n");
        }
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
//  The solver (Solver_Mod) for a model: JCP2016 or Parker<NS>
//==============================================================================
//== device selection: a GPU if there is one, else a CPU (e.g. the OpenCL CPU
//   runtime, for tests without a GPU); ONEAPI_DEVICE_SELECTOR narrows the choice
inline int prefer_gpu(const sycl::device &d) { return d.is_gpu() ? 2 : d.is_cpu() ? 1 : 0; }

template <class Model>
class Solver {
public:
    using Fields = typename Model::Fields;
    using Params = typename Model::Params;

    explicit Solver(std::unique_ptr<Case<Model>> c)
        : q_(prefer_gpu, sycl::property::queue::in_order{}),
          fields_(q_, MaxBlks), mesh_(q_, MaxBlks),
          err_(q_, MaxBlks), code_(q_, MaxBlks), changed_(q_, 1),
          list_(q_, 4 * size_t(MaxBlks)), list2_(q_, 4 * size_t(MaxBlks)),
          case_(std::move(c)), err_slot_(MaxBlks + 1, 0) {}

    void run(int argc, char **argv);

private:
    sycl::queue q_;
    DeviceFields<Fields> fields_;
    DeviceMesh mesh_;
    DeviceArray<real> err_;
    DeviceArray<std::uint8_t> code_;    // criterion bits per leaf (store_criterion)
    DeviceArray<int> changed_;
    bool reflag_ = true;                // the bits are stale: flag on the host anyway
    DeviceArray<int> list_, list2_;     // block lists for the kernels
    Quadtree tree_;
    std::unique_ptr<Case<Model>> case_;
    Params p_{};
    std::vector<real> err_leaf_, err_slot_;
    HostFields<Fields::NC> host_;
    std::optional<VtkSeries> vtk_;     // empty: AMR_OUTPUT=none
    std::unique_ptr<FrameStream> render_;   // AMR_RENDER

    Fields F() const { return fields_.view(); }
    Mesh M() const { return mesh_.view(); }
    int nLeaf() const { return tree_.nLeaf(); }
    const BlockRef *put_refs(const std::vector<BlockRef> &r) {
        return reinterpret_cast<const BlockRef *>(list_.put(r));
    }

    //== Solver_Mod: Initialize_State
    void initialize_state() {
        const auto r = tree_.leaf_refs();
        case_->initial_state(q_, F(), put_refs(r), int(r.size()));
    }

    void ghost_cells_C() {
        ghost_C_pass<Model::NG>(q_, F(), M(), mesh_.leaves(), nLeaf(), false);
        if (mesh_.n_coarse() > 0) ghost_C_pass<Model::NG>(q_, F(), M(), mesh_.coarse(), mesh_.n_coarse(), true);
    }

    void ghost_cells_hg() {
        if (mesh_.n_hg1() > 0) ghost_hg_pass(q_, F(), M(), mesh_.hg1(), mesh_.n_hg1(), false);
        if (mesh_.n_coarse() > 0) ghost_hg_pass(q_, F(), M(), mesh_.coarse(), mesh_.n_coarse(), true);
    }

    //== criterion on the device, flagging on the host, field work on the device.
    //   The flagging is a function of the criterion bits of the leaves: if
    //   none changed since a flagging without mesh change, it is skipped
    bool adapt_mesh() {
        const real rt = case_->refine_tol, dt = case_->deref_tol;
        q_.fill(changed_.get(), 0, 1);
        Model::criterion(q_, F(), M(), p_, nLeaf(), rt, dt, err_.get(), code_.get(), changed_.get());
        int changed = 0;
        changed_.get(&changed, 1);
        if (!changed && !reflag_) return false;

        err_leaf_.resize(nLeaf());
        err_.get(err_leaf_.data(), err_leaf_.size());
        for (int m = 0; m < nLeaf(); ++m) err_slot_[tree_.leaves()[m]] = err_leaf_[m];

        const auto change = tree_.adapt(err_slot_, rt, dt, case_->max_island);
        reflag_ = change.has_value();       // a new mesh: the bits belong to the old leaves
        if (!change) return false;
        if (!change->restrict_pairs.empty())
            restrict_children(q_, F(), list_.put(change->restrict_pairs), int(change->restrict_pairs.size() / 2));
        if (!change->prolong_pairs.empty()) {
            clear_slots(q_, F(), list_.put(change->cleared), int(change->cleared.size()));
            prolong(q_, F(), list2_.put(change->prolong_pairs), int(change->prolong_pairs.size() / 2));
            Model::moments(q_, F(), p_, put_refs(change->new_blocks), int(change->new_blocks.size()));
        }
        mesh_.upload(tree_);
        return true;
    }
};

template <class Model>
void Solver<Model>::run(int argc, char **argv) {
    using clock = std::chrono::steady_clock;
    const auto t_start = clock::now();
    const char *s = std::getenv("AMR_GPU_SYNC");
    const bool sync = s && *s && *s != '0';
    const sycl::device dev = q_.get_device();
    std::printf("\n%s: %s\n", dev.is_gpu() ? "GPU" : "Device", dev.get_info<sycl::info::device::name>().c_str());
    std::printf("\nCase: %s\n", case_->name().c_str());
    p_ = case_->setup();
    const char *o = std::getenv("AMR_OUTPUT");
    const std::string writer = o && *o ? o : "vtk";
    if (writer == "vtk") vtk_.emplace();
    else if (writer != "none") throw std::runtime_error("AMR_OUTPUT: unknown writer " + writer + " (vtk, none)");
    std::printf("\nOutput writer: %s\n", writer.c_str());
    if (const char *r = std::getenv("AMR_RENDER"); r && *r) render_ = std::make_unique<FrameStream>(r);
    std::printf("\nBlock pool: %d blocks (%.0f MiB on the device)%s\n", MaxBlks,
                MaxBlks * DeviceFields<Fields>::mib_per_block(), periodicY ? ", periodic in x and y" : "");
    const auto vtk_fields = case_->vtk_fields(), render_fields = case_->render_fields();

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
    const int Tf = argc >= 2 ? std::atoi(argv[1]) : int(std::lround(case_->t_end / case_->time_per_step / 20)) * 20;
    const int nOutStep = argc >= 3 ? std::atoi(argv[2]) : std::max(Tf / 20, 1);
    if (Tf < 1 || nOutStep < 1) throw std::runtime_error("usage: amr_lbm_sycl [steps [output_interval]]");
    const int nDiagStep = std::max(Tf / 120, 1), nProgStep = std::max(Tf / 20, 1);

    std::printf("\nModel: %s\n", Model::name());
    case_->print_parameters(Tf);

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
                    Tdone, Tf, 100.0 * Tdone / Tf, double(Tdone * case_->time_per_step), nLeaf() * NX * NY,
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
                    Model::moments(q_, F(), p_, put_refs(r), int(r.size()));
                    moments_done = true;
                }
                host_.fetch(fields_, tree_.nblks(), (is_output && (vtk_ || render_)) || case_->diagnostics_need_all());
                const real time = real(Tn) * case_->time_per_step;
                if (is_output) {
                    if (vtk_) vtk_->step(time, tree_, host_, vtk_fields);
                    if (render_) render_->push(time, tree_, host_, render_fields);
                    case_->print_line(time, tree_, host_);
                }
                case_->diagnostics(time, tree_, host_);
            });
            if (Tn == Tf) break;
        }

        timed("collision", [&] { Model::collide(q_, F(), M(), p_, nLeaf(), !moments_done); });
        timed("Ghost_hg", [&] { ghost_cells_hg(); });
        timed("advection", [&] {
            stream(q_, F(), fields_.next(), M(), Levels(p_), nLeaf());
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

template <class Model>
void run(std::unique_ptr<Case<Model>> c, int argc, char **argv) {
    Solver<Model> solver(std::move(c));
    solver.run(argc, argv);
}

//== the case (and with it the model) of the namelist group in the AMR_CASE
//   file (default: rising bubble)
int main(int argc, char **argv) {
    try {
        const std::string g = namelist::group();
        if (g.empty() || g == "rising_bubble") run<JCP2016>(std::make_unique<RisingBubble>(), argc, argv);
        else if (g == "rayleigh_taylor") run<JCP2016>(std::make_unique<RayleighTaylor>(), argc, argv);
        else if (g == "parker") run<Parker<2>>(std::make_unique<ParkerBinary>(), argc, argv);
        else if (g == "ouzo") run<Parker<3>>(std::make_unique<OuzoCase>(), argc, argv);
        else throw std::runtime_error("AMR_CASE: unknown case &" + g + " (rising_bubble, rayleigh_taylor, parker, ouzo)");
    } catch (const sycl::exception &e) {
        std::fprintf(stderr, "SYCL error: %s\n", e.what());
        return 1;
    } catch (const std::exception &e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
