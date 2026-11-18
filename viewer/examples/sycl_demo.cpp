// amrviz with fields in GPU memory (SYCL)
//
// A synthetic AMR run: the host keeps a quadtree of 8x8-cell blocks refined
// around a moving ring; the GPU holds two fields in a solver-like pool
// (one block per slot, one ghost layer, float). Every frame a GPU kernel
// packs the interiors of the leaf blocks into one array per field, which
// reaches the viewer in the two ways amrviz supports:
//   field 0: q.memcpy straight into amrviz_field_buffer(0)
//   field 1: amrviz_field_packed(1, device_ptr) through a copy hook
//
// Build with viewer/examples/CMakeLists.txt (icpx), run ./sycl_demo [frames].
#include "amrviz.h"

#include <sycl/sycl.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

constexpr int NX = 8, NY = 8, NG = 1;                 // block cells, ghost layers
constexpr int SX = NX + 2 * NG, SY = NY + 2 * NG;     // stored block size
constexpr int L0 = 256, LMAX = 5;                     // domain L0 x L0, finest level
constexpr int NFIELDS = 2;

struct Leaf { float x0, y0, dx, level; };

// leaf blocks of a quadtree refined to LMAX where a block touches the ring
static void build_mesh(std::vector<Leaf> &leaves, float cx, float cy, float R)
{
    leaves.clear();
    std::vector<Leaf> stack{ { 0.0f, 0.0f, float(L0 / NX), 0.0f } };
    while (!stack.empty()) {
        Leaf b = stack.back();
        stack.pop_back();
        const float w = b.dx * NX;
        const float px = std::fmin(std::fmax(cx, b.x0), b.x0 + w);      // nearest point of the block
        const float py = std::fmin(std::fmax(cy, b.y0), b.y0 + w);
        const float dmin = std::hypot(px - cx, py - cy);
        const float dmax = std::hypot(std::fmax(std::fabs(b.x0 - cx), std::fabs(b.x0 + w - cx)),
                                      std::fmax(std::fabs(b.y0 - cy), std::fabs(b.y0 + w - cy)));
        const bool near = dmin <= R + 6.0f && dmax >= R - 6.0f;
        if (b.level < LMAX && (near || b.level < 1)) {
            const float h = 0.5f * b.dx, hw = h * NX;
            for (int k = 0; k < 4; ++k)
                stack.push_back({ b.x0 + (k & 1) * hw, b.y0 + (k >> 1) * hw, h, b.level + 1 });
        } else {
            leaves.push_back(b);
        }
    }
}

// copy hook for amrviz_field_packed: device (or USM) pointer -> host buffer
static int sycl_copy(void *dst, const void *src, size_t bytes, void *user)
{
    try {
        static_cast<sycl::queue *>(user)->memcpy(dst, src, bytes).wait();
        return 0;
    } catch (const sycl::exception &e) {
        std::fprintf(stderr, "sycl_copy: %s\n", e.what());
        return 1;
    }
}

int main(int argc, char **argv)
{
    const int frames = argc > 1 ? std::atoi(argv[1]) : 2000;
    sycl::queue q{ sycl::gpu_selector_v };
    std::printf("device: %s\n", q.get_device().get_info<sycl::info::device::name>().c_str());

    const int cap = 4096;                                  // block slots
    float *pool  = sycl::malloc_device<float>(size_t(NFIELDS) * cap * SX * SY, q);
    float *pack  = sycl::malloc_device<float>(size_t(NFIELDS) * cap * NX * NY, q);
    Leaf  *dmesh = sycl::malloc_device<Leaf>(cap, q);

    if (amrviz_open(800, 800, NX, NY, L0, L0, NFIELDS) != 0) return 1;
    amrviz_set_field_name(0, "C");
    amrviz_set_field_name(1, "wave");
    amrviz_set_copy(sycl_copy, &q);

    std::vector<Leaf> leaves;
    for (int frame = 0; frame < frames && amrviz_is_open(); ++frame) {
        const float t = 0.01f * frame;
        const float cx = 128.0f + 60.0f * std::cos(t), cy = 128.0f + 60.0f * std::sin(t), R = 40.0f;
        build_mesh(leaves, cx, cy, R);
        const int n = int(leaves.size());
        if (n > cap) { std::fprintf(stderr, "too many blocks\n"); break; }
        q.memcpy(dmesh, leaves.data(), n * sizeof(Leaf)).wait();

        // "solver": both fields in the pool, ghost layers included
        q.parallel_for(sycl::range<3>(n, SY, SX), [=](sycl::id<3> id) {
            const int b = id[0], j = id[1] - NG, i = id[2] - NG;
            const Leaf m = dmesh[b];
            const float x = m.x0 + (i + 0.5f) * m.dx, y = m.y0 + (j + 0.5f) * m.dx;
            const float r = sycl::hypot(x - cx, y - cy);
            const size_t k = (size_t(b) * SY + id[1]) * SX + id[2];
            pool[k] = 0.5f * (1.0f + sycl::tanh((r - R) / 3.0f));
            pool[size_t(cap) * SX * SY + k] = sycl::sin(0.1f * x - 3.0f * t) * sycl::cos(0.07f * y);
        });

        // pack the interiors: pack[f][b][j][i]
        q.parallel_for(sycl::range<3>(size_t(NFIELDS) * n, NY, NX), [=](sycl::id<3> id) {
            const int f = id[0] / n, b = id[0] % n, j = id[1], i = id[2];
            pack[((size_t(f) * n + b) * NY + j) * NX + i] =
                pool[((size_t(f) * cap + b) * SY + j + NG) * SX + i + NG];
        }).wait();

        amrviz_begin_snapshot(t, frame);
        float *meta = amrviz_blocks_buffer(n);
        for (int b = 0; b < n; ++b) {
            meta[4 * b + 0] = leaves[b].x0;
            meta[4 * b + 1] = leaves[b].y0;
            meta[4 * b + 2] = leaves[b].dx;
            meta[4 * b + 3] = leaves[b].level;
        }
        q.memcpy(amrviz_field_buffer(0), pack, size_t(n) * NX * NY * sizeof(float)).wait();
        amrviz_field_packed(1, pack + size_t(n) * NX * NY);
        amrviz_end_snapshot();
        amrviz_wait_if_paused();
    }

    amrviz_close(1);
    sycl::free(pool, q);
    sycl::free(pack, q);
    sycl::free(dmesh, q);
    return 0;
}
