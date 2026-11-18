/*
 * amrviz - real-time view of a block-structured AMR mesh and its fields
 *
 * A small C API (callable from C, C++ and Fortran through ISO_C_BINDING).
 * amrviz_open starts a render thread that owns an SDL2 window with an
 * OpenGL 3.3 core context. The solver hands over snapshots: a list of
 * square blocks of nx x ny cells, each with the values of every field.
 * The library copies them (as float) into a back buffer; amrviz_end_snapshot
 * publishes it to the render thread without waiting for it (triple
 * buffering), so the solver never blocks on drawing or vsync -- only while
 * the user has paused the run (amrviz_wait_if_paused).
 *
 * On the GPU every block is one instanced quad; the cell values of the
 * shown field are one flat float buffer read through a texture buffer
 * (cell k of block b at b*nx*ny + k), so the upload is proportional to the
 * number of leaf cells and the picture shows the true AMR layout.
 *
 * Indices: fields are numbered 0..nfields-1. Field nfields (the last entry
 * of the field menu) is the block level, taken from the block list.
 *
 * Controls (inspired by FluidX3D, https://github.com/ProjectPhysX/FluidX3D;
 * H shows them in the window):
 *   P, Space                pause / resume the solver (at its next snapshot)
 *   N                       one more snapshot while paused
 *   Z / Shift+Z, Tab        next / previous field
 *   C / Shift+C             next / previous colour map (viridis, magma,
 *                           coolwarm with a range symmetric about zero, gray)
 *   1 2 3 4                 field colouring, block outlines, cell outlines,
 *                           status panel on/off
 *   F                       freeze / unfreeze the colour range
 *   +, -, mouse wheel       zoom (the wheel about the cursor)
 *   I J K L, arrows, drag   pan
 *   R                       reset the view
 *   G                       print the view settings to stdout
 *   H                       help overlay
 *   Esc, close button       close the viewer (the solver runs on)
 *
 * The status panel shows the field and colour map with a colour bar, time,
 * step, number of blocks and cells, levels, the solver speed in steps/s and
 * the value under the cursor.
 *
 * Packed snapshots (GPU solvers): instead of amrviz_add_block and
 * amrviz_block_field_* per block, a snapshot can be handed over in one go:
 * amrviz_blocks_buffer(n) returns the block table to fill, and each field is
 * one float array of n*nx*ny cells (block b, cell (i,j) at b*nx*ny + j*nx + i),
 * packed by the caller, e.g. on the GPU. It reaches the viewer either
 *   - by a device-to-host copy straight into amrviz_field_buffer(field)
 *     (omp_target_memcpy, sycl::queue::memcpy, cudaMemcpy), or
 *   - through amrviz_field_packed(field, src) with any pointer the copy
 *     function registered with amrviz_set_copy can read (default: memcpy,
 *     i.e. host memory).
 * The library itself never touches device memory and links no GPU runtime.
 *
 * All functions except amrviz_open/amrviz_close are no-ops when the viewer
 * is not open (never opened, failed, or closed by the user), so the solver
 * needs to test amrviz_is_open() only to skip its own work.
 */
#ifndef AMRVIZ_H
#define AMRVIZ_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define AMRVIZ_MAX_FIELDS 16

/* Open the window: win_w x win_h pixels, blocks of nx x ny cells, a domain
 * of dom_w x dom_h (the unit of the block coordinates), nfields fields.
 * Returns 0 on success, -1 if the window or the GL context could not be
 * created (the reason is printed to stderr). */
int  amrviz_open(int win_w, int win_h, int nx, int ny,
                 double dom_w, double dom_h, int nfields);

/* Name of field 0..nfields-1, shown in the window title (NUL-terminated). */
void amrviz_set_field_name(int field, const char *name);

/* 1 while the window is open, 0 otherwise. */
int  amrviz_is_open(void);

/* Start a snapshot at the given time and step: clears the block list. */
void amrviz_begin_snapshot(double time, long step);

/* Append a block with lower-left corner (x0, y0) and cell size dx; returns
 * its slot (0-based) for amrviz_block_field_*, or -1 if the viewer is not
 * open or out of memory. */
int  amrviz_add_block(int level, double x0, double y0, double dx);

/* Copy the nx x ny cell values of a field of block `slot`. The values are
 * read column-major (Fortran order): cell (i,j), 0 <= i < nx, 0 <= j < ny,
 * is a[offset + i + j*ld]. ld and offset skip ghost layers, e.g. for
 * C(-1:nx+2,-1:ny+2): ld = nx+4, offset = 2 + 2*(nx+4). */
void amrviz_block_field_d(int slot, int field, const double *a, int ld, int offset);
void amrviz_block_field_f(int slot, int field, const float  *a, int ld, int offset);

/* --- packed snapshots (between amrviz_begin_snapshot and amrviz_end_snapshot) --- */

/* Set the number of blocks to n (replacing any added so far) and return the
 * block table to fill: 4 floats per block, (x0, y0, dx, level). Returns NULL
 * if the viewer is not open or out of memory. */
float *amrviz_blocks_buffer(int nblocks);

/* The host array of field `field` for the current blocks: nblocks*nx*ny
 * floats, block-major, cells in x-fastest order. Valid until
 * amrviz_end_snapshot. NULL if the viewer is not open or field is invalid. */
float *amrviz_field_buffer(int field);

/* Copy function for amrviz_field_packed: copy `bytes` bytes from src (the
 * caller's memory space) to the host pointer dst; return 0 on success. */
typedef int (*amrviz_copy_fn)(void *dst, const void *src, size_t bytes, void *user);

/* Register the copy function and its user pointer (NULL: memcpy). */
void amrviz_set_copy(amrviz_copy_fn copy, void *user);

/* Fill field `field` from a packed array of nblocks*nx*ny floats. */
void amrviz_field_packed(int field, const float *src);

/* Publish the snapshot to the render thread (never waits for drawing). */
void amrviz_end_snapshot(void);

/* Block while the user has paused the run (returns at once otherwise). */
void amrviz_wait_if_paused(void);

/* Close the viewer and free it. wait_for_user != 0: keep the window open
 * until the user closes it (e.g. at the end of the run). */
void amrviz_close(int wait_for_user);

#ifdef __cplusplus
}
#endif

#endif /* AMRVIZ_H */
