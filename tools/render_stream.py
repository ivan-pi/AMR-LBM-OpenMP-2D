#!/usr/bin/env python3
"""Render the frames streamed by the GPU solver into a movie, as they arrive.

amr_lbm_sycl writes a frame to this script's standard input at every output
step when started with AMR_RENDER set to this command (in-situ rendering:
no field files on disk, drawing overlaps the run):

  AMR_OUTPUT=none AMR_RENDER="python3 tools/render_stream.py -o rt.mp4" \\
      ./build-gpu/amr_lbm_sycl 200000 1000

Frame layout (little endian): "AMRF", int32 version (1), float64 time,
int32 L0, nBaseX, nBaseY, block size n, number of leaves m; m x int32
(x0, y0, dx) in finest-level lattice units; m x n x n float32 C (x fastest).
Each frame is painted onto a uniform image at the finest cell size and
written to the movie at once, so memory does not grow with the run.
Needs numpy, matplotlib and ffmpeg.
"""
import argparse
import struct
import sys

import numpy as np


def read_exact(f, n):
    b = f.read(n)
    if len(b) < n:
        if b:
            raise EOFError('truncated frame')
        return None
    return b


def read_frame(f):
    magic = read_exact(f, 4)
    if magic is None:
        return None
    if magic != b'AMRF':
        raise ValueError(f'bad frame marker {magic!r}')
    (version,) = struct.unpack('<i', read_exact(f, 4))
    if version != 1:
        raise ValueError(f'frame version {version}')
    (time,) = struct.unpack('<d', read_exact(f, 8))
    L0, nbx, nby, n, m = struct.unpack('<5i', read_exact(f, 20))
    blocks = np.frombuffer(read_exact(f, 12 * m), '<i4').reshape(m, 3)
    C = np.frombuffer(read_exact(f, 4 * m * n * n), '<f4').reshape(m, n, n)     # (block, y, x)
    return dict(time=time, L0=L0, width=nbx * L0, height=nby * L0, n=n, blocks=blocks, C=C)


def rasterize(fr):
    """The leaf cells painted onto the finest-level grid (rows = y)."""
    img = np.empty((fr['height'], fr['width']), np.float32)
    n, b = fr['n'], fr['blocks']
    for dx in np.unique(b[:, 2]):
        sel = b[:, 2] == dx
        v = fr['C'][sel]
        if dx > 1:
            v = v.repeat(dx, axis=1).repeat(dx, axis=2)
        s = n * dx
        rows = b[sel, 1][:, None, None] + np.arange(s)[None, :, None]
        cols = b[sel, 0][:, None, None] + np.arange(s)[None, None, :]
        img[rows, cols] = v
    return img


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument('-o', '--output', required=True, help='movie file (.mp4)')
    ap.add_argument('--fps', type=int, default=20)
    ap.add_argument('--height', type=float, default=9.6, help='figure height in inches (default 9.6)')
    ap.add_argument('--dpi', type=int, default=200)
    ap.add_argument('--clim', type=float, nargs=2, default=(0.0, 1.0))
    ap.add_argument('--cmap', default='RdBu_r')
    ap.add_argument('--label', default='C', help='colour bar label')
    ap.add_argument('--title', default='', help='text before the time in the title')
    args = ap.parse_args()

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import animation

    src = sys.stdin.buffer
    fr = read_frame(src)
    if fr is None:
        sys.exit('render_stream: no frames')

    aspect = fr['height'] / fr['width']
    axh = args.height - 0.9
    w_in = axh / aspect + 1.6
    # even pixel sizes (H.264)
    W = 2 * round(w_in * args.dpi / 2) / args.dpi
    H = 2 * round(args.height * args.dpi / 2) / args.dpi
    fig, ax = plt.subplots(figsize=(W, H), layout='constrained')
    ext = (0, fr['width'] / fr['L0'], 0, fr['height'] / fr['L0'])
    im = ax.imshow(rasterize(fr), origin='lower', extent=ext, cmap=args.cmap,
                   vmin=args.clim[0], vmax=args.clim[1], interpolation='antialiased')
    fig.colorbar(im, ax=ax, label=args.label, shrink=0.5)
    ax.set_xlabel('x / L0')
    ax.set_ylabel('y / L0')
    title = ax.set_title('')

    writer = animation.FFMpegWriter(fps=args.fps, bitrate=8000)
    count = 0
    with writer.saving(fig, args.output, args.dpi):
        while fr is not None:
            im.set_data(rasterize(fr))
            title.set_text(f'{args.title}t* = {fr["time"]:.3f}   ({len(fr["blocks"])} blocks)')
            writer.grab_frame()
            count += 1
            print(f'render_stream: frame {count}, t* = {fr["time"]:.3f}', file=sys.stderr, flush=True)
            fr = read_frame(src)
    print(f'render_stream: wrote {args.output} ({count} frames)', file=sys.stderr)


if __name__ == '__main__':
    main()
