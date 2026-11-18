#!/usr/bin/env python3
"""Render the frames streamed by the GPU solver into a movie, as they arrive.

amr_lbm_sycl writes a frame to this script's standard input at every output
step when started with AMR_RENDER set to this command (in-situ rendering:
no field files on disk, drawing overlaps the run):

  AMR_OUTPUT=none AMR_RENDER="python3 tools/render_stream.py -o rt.mp4" \\
      ./build-gpu/amr_lbm_sycl 200000 1000

One column per field (--fields, default rho,speed,level), each painted onto a
uniform image at the finest cell size; "level" is the refinement level of
each block. Each frame is written to the movie at once, so memory does not
grow with the run. The frames carry no text but the time.

Frame layout (little endian): "AMRF", int32 version (2), float64 time,
int32 L0, nBaseX, nBaseY, block size n, number of leaves m, number of fields
k; k names of 16 bytes; m x int32 (x0, y0, dx) in finest-level lattice
units; per field m x n x n float32 (x fastest).
Needs numpy, matplotlib and ffmpeg.
"""
import argparse
import struct
import sys

import numpy as np

LABELS = {'C': 'C', 'rho': 'density', 'speed': '|u|', 'level': 'level'}
CMAPS = {'C': 'RdBu_r', 'rho': 'RdBu_r', 'speed': 'magma', 'level': 'viridis'}


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
    if version != 2:
        raise ValueError(f'frame version {version} (expected 2)')
    (time,) = struct.unpack('<d', read_exact(f, 8))
    L0, nbx, nby, n, m, k = struct.unpack('<6i', read_exact(f, 24))
    names = [read_exact(f, 16).decode().strip() for _ in range(k)]
    blocks = np.frombuffer(read_exact(f, 12 * m), '<i4').reshape(m, 3)
    fields = {name: np.frombuffer(read_exact(f, 4 * m * n * n), '<f4').reshape(m, n, n) for name in names}
    return dict(time=time, L0=L0, width=nbx * L0, height=nby * L0, n=n, blocks=blocks, fields=fields)


def rasterize(fr, name):
    """A cell field (or the block level) painted onto the finest-level grid (rows = y)."""
    n, b = fr['n'], fr['blocks']
    if name == 'level':
        top = int(np.log2(fr['L0'] // n))                     # finest level: dx = 1
        vals = (top - np.log2(b[:, 2])).astype(np.float32)[:, None, None] * np.ones((1, n, n), np.float32)
    else:
        vals = fr['fields'][name]
    img = np.empty((fr['height'], fr['width']), np.float32)
    for dx in np.unique(b[:, 2]):
        sel = b[:, 2] == dx
        v = vals[sel]
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
    ap.add_argument('--fields', default='rho,speed,level', help='comma-separated: C, rho, speed, level')
    ap.add_argument('--clim', action='append', default=[], metavar='NAME=LO,HI',
                    help='colour range of a field, e.g. speed=0,0.05 (default: from the first frame; '
                         'speed: 0 to 2 sqrt(g L0) is a good choice)')
    ap.add_argument('--fps', type=int, default=20)
    ap.add_argument('--height', type=float, default=10.8, help='figure height in inches (default 10.8)')
    ap.add_argument('--dpi', type=int, default=200)
    args = ap.parse_args()
    names = args.fields.split(',')
    clims = {k: tuple(map(float, v.split(','))) for k, v in (c.split('=') for c in args.clim)}

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import animation, colors

    src = sys.stdin.buffer
    fr = read_frame(src)
    if fr is None:
        sys.exit('render_stream: no frames')
    for name in names:
        if name != 'level' and name not in fr['fields']:
            sys.exit(f'render_stream: no field {name!r} (have {sorted(fr["fields"])} and level)')

    aspect = fr['height'] / fr['width']
    axh = args.height - 1.5                                    # axes height; colour bars below
    w_in = len(names) * (axh / aspect + 0.45) + 0.6
    W = 2 * round(w_in * args.dpi / 2) / args.dpi              # even pixel sizes (H.264)
    H = 2 * round(args.height * args.dpi / 2) / args.dpi
    fig, axes = plt.subplots(1, len(names), figsize=(W, H), sharey=True, layout='constrained')
    axes = np.atleast_1d(axes)
    ext = (0, fr['width'] / fr['L0'], 0, fr['height'] / fr['L0'])
    images = []
    for ax, name in zip(axes, names):
        img = rasterize(fr, name)
        kw = dict(cmap=CMAPS.get(name, 'viridis'))
        if name == 'level':                                    # from the coarsest level of the first frame
            top = int(np.log2(fr['L0'] // fr['n']))
            low = int(np.nanmin(img))
            kw['norm'] = colors.BoundaryNorm(np.arange(low - 0.5, top + 1), plt.get_cmap(kw['cmap']).N)
        else:
            lo, hi = clims.get(name, (float(np.nanmin(img)), float(np.nanmax(img))))
            if hi <= lo:
                hi = lo + 1
            kw.update(vmin=lo, vmax=hi)
        im = ax.imshow(img, origin='lower', extent=ext, interpolation='antialiased', **kw)
        cb = fig.colorbar(im, ax=ax, location='bottom', shrink=0.9, aspect=12, label=LABELS.get(name, name))
        if name == 'level':
            cb.set_ticks(range(low, top + 1))
        ax.set_xlabel('x / L0')
        images.append((im, name))
    axes[0].set_ylabel('y / L0')
    stamp = fig.suptitle('', fontsize=9)

    writer = animation.FFMpegWriter(fps=args.fps, bitrate=8000)
    count = 0
    with writer.saving(fig, args.output, args.dpi):
        while fr is not None:
            for im, name in images:
                im.set_data(rasterize(fr, name))
            stamp.set_text(f't* = {fr["time"]:.2f}')
            writer.grab_frame()
            count += 1
            print(f'render_stream: frame {count}, t* = {fr["time"]:.3f}', file=sys.stderr, flush=True)
            fr = read_frame(src)
    print(f'render_stream: wrote {args.output} ({count} frames)', file=sys.stderr)


if __name__ == '__main__':
    main()
