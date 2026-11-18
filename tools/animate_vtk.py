#!/usr/bin/env python3
"""Plot or animate the VTK output of the AMR-LBM rising-bubble solver.

Reads either the VTK XML time series written with AMR_OUTPUT=vtk
(Rising_Bubble_2D.pvd + AMR_Blocks.pvd, raw-appended .vtu files) or the
VTKHDF files written with AMR_OUTPUT=vtkhdf (Rising_Bubble_2D.vtkhdf +
AMR_Blocks.vtkhdf, needs h5py). Only numpy and matplotlib are required;
VTK/ParaView is not.

The cell field is painted onto a uniform image at the finest cell size, and
the AMR block outlines (coloured by refinement level) are drawn on top.

Examples
--------
  # animation of the phase field with the block structure (mp4 needs ffmpeg)
  python3 tools/animate_vtk.py Rising_Bubble_2D.pvd -o bubble.mp4

  # velocity magnitude, no block outlines, as a GIF
  python3 tools/animate_vtk.py Rising_Bubble_2D.pvd -f speed --no-blocks -o speed.gif

  # a single frame (the last one) as PNG
  python3 tools/animate_vtk.py Rising_Bubble_2D.vtkhdf --frame -1 -o last.png

  # interactive window
  python3 tools/animate_vtk.py Rising_Bubble_2D.pvd
"""
import argparse
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np

VTK_TYPES = {
    'Float32': np.float32, 'Float64': np.float64,
    'Int8': np.int8, 'UInt8': np.uint8, 'Int32': np.int32, 'UInt32': np.uint32,
    'Int64': np.int64, 'UInt64': np.uint64,
}

FIELDS = ('C', 'P', 'Ux', 'Uy', 'speed', 'Level', 'Block')


# ----------------------------------------------------------------------------
# readers: both return a list of frames, each a dict with
#   time, quads (ncell,4,2) corner coordinates, cell (dict of cell arrays)

def read_vtu(path):
    """Read an UnstructuredGrid .vtu file with raw appended data."""
    with open(path, 'rb') as f:
        raw = f.read()
    marker = b'<AppendedData encoding="raw">'
    k = raw.find(marker)
    if k < 0:
        raise ValueError(f'{path}: only raw appended .vtu files are supported')
    start = raw.index(b'_', k + len(marker)) + 1
    xml = raw[:k].decode() + '</VTKFile>'      # the binary section is cut off
    root = ET.fromstring(xml)
    hdr = np.dtype(VTK_TYPES[root.get('header_type', 'UInt32')])
    endian = '<' if root.get('byte_order', 'LittleEndian') == 'LittleEndian' else '>'

    def array(da):
        off = start + int(da.get('offset'))
        nbytes = int(np.frombuffer(raw, hdr.newbyteorder(endian), 1, off)[0])
        dt = np.dtype(VTK_TYPES[da.get('type')]).newbyteorder(endian)
        a = np.frombuffer(raw, dt, nbytes // dt.itemsize, off + hdr.itemsize)
        nc = int(da.get('NumberOfComponents', '1'))
        return a.reshape(-1, nc) if nc > 1 else a

    piece = root.find('.//Piece')
    pts = array(piece.find('Points/DataArray'))
    cells = {da.get('Name'): array(da) for da in piece.findall('Cells/DataArray')}
    cdata = {da.get('Name'): array(da) for da in piece.findall('CellData/DataArray')}
    conn = cells['connectivity'].reshape(-1, 4)            # all cells are quads
    return pts[conn][:, :, :2], cdata


def read_pvd(path):
    tree = ET.parse(path)
    base = os.path.dirname(path)
    frames = []
    for ds in tree.getroot().iter('DataSet'):
        quads, cdata = read_vtu(os.path.join(base, ds.get('file')))
        frames.append(dict(time=float(ds.get('timestep')), quads=quads, cell=cdata))
    return frames


def read_vtkhdf(path):
    import h5py
    frames = []
    with h5py.File(path, 'r') as f:
        r = f['VTKHDF']
        st = r['Steps']
        nsteps = int(st.attrs['NSteps'])
        for k in range(nsteps):
            po, npt = int(st['PointOffsets'][k]), int(r['NumberOfPoints'][k])
            co, nc = int(np.ravel(st['CellOffsets'][k])[0]), int(r['NumberOfCells'][k])
            io = int(np.ravel(st['ConnectivityIdOffsets'][k])[0])
            nid = int(r['NumberOfConnectivityIds'][k])
            pts = r['Points'][po:po + npt]
            conn = r['Connectivity'][io:io + nid].reshape(-1, 4)
            cdata = {}
            for name in r['CellData']:
                o = int(st['CellDataOffsets'][name][k])
                cdata[name] = r['CellData'][name][o:o + nc]
            frames.append(dict(time=float(st['Values'][k]),
                               quads=pts[conn][:, :, :2], cell=cdata))
    return frames


def read_series(path):
    if path.endswith('.pvd'):
        return read_pvd(path)
    if path.endswith('.vtkhdf') or path.endswith('.hdf'):
        return read_vtkhdf(path)
    if path.endswith('.vtu'):
        quads, cdata = read_vtu(path)
        return [dict(time=0.0, quads=quads, cell=cdata)]
    raise ValueError(f'unknown file type: {path}')


def blocks_path(cells_path):
    """AMR_Blocks.* next to Rising_Bubble_2D.* (same format)."""
    d, b = os.path.split(cells_path)
    ext = os.path.splitext(b)[1]
    p = os.path.join(d, 'AMR_Blocks' + ext)
    return p if os.path.exists(p) else None


# ----------------------------------------------------------------------------
# painting the AMR cells onto a uniform image

def cell_field(frame, name):
    c = frame['cell']
    if name in ('Ux', 'Uy', 'speed'):
        v = c['Velocity']
        return {'Ux': v[:, 0], 'Uy': v[:, 1], 'speed': np.hypot(v[:, 0], v[:, 1])}[name]
    if name not in c:
        raise KeyError(f'field {name!r} not in file (have {sorted(c)})')
    return np.asarray(c[name], dtype=float)


def rasterize(quads, values):
    """Paint axis-aligned square cells into an image with the finest cell size."""
    x0, y0 = quads[:, :, 0].min(1), quads[:, :, 1].min(1)
    size = quads[:, :, 0].max(1) - x0
    h = size.min()
    xmin, ymin = x0.min(), y0.min()
    xmax, ymax = (x0 + size).max(), (y0 + size).max()
    nxp, nyp = int(round((xmax - xmin) / h)), int(round((ymax - ymin) / h))
    img = np.full((nyp, nxp), np.nan)
    ix = np.rint((x0 - xmin) / h).astype(int)
    iy = np.rint((y0 - ymin) / h).astype(int)
    s = np.rint(size / h).astype(int)
    for sz in np.unique(s):                                   # one pass per level
        sel = s == sz
        r = np.arange(sz)
        rows = (iy[sel][:, None, None] + r[None, :, None])
        cols = (ix[sel][:, None, None] + r[None, None, :])
        img[rows, cols] = values[sel][:, None, None]
    return img, (xmin, xmax, ymin, ymax)


def block_segments(quads):
    """Outline segments of the blocks (closed quads)."""
    return np.concatenate([quads, quads[:, :1]], axis=1)


# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0],
                                 epilog=__doc__[__doc__.index('Examples'):],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input', nargs='?', default='Rising_Bubble_2D.pvd',
                    help='Rising_Bubble_2D.pvd, .vtkhdf or a single .vtu (default: %(default)s)')
    ap.add_argument('-f', '--field', default='C', choices=FIELDS, help='cell field (default: C)')
    ap.add_argument('-o', '--output', help='output .mp4/.gif/.png (default: show a window)')
    ap.add_argument('--frame', type=int, help='plot only this frame (negative: from the end)')
    ap.add_argument('--blocks', help='block file (default: AMR_Blocks.* next to input)')
    ap.add_argument('--no-blocks', action='store_true', help='do not draw the block outlines')
    ap.add_argument('--ylim', type=float, nargs=2, help='vertical range, e.g. 0 1.5')
    ap.add_argument('--clim', type=float, nargs=2, help='colour range (default: from all frames)')
    ap.add_argument('--cmap', default=None, help='matplotlib colormap')
    ap.add_argument('--fps', type=int, default=5, help='frames per second (default: 5)')
    ap.add_argument('--dpi', type=int, default=150)
    args = ap.parse_args()

    import matplotlib
    if args.output:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import animation
    from matplotlib.collections import LineCollection

    frames = read_series(args.input)
    blocks = None
    if not args.no_blocks:
        bp = args.blocks or blocks_path(args.input)
        if bp:
            blocks = read_series(bp)
            if len(blocks) != len(frames):
                print(f'warning: {len(blocks)} block frames vs {len(frames)} data frames',
                      file=sys.stderr)
                blocks = None
    if args.frame is not None:
        sel = [range(len(frames))[args.frame]]
    else:
        sel = list(range(len(frames)))
    print(f'{args.input}: {len(frames)} frames, plotting {len(sel)}; '
          f'{len(frames[sel[-1]]["quads"])} cells in the last one')

    images = [rasterize(frames[k]['quads'], cell_field(frames[k], args.field)) for k in sel]
    if args.clim:
        vmin, vmax = args.clim
    else:
        vmin = min(np.nanmin(im) for im, _ in images)
        vmax = max(np.nanmax(im) for im, _ in images)
        if args.field in ('Ux', 'Uy', 'P'):                    # symmetric for signed fields
            vmax = max(abs(vmin), abs(vmax)); vmin = -vmax
    cmap = args.cmap or {'C': 'RdBu_r', 'Level': 'viridis', 'Block': 'tab20'}.get(args.field, 'coolwarm' if args.field in ('Ux', 'Uy', 'P') else 'magma')

    img0, ext = images[0]
    aspect = (ext[3] - ext[2]) / (ext[1] - ext[0])
    if args.ylim:
        aspect = (args.ylim[1] - args.ylim[0]) / (ext[1] - ext[0])
    h = min(9.0, 3.0 * aspect)                      # axes height in inches
    fig, ax = plt.subplots(figsize=(h / aspect + 1.7, h + 0.9), layout='constrained')
    im = ax.imshow(img0, origin='lower', extent=ext, cmap=cmap, vmin=vmin, vmax=vmax,
                   interpolation='nearest')
    fig.colorbar(im, ax=ax, label=args.field, shrink=min(1.0, 4.0 / h + 0.2))
    lc = None
    if blocks is not None:
        lev = blocks[sel[0]]['cell']['Level']
        lc = LineCollection(block_segments(blocks[sel[0]]['quads']), linewidths=0.4,
                            cmap='Greys', norm=plt.Normalize(-3, max(lev.max(), 1)))
        lc.set_array(np.asarray(lev, float))
        ax.add_collection(lc)
    if args.ylim:
        ax.set_ylim(*args.ylim)
    ax.set_xlabel('x / L0'); ax.set_ylabel('y / L0')
    title = ax.set_title('')

    def draw(i):
        k = sel[i]
        img, e = images[i]
        im.set_data(img); im.set_extent(e)
        if lc is not None:
            lc.set_segments(block_segments(blocks[k]['quads']))
            lc.set_array(np.asarray(blocks[k]['cell']['Level'], float))
        nb = len(blocks[k]['quads']) if blocks is not None else None
        title.set_text(f't = {frames[k]["time"]:.3f}' + (f'   ({nb} blocks)' if nb else ''))
        return (im, title) + ((lc,) if lc is not None else ())

    draw(0)
    if len(sel) == 1:
        if args.output:
            fig.savefig(args.output, dpi=args.dpi)
            print('wrote', args.output)
        else:
            plt.show()
        return

    ani = animation.FuncAnimation(fig, draw, frames=len(sel), interval=1000 / args.fps, blit=False)
    if args.output:
        if args.output.endswith('.gif'):
            writer = animation.PillowWriter(fps=args.fps)
        else:
            writer = animation.FFMpegWriter(fps=args.fps, bitrate=4000)
        ani.save(args.output, writer=writer, dpi=args.dpi)
        print('wrote', args.output)
    else:
        plt.show()


if __name__ == '__main__':
    main()
