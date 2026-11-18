#!/usr/bin/env python3
"""Plot the bubble front and spike tip of Rayleigh-Taylor runs (rt_fronts.dat).

Each file (written by the Rayleigh-Taylor case) has the columns t*, the
bubble front at x = 0 and the spike tip at x = L0/2, both relative to the
unperturbed interface in units of L0. With --paper the values read off
Fig. 6 of Fakhari et al., PRE 96, 053301 (2017) (At = 0.5, Re = 3000) are
shown as symbols (read off the figure, about +-0.01).

    python3 tools/plot_rt.py jcp/rt_fronts.dat:JCP2016 pre/rt_fronts.dat:PRE2017 --paper -o rt.png
"""
import argparse

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# PRE 2017, Fig. 6 (current model), read off the figure at 300 dpi
PAPER_T      = [0.0, 1.0, 1.5, 2.0, 2.5, 3.0]
PAPER_BUBBLE = [0.10, 0.29, 0.425, 0.555, 0.685, 0.835]
PAPER_SPIKE  = [-0.10, -0.36, -0.61, -0.86, -1.11, -1.43]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('files', nargs='+', help='rt_fronts.dat files, optionally as path:label')
    ap.add_argument('--paper', action='store_true', help='add the values of PRE 2017, Fig. 6')
    ap.add_argument('-o', '--output', default='rt_fronts.png')
    args = ap.parse_args()

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for spec in args.files:
        path, _, label = spec.partition(':')
        d = np.loadtxt(path, ndmin=2)
        line, = ax.plot(d[:, 0], d[:, 1], label=label or path)
        ax.plot(d[:, 0], d[:, 2], color=line.get_color())
    if args.paper:
        ax.plot(PAPER_T, PAPER_BUBBLE, 'ko', mfc='none', label='PRE 2017, Fig. 6')
        ax.plot(PAPER_T, PAPER_SPIKE, 'ko', mfc='none')
    ax.axhline(0, color='0.6', lw=0.8)
    ax.set_xlabel('t*')
    ax.set_ylabel('bubble front (> 0) and spike tip (< 0) / L0')
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output, dpi=120)
    print(args.output)


if __name__ == '__main__':
    main()
