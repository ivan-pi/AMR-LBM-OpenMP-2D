#!/usr/bin/env python3
"""Plot the rise of the bubble from one or more bubble_rise.dat files.

Each file (written by the rising-bubble case) has the columns
t*, y_c (centroid height, lattice units), u_r (rise velocity) and
Re_r = rho_h u_r D / mu_h.

    python3 tools/plot_rise.py run1/bubble_rise.dat:JCP run2/bubble_rise.dat:PRE -o rise.png
    python3 tools/plot_rise.py a/bubble_rise.dat --quantity y_c
"""
import argparse

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

COLUMNS = {'y_c': 1, 'u_r': 2, 'Re_r': 3}
LABELS = {'y_c': r'$y_c$ (lattice units)', 'u_r': r'$u_r$ (lattice units)', 'Re_r': r'Re$_r$'}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('files', nargs='+', help='bubble_rise.dat files, optionally as path:label')
    ap.add_argument('-q', '--quantity', choices=sorted(COLUMNS), default='Re_r')
    ap.add_argument('-o', '--output', default='rise.png')
    ap.add_argument('--title', default='')
    args = ap.parse_args()

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for spec in args.files:
        path, _, label = spec.partition(':')
        d = np.loadtxt(path, ndmin=2)
        ax.plot(d[:, 0], d[:, COLUMNS[args.quantity]], label=label or path)
    ax.set_xlabel('t*')
    ax.set_ylabel(LABELS[args.quantity])
    if args.title:
        ax.set_title(args.title)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output, dpi=120)
    print(args.output)


if __name__ == '__main__':
    main()
