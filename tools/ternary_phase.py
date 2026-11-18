#!/usr/bin/env python3
"""Phase diagram of the ternary Margules mixture of the Parker model (gpu/).

The excess Gibbs energy (per mole, in units of RT) is the three-suffix
Margules / Wohl form used by the solver,

    g = sum_{i<j} x_i x_j (A_ji x_i + A_ij x_j) + C x_1 x_2 x_3,

and, as g is a homogeneous cubic of the mole fractions,

    ln gamma_i = dg/dx_i - 2 g.

The script computes the spinodal (the Hessian of g_mix = sum x ln x + g on the
simplex is not positive definite) and the binodal by continuation of the tie
lines from the immiscible binary edge (species 1-3) until the plait point. It
prints the coexisting compositions and, for a stock solution diluted with pure
species 3 along a straight mixing line, the stretch inside the spinodal.

    python3 tools/ternary_phase.py                       # the ouzo defaults
    python3 tools/ternary_phase.py --A13 3.5 --A12 0.3 --A23 0.6 -o phase.png
    python3 tools/ternary_phase.py --stock 0.25 0.3 0.45
"""
import argparse

import numpy as np
from scipy.optimize import fsolve


def margules(A, C):
    """ln gamma(x) for the matrix A (A[i][j] = A_ij) and the ternary constant C."""
    A = np.asarray(A, float)
    n = A.shape[0]

    def g_and_grad(x):
        g = 0.0
        dg = np.zeros(n)
        for i in range(n):
            for j in range(i + 1, n):
                g += x[i] * x[j] * (A[j, i] * x[i] + A[i, j] * x[j])
                dg[i] += x[j] * (2 * A[j, i] * x[i] + A[i, j] * x[j])
                dg[j] += x[i] * (A[j, i] * x[i] + 2 * A[i, j] * x[j])
        if n == 3:
            g += C * x[0] * x[1] * x[2]
            dg += C * np.array([x[1] * x[2], x[0] * x[2], x[0] * x[1]])
        return g, dg

    def ln_gamma(x):
        g, dg = g_and_grad(x)
        return dg - 2 * g

    def g_mix(x):
        x = np.clip(x, 1e-300, None)
        return float(np.sum(x * np.log(x)) + g_and_grad(x)[0])

    return ln_gamma, g_mix


def mu(ln_gamma, x):
    return np.log(np.clip(x, 1e-300, None)) + ln_gamma(x)


def hessian_det(g_mix, x1, x2, h=1e-5):
    """determinant of the Hessian of g_mix(x1, x2, 1-x1-x2) (central differences)"""
    f = lambda a, b: g_mix(np.array([a, b, 1 - a - b]))
    fxx = (f(x1 + h, x2) - 2 * f(x1, x2) + f(x1 - h, x2)) / h**2
    fyy = (f(x1, x2 + h) - 2 * f(x1, x2) + f(x1, x2 - h)) / h**2
    fxy = (f(x1 + h, x2 + h) - f(x1 + h, x2 - h) - f(x1 - h, x2 + h) + f(x1 - h, x2 - h)) / (4 * h * h)
    return fxx * fyy - fxy**2, fxx


def unstable(g_mix, x):
    if min(x) < 2e-5:
        return False
    d, fxx = hessian_det(g_mix, x[0], x[1])
    return d < 0 or fxx < 0


def binary_coexistence(ln_gamma):
    """species 1-3 edge: two phases with equal mu_1 and mu_3"""
    def res(v):
        a = np.array([v[0], 0, 1 - v[0]])
        b = np.array([v[1], 0, 1 - v[1]])
        d = mu(ln_gamma, a) - mu(ln_gamma, b)
        return [d[0], d[2]]
    v = fsolve(res, [0.02, 0.98], xtol=1e-14)
    return np.sort(v)


def tie_lines(ln_gamma, steps=400):
    """tie lines from the 1-3 edge, x_2 of the water-rich phase as the parameter"""
    lo, hi = binary_coexistence(ln_gamma)
    a0, b0 = np.array([lo, 0, 1 - lo]), np.array([hi, 0, 1 - hi])
    s0 = 1e-5                            # x2 = 0 exactly makes mu_2 singular: start just off the edge
    K = np.exp(ln_gamma(a0)[1] - ln_gamma(b0)[1])      # partition coefficient of species 2
    v = np.array([lo, hi, K * s0])       # x1 (water-rich), x1 (oil-rich), x2 (oil-rich)
    lines = []
    for s in np.linspace(s0, 0.6, steps):
        def res(w):
            a = np.array([w[0], s, 1 - w[0] - s])
            b = np.array([w[1], w[2], 1 - w[1] - w[2]])
            return mu(ln_gamma, a) - mu(ln_gamma, b)
        w = fsolve(res, v, xtol=1e-13)
        a = np.array([w[0], s, 1 - w[0] - s])
        b = np.array([w[1], w[2], 1 - w[1] - w[2]])
        if np.abs(res(w)).max() > 1e-9 or np.linalg.norm(a - b) < 2e-3 or min(a.min(), b.min()) < 0:
            break
        lines.append((a, b))
        v = w
    return lines


def one_phase(ln_gamma, x, n=300):
    """global stability: the tangent plane distance sum_i y_i (mu_i(y) - mu_i(x)) >= 0 on a grid of y"""
    g = np.linspace(1e-4, 1 - 1e-4, n)
    Y = np.array([[a, b, 1 - a - b] for a in g for b in g if 1 - a - b > 1e-4])
    MU = np.log(Y) + np.array([ln_gamma(y) for y in Y])
    return (Y * (MU - mu(ln_gamma, x))).sum(1).min() > -1e-10


def tern_xy(x):
    """triangle coordinates: species 3 (water) at (0,0), 1 (oil) at (1,0), 2 (ethanol) at the top"""
    x = np.atleast_2d(x)
    return x[:, 0] + 0.5 * x[:, 1], np.sqrt(3) / 2 * x[:, 1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--A13', type=float, default=4.0, help='Margules A_13 (oil-water)')
    ap.add_argument('--A31', type=float, default=None, help='A_31 (default: A_13)')
    ap.add_argument('--A12', type=float, default=0.3, help='A_12 (oil-ethanol)')
    ap.add_argument('--A21', type=float, default=None)
    ap.add_argument('--A23', type=float, default=0.6, help='A_23 (ethanol-water)')
    ap.add_argument('--A32', type=float, default=None)
    ap.add_argument('--C', type=float, default=0.0, help='ternary constant')
    ap.add_argument('--stock', type=float, nargs=3, default=[0.25, 0.3, 0.45], metavar=('X1', 'X2', 'X3'),
                    help='stock composition (oil, ethanol, water), diluted with pure species 3')
    ap.add_argument('--names', default='anethole,ethanol,water')
    ap.add_argument('-o', '--output', help='figure file')
    a = ap.parse_args()

    A = np.zeros((3, 3))
    A[0, 2], A[2, 0] = a.A13, a.A31 if a.A31 is not None else a.A13
    A[0, 1], A[1, 0] = a.A12, a.A21 if a.A21 is not None else a.A12
    A[1, 2], A[2, 1] = a.A23, a.A32 if a.A32 is not None else a.A23
    ln_gamma, g_mix = margules(A, a.C)
    names = a.names.split(',')

    lo, hi = binary_coexistence(ln_gamma)
    print(f'binary {names[0]}-{names[2]} coexistence x_{names[0]} = {lo:.5f}, {hi:.5f}')
    lines = tie_lines(ln_gamma)
    if lines:
        pa, pb = lines[-1]
        print(f'plait point near x = {0.5 * (pa + pb)} ({len(lines)} tie lines)')
    for a_, b_ in lines[:: max(len(lines) // 8, 1)]:
        print('  tie line', np.round(a_, 4), np.round(b_, 4))

    stock = np.array(a.stock) / np.sum(a.stock)
    state = 'unstable (inside the spinodal)' if unstable(g_mix, stock) else \
        'one phase' if one_phase(ln_gamma, stock) else 'metastable (between binodal and spinodal)'
    print(f'stock {np.round(stock, 4)}: {state}')
    s = np.linspace(0, 1, 2001)
    path = np.outer(s, stock) + np.outer(1 - s, [0, 0, 1])
    inside = np.array([unstable(g_mix, x) for x in path])
    if inside.any():
        print(f'dilution line inside the spinodal for stock fractions {s[inside].min():.3f} .. {s[inside].max():.3f}'
              f' (x_{names[0]} {path[inside][:, 0].min():.4f} .. {path[inside][:, 0].max():.4f})')
    else:
        print('dilution line does not enter the spinodal')

    if a.output:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        n = 300
        X1, X2 = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n))
        U = np.zeros_like(X1, dtype=bool)
        for i in range(n):
            for j in range(n):
                x = np.array([X1[i, j], X2[i, j], 1 - X1[i, j] - X2[i, j]])
                if x[2] > 0:
                    U[i, j] = unstable(g_mix, x)
        fig, ax = plt.subplots(figsize=(6, 5.4))
        tri = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
        ax.plot(*tern_xy(tri), 'k-', lw=1)
        pts = np.array([[X1[i, j], X2[i, j], 1 - X1[i, j] - X2[i, j]] for i in range(n) for j in range(n) if U[i, j]])
        if len(pts):
            ax.scatter(*tern_xy(pts), s=1, c='#f4b9a6', label='spinodal (unstable)')
        if lines:
            A_ = np.array([l[0] for l in lines]); B_ = np.array([l[1] for l in lines])
            bin_ = np.vstack([A_, B_[::-1]])
            ax.plot(*tern_xy(bin_), '-', c='#b8402a', lw=1.5, label='binodal')
            for a_, b_ in lines[:: max(len(lines) // 12, 1)]:
                ax.plot(*tern_xy(np.array([a_, b_])), '-', c='#777', lw=0.6)
        ax.plot(*tern_xy(np.array([stock, [0, 0, 1]])), '--', c='#2a62b8', lw=1.2, label='dilution with water')
        ax.plot(*tern_xy(stock), 'o', c='#2a62b8')
        for v, nm, (dx, dy, ha, va) in zip(np.eye(3), names, [(0.01, -0.03, 'left', 'top'), (0, 0.03, 'center', 'bottom'),
                                                             (-0.01, -0.03, 'right', 'top')]):
            x, y = tern_xy(v)
            ax.text(x[0] + dx, y[0] + dy, nm, ha=ha, va=va)
        ax.set_xlim(-0.12, 1.12); ax.set_ylim(-0.1, 0.96)
        ax.set_aspect('equal'); ax.axis('off')
        ax.legend(loc='upper right', fontsize=8, frameon=False)
        ax.set_title(f'Margules A13={A[0,2]:g}, A12={A[0,1]:g}, A23={A[1,2]:g}; '
                     f'stock {np.round(stock, 3).tolist()}', fontsize=9)
        fig.savefig(a.output, dpi=150, bbox_inches='tight')
        print('wrote', a.output)


if __name__ == '__main__':
    main()
