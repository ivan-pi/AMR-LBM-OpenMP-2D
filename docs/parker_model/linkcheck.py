"""Check the exact discrete steady balance of the LB on links:
   1/2 [(c u')_x + (c u')_{x+1}] = (D/cs2) (Pi^eq_{x+1} - Pi^eq_x),  Pi^eq = c (cs2 + kappa u'^2)"""
import sys
import numpy as np
sys.argv = sys.argv[:1]
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import parker_check as P

for nu in (1 / 6, 0.1):
    D, kap = 0.1, 1.0
    xa, xb = P.coexistence()
    X = np.arange(80) + 0.5
    prof = 0.5 * (np.tanh((X - 20) / 2) - np.tanh((X - 60) / 2))
    s = P.Parker((xa + (xb - xa) * prof)[None, :].repeat(2, 0), D, nu, kap)
    for t in range(40000):
        s.step()
    ci, c, ux, uy, px, py, x1, lg = s.macro()
    for i in range(2):
        cu = (ci[i] * px[i])[0]
        Pi = (ci[i] * (P.cs2 + kap * px[i]**2))[0]
        lhs = 0.5 * (cu + np.roll(cu, -1))
        rhs = D / P.cs2 * (np.roll(Pi, -1) - Pi)
        rhs_wide = D / P.cs2 * 0.5 * (np.roll(Pi, -1) - np.roll(Pi, 1))
        print(f'nu={nu:.4f} (omega_nu={1/(nu/P.cs2+.5):.3f}) species {i+1}: max|c u\'|={np.abs(cu).max():.3e}'
              f'  link residual={np.abs(lhs - rhs).max():.2e}   nodal central-difference residual={np.abs(cu - rhs_wide).max():.2e}')
