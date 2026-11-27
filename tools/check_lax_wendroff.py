#!/usr/bin/env python3
"""Check of the Lax-Wendroff streaming stencils of src/LBM_Common_Mod.f90:
Taylor analysis (SymPy) and von Neumann stability (NumPy).

    python3 tools/check_lax_wendroff.py

A population moving with velocity (ex, ey) is advanced by one time step dt on a
level with cell size h; c = dt/h is the CFL number passed to the kernels (0.5,
0.25, ... on the levels below the finest). The exact update is
f(x - c h ex, y - c h ey). The stencils are written here with the indices and
weights of the Fortran code:

    3-point          Streaming_LaxWendroff_AMR: upwind, central, downwind value
                     along the direction
    9-point tensor   Streaming_LaxWendroff9Pt_AMR, diagonal directions: the
                     1-D scheme in x, then in y
    9-point classic  the textbook 2-D Lax-Wendroff scheme with central
                     differences and a 4-point cross derivative; not in the
                     solver, for comparison

1. Taylor analysis: every stencil must reproduce, through h^2,

       f - dt (ex f_x + ey f_y) + dt^2/2 (ex^2 f_xx + 2 ex ey f_xy + ey^2 f_yy)

   The leading error (h^3) and the difference between the two 9-point
   stencils are printed.
2. von Neumann analysis: the largest amplification factor over all
   wavenumbers for the velocity (1,1). The 3-point and the tensor stencil
   stay at 1 up to c = 1; the classic one is unstable above c = 0.35, i.e. on
   the level below the finest (c = 0.5).

Exit status 0 if the Taylor check passes and the two stencils of the solver
are stable at the CFL numbers of the levels."""
import sys

import numpy as np
import sympy as sp

c, h, x, y = sp.symbols('c h x y')
cU = c*(1 + c)/2            # upwind weight
c0 = 1 - c**2               # central
cD = c*(1 - c)/2            # downwind (enters with a minus sign)
ORDER = 4                   # expand through h**ORDER

# derivatives of f at the cell as symbols: D[(a, b)] = d^(a+b) f / dx^a dy^b
D = {(a, b): sp.Symbol('f_' + 'x'*a + 'y'*b if a + b else 'f') for a in range(ORDER + 2) for b in range(ORDER + 2)}

def f(i, j):
    """Taylor series of f at the neighbour (X+i, Y+j), through h**ORDER."""
    return sum(D[(a, b)]*(i*h)**a*(j*h)**b/(sp.factorial(a)*sp.factorial(b))
               for a in range(ORDER + 1) for b in range(ORDER + 1 - a))

def characteristic(ex, ey):
    """Streaming_LaxWendroff_AMR: three points along the direction."""
    return cU*f(-ex, -ey) + c0*f(0, 0) - cD*f(ex, ey)

def tensor(ex, ey):
    """Streaming_LaxWendroff9Pt_AMR, diagonal directions: pass in x, then pass in y."""
    xpass = lambda j: cU*f(-ex, j) + c0*f(0, j) - cD*f(ex, j)
    return cU*xpass(-ey) + c0*xpass(0) - cD*xpass(ey)

def classic(ex, ey):
    """Textbook 2-D Lax-Wendroff: central differences, 4-point cross derivative."""
    fx  = (f(1, 0) - f(-1, 0))/(2*h);   fy  = (f(0, 1) - f(0, -1))/(2*h)
    fxx = (f(1, 0) - 2*f(0, 0) + f(-1, 0))/h**2
    fyy = (f(0, 1) - 2*f(0, 0) + f(0, -1))/h**2
    fxy = (f(1, 1) - f(1, -1) - f(-1, 1) + f(-1, -1))/(4*h**2)
    dt = c*h
    return f(0, 0) - dt*(ex*fx + ey*fy) + dt**2/2*(ex**2*fxx + 2*ex*ey*fxy + ey**2*fyy)

def exact(ex, ey):
    return f(-c*ex, -c*ey)

def lax_wendroff_terms(ex, ey):
    """f - dt (ex f_x + ey f_y) + dt^2/2 (ex^2 f_xx + 2 ex ey f_xy + ey^2 f_yy), dt = c h"""
    dt = c*h
    return (D[(0, 0)] - dt*(ex*D[(1, 0)] + ey*D[(0, 1)])
            + dt**2/2*(ex**2*D[(2, 0)] + 2*ex*ey*D[(1, 1)] + ey**2*D[(0, 2)]))

def order(expr, n):
    return sp.factor(sp.expand(expr).coeff(h, n))

ok = True
for i, (ex, ey) in {2: (1, 0), 3: (0, 1), 4: (-1, 0), 5: (0, -1), 6: (1, 1), 7: (-1, 1), 8: (-1, -1), 9: (1, -1)}.items():
    schemes = {'3-point': characteristic(ex, ey)}
    if ex and ey:
        schemes['9-point tensor'] = tensor(ex, ey)
        schemes['9-point classic'] = classic(ex, ey)
    for name, s in schemes.items():
        low = sp.expand(sp.series(sp.expand(s - lax_wendroff_terms(ex, ey)), h, 0, 3).removeO())
        ok = ok and low == 0
        print('direction %d (%2d,%2d) %-16s terms through h^2 equal Lax-Wendroff: %s' % (i, ex, ey, name, low == 0))
        if (ex, ey) == (1, 1):
            print('      error at h^3:', order(s - exact(ex, ey), 3))

print()
print('difference tensor - classic, direction (1,1):')
d = sp.expand(tensor(1, 1) - classic(1, 1))
for n in range(0, ORDER + 1):
    print('   h^%d: %s' % (n, order(d, n)))
print()
print('weights of the tensor stencil for direction (1,1), rows y = +1, 0, -1:')
w = {(i, j): sp.factor({-1: cU, 0: c0, 1: -cD}[i]*{-1: cU, 0: c0, 1: -cD}[j]) for i in (-1, 0, 1) for j in (-1, 0, 1)}
for j in (1, 0, -1):
    print('   ', [w[(i, j)] for i in (-1, 0, 1)])
print('   sum of the weights:', sp.simplify(sum(w.values())))
print('\nall stencils reproduce the Lax-Wendroff terms:', ok)


# ---------------------------------------------------------------------------
# von Neumann analysis: amplification factor g of a Fourier mode with phase
# angles (tx, ty) per cell, velocity (1,1)

def amplification(kind, cfl, n=721):
    wU, w0, wD = cfl*(1 + cfl)/2, 1 - cfl**2, cfl*(1 - cfl)/2
    t = np.linspace(-np.pi, np.pi, n)
    tx, ty = np.meshgrid(t, t)
    g1 = lambda th: wU*np.exp(-1j*th) + w0 - wD*np.exp(1j*th)       # the 1-D scheme
    if kind == '3-point':
        g = g1(tx + ty)
    elif kind == 'tensor':
        g = g1(tx)*g1(ty)
    else:                       # classic
        g = (1 - 1j*cfl*(np.sin(tx) + np.sin(ty)) + cfl**2*((np.cos(tx) - 1) + (np.cos(ty) - 1))
             - cfl**2*np.sin(tx)*np.sin(ty))
    return np.abs(g).max()

print()
print('largest amplification factor over all wavenumbers, velocity (1,1):')
print('    c      3-point   tensor   classic')
stable = True
for cfl in (0.125, 0.25, 0.35, 0.4, 0.5, 0.75, 1.0):
    g = [amplification(kind, cfl) for kind in ('3-point', 'tensor', 'classic')]
    stable = stable and max(g[:2]) < 1 + 1e-12
    print('    %-5g  %7.4f  %7.4f  %7.4f' % (cfl, *g))
print('\n3-point and tensor stencils stable up to c = 1:', stable)

sys.exit(0 if ok and stable else 1)
