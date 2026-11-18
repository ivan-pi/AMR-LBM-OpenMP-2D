"""
Symbolic analysis of Parker's multicomponent LBM (Parker 2008; Parker &
Jovanovic 2009) with the scaled non-conserved equilibrium moments (kappa)
of the central-moment draft (latex/moments.tex, scaled_raw_moments.m).

Goal: the surface tension in terms of the model moments and parameters.

  1. Parker's MRT basis and equilibria (TwoDSampleCode.f90) -> raw moments;
     identify the scaling k4..k9 they imply
  2. the mixture's second moments: conserved part + kappa sum_i c_i u'_i u'_i
  3. second-order (Chapman-Enskog / Taylor) analysis of a steady flat
     interface at rest: species flux J_i and mixture momentum flux Pi
  4. sigma = int (Pi_nn - Pi_tt) dn in terms of kappa, D and the profile
  5. the diffusion velocity of the code, the dead grad(c) term and the
     Gibbs-Duhem relation; the binary closed form

Run: python3 parker_surface_tension.py   (needs sympy)
"""
import sympy as sp

R = sp.Rational
cs2 = R(1, 3)

# D2Q9 in Parker's ordering (TwoDSampleCode.f90: 1 rest, 2 E, 3 N, 4 W, 5 S, 6 NE, 7 NW, 8 SW, 9 SE)
cx = [0, 1, 0, -1, 0, 1, -1, -1, 1]
cy = [0, 0, 1, 0, -1, 1, 1, -1, -1]
Q = 9

# Parker's (Lallemand-Luo) basis: rho, e, eps, jx, qx, jy, qy, pxx, pxy
Mp = sp.Matrix([
    [1] * Q,
    [-4 + 3 * (cx[k]**2 + cy[k]**2) for k in range(Q)],
    [4 - R(21, 2) * (cx[k]**2 + cy[k]**2) + R(9, 2) * (cx[k]**2 + cy[k]**2)**2 for k in range(Q)],
    cx,
    [(-5 + 3 * (cx[k]**2 + cy[k]**2)) * cx[k] for k in range(Q)],
    cy,
    [(-5 + 3 * (cx[k]**2 + cy[k]**2)) * cy[k] for k in range(Q)],
    [cx[k]**2 - cy[k]**2 for k in range(Q)],
    [cx[k] * cy[k] for k in range(Q)],
])
# raw basis: 1, x, y, x^2, y^2, xy, x^2 y, x y^2, x^2 y^2
mono = [(0, 0), (1, 0), (0, 1), (2, 0), (0, 2), (1, 1), (2, 1), (1, 2), (2, 2)]
Mr = sp.Matrix([[cx[k]**a * cy[k]**b for k in range(Q)] for (a, b) in mono])
names_r = ['m00', 'm10', 'm01', 'm20', 'm02', 'm11', 'm21', 'm12', 'm22']

print('=' * 78)
print('1. Parker basis vs the code')
print('=' * 78)
# rows of the code (n2, n3, n5, n7, n8, n9) and its inverse (link formulas, /36)
code_rows = {1: [-4, -1, -1, -1, -1, 2, 2, 2, 2], 2: [4, -2, -2, -2, -2, 1, 1, 1, 1],
             4: [0, -2, 0, 2, 0, 1, -1, -1, 1], 6: [0, 0, -2, 0, 2, 1, 1, -1, -1],
             7: [0, 1, -1, 1, -1, 0, 0, 0, 0], 8: [0, 0, 0, 0, 0, 1, -1, 1, -1]}
assert all(list(Mp.row(r)) == v for r, v in code_rows.items())
Minv_code = sp.Matrix([
    [4, -4, 4, 0, 0, 0, 0, 0, 0],
    [4, -1, -2, 6, -6, 0, 0, 9, 0],
    [4, -1, -2, 0, 0, 6, -6, -9, 0],
    [4, -1, -2, -6, 6, 0, 0, 9, 0],
    [4, -1, -2, 0, 0, -6, 6, -9, 0],
    [4, 2, 1, 6, 3, 6, 3, 0, 9],
    [4, 2, 1, -6, -3, 6, 3, 0, -9],
    [4, 2, 1, -6, -3, -6, -3, 0, 9],
    [4, 2, 1, 6, 3, -6, -3, 0, -9]]) / 36
assert sp.simplify(Minv_code - Mp.inv()) == sp.zeros(Q, Q)
print('M (rows n1..n9) and M^-1 (link formulas) of TwoDSampleCode.f90: consistent')

# ---------------------------------------------------------------------------
c = sp.Symbol('c_i', positive=True)
ux, uy = sp.symbols('u_x u_y', real=True)          # bulk (mixture) velocity
px, py = sp.symbols("u'_x u'_y", real=True)        # diffusion velocity of species i
kap = sp.Symbol('kappa', real=True)

U2 = ux**2 + uy**2 + 2 * (ux * px + uy * py) + kap * (px**2 + py**2)
eqP = sp.Matrix([c,
                 c * (-2 + 3 * U2),
                 c * (1 - 3 * U2),
                 c * (ux + px), -c * (ux + px),
                 c * (uy + py), -c * (uy + py),
                 c * (ux**2 - uy**2 + 2 * px * ux - 2 * py * uy + kap * (px**2 - py**2)),
                 c * (ux * uy + px * uy + py * ux + kap * px * py)])
eqR = sp.expand(Mr * Mp.inv() * eqP)

print('\nParker equilibrium as raw moments:')
for n, e in zip(names_r, eqR):
    print(f'  {n} = {sp.factor(e)}')

# scaled raw moments of the draft (scaled_raw_moments.m)
k = sp.symbols('k4:10', real=True)
k4, k5, k6, k7, k8, k9 = k
vx, vy = ux + px, uy + py
srm = sp.Matrix([
    c, c * vx, c * vy,
    c * (cs2 + ux**2 + 2 * ux * px + k4 * px**2),
    c * (cs2 + uy**2 + 2 * uy * py + k5 * py**2),
    c * (ux * uy + ux * py + px * uy + k6 * px * py),
    c * ((cs2 + ux**2) * vy + 2 * ux * uy * px + k7 * (2 * ux * px * py + px**2 * vy)),
    c * ((cs2 + uy**2) * vx + 2 * uy * ux * py + k8 * (2 * uy * py * px + py**2 * vx)),
    c * ((cs2 + ux**2) * (cs2 + uy**2) + 2 * cs2 * (ux * px + uy * py) + 2 * ux**2 * uy * py
         + 2 * uy**2 * ux * px + k9 * cs2 * (px**2 + py**2) + k9 * ux**2 * py**2 + k9 * uy**2 * px**2
         + k9 * (2 * ux * px + px**2) * (2 * uy * py + py**2))])
d2 = sp.expand(eqR[3:6, 0] - srm[3:6, 0].subs({k4: kap, k5: kap, k6: kap}))
print('\nsecond-order rows equal the scaled raw moments with k4 = k5 = k6 = kappa:', d2 == sp.zeros(3, 1))
print('third/fourth-order rows: Parker keeps only the terms of the second-order')
print('(Lallemand-Luo) equilibrium, i.e. no k7, k8, k9 scaling exists; e.g.')
print('  m21 - srm21(k7=0) =', sp.factor(sp.expand(eqR[6] - srm[6].subs(k7, 0))))

print('\n' + '=' * 78)
print('2. Mixture second moments (two species, sum_i c_i u\'_i = 0)')
print('=' * 78)
c1, c2 = sp.symbols('c_1 c_2', positive=True)
p1x, p1y, p2x, p2y = sp.symbols("u'_1x u'_1y u'_2x u'_2y", real=True)
def eq_species(ci, pix, piy):
    return srm.subs({c: ci, px: pix, py: piy, k4: kap, k5: kap, k6: kap})
S = eq_species(c1, p1x, p1y) + eq_species(c2, p2x, p2y)
S = S.subs({p2x: -c1 * p1x / c2, p2y: -c1 * p1y / c2})
C = c1 + c2
conserved = sp.Matrix([C * (cs2 + ux**2), C * (cs2 + uy**2), C * ux * uy])
extra = sp.simplify(sp.expand(S[3:6, 0] - conserved))
print('sum_i Pi^eq_i - (c cs^2 I + c u u) = ')
for n, e in zip(['xx', 'yy', 'xy'], extra):
    print(f'  {n}: {sp.factor(e)}')
check = sp.simplify(extra - kap * sp.Matrix([c1 * p1x**2 + c2 * (c1 * p1x / c2)**2,
                                               c1 * p1y**2 + c2 * (c1 * p1y / c2)**2,
                                               c1 * p1x * p1y + c2 * (c1 / c2)**2 * p1x * p1y]))
print('  = kappa * sum_i c_i u\'_i u\'_i :', check == sp.zeros(3, 1))

print('\n' + '=' * 78)
print('3. Second-order analysis of a steady flat interface at rest (normal x)')
print('=' * 78)
# Taylor expansion of the steady LB equation f_k(x) = f*_k(x - c_k):
#   f1 = -Omega^-1 T feq,  Omega f2 = -T (I - Omega/2) f1,  T = c_k d/dx.
# Fluxes to second order: moments of feq + (I - Omega/2) f1.
# Diffusive ordering: u' = O(d/dx), so the mass row of T feq, d(c u')/dx, is
# second order and the mass row of f1 is zero.
x = sp.Symbol('x', real=True)
w2, w3, wD, wV = sp.symbols('omega_e omega_eps omega_D omega_nu', positive=True)
K = sp.diag(1, w2, w3, wD, wD, wD, wD, wV, wV)     # mass row: placeholder, its f1 row is 0
A = Mp * sp.diag(*cx) * Mp.inv()                     # flux matrix d/dx in Parker moments

def fluxes(ci, pi):
    """raw moments of feq + (I-K/2) f1 for species with c(x), u'_x(x), u=0, u'_y=0"""
    meqP = eqP.subs({ux: 0, uy: 0, py: 0}).subs({c: ci, px: pi})
    Tm = A * meqP.diff(x)
    Tm[0] = 0                       # O(d^2): goes to the next order
    m1 = -K.inv() * Tm
    m1[0] = 0
    return sp.expand(Mr * Mp.inv() * (meqP + (sp.eye(Q) - K / 2) * m1))

c1f, c2f = sp.Function('c_1')(x), sp.Function('c_2')(x)
p1f = sp.Function("u'_1")(x)
p2f = -c1f * p1f / c2f              # sum_i c_i u'_i = 0 (mass-average reference)
F1, F2 = fluxes(c1f, p1f), fluxes(c2f, p2f)
D = sp.Symbol('D', positive=True)
J1 = sp.simplify(F1[1])
print('species flux J_1 =', J1)
J1_expected = c1f * p1f - (1 / wD - R(1, 2)) * sp.diff(c1f * (cs2 + kap * p1f**2), x)
print('  = c_1 u\'_1 - (1/omega_D - 1/2) d/dx[ c_1 (cs^2 + kappa u\'_1^2) ] :',
      sp.simplify(J1 - J1_expected) == 0)
print('  -> D = cs^2 (1/omega_D - 1/2) = diffco of the code')

Pi = (F1 + F2)
Pxx, Pyy, Pxy = [sp.simplify(Pi[r]) for r in (3, 4, 5)]
print('\nmixture momentum flux (second order, steady, u = 0):')
print('  Pi_xx =', sp.factor(Pxx))
print('  Pi_yy =', sp.factor(Pyy))
print('  Pi_xy =', Pxy)
N = sp.simplify(Pxx - Pyy - kap * (c1f * p1f**2 + c2f * p2f**2))
print('  Pi_xx - Pi_yy - kappa sum_i c_i u\'_i^2 =', N)
print('  (no omega_nu, omega_e, omega_eps: the viscous/bulk terms vanish at rest)')
print('  mechanical equilibrium d/dx Pi_xx = 0:  c cs^2 + kappa sum_i c_i u\'_i^2 = p_0')

print('\n' + '=' * 78)
print('4. Surface tension')
print('=' * 78)
print('sigma = int (Pi_nn - Pi_tt) dn = kappa * int sum_i c_i u\'_{i,n}^2 dn          (exact, 2nd order)')
print('zero flux J_i = 0:  c_i u\'_i = D d/dn c_i + (D/cs^2) kappa d/dn (c_i u\'_i^2)')
print('leading order:  u\'_i = D d/dn ln c_i  ->')
print('sigma = kappa D^2 sum_i int (d c_i/dn)^2 / c_i dn')
xm = sp.Function('x_1')(x)
cc = sp.Symbol('c', positive=True)
integrand = sp.simplify(((cc * xm).diff(x))**2 / (cc * xm) + ((cc * (1 - xm)).diff(x))**2 / (cc * (1 - xm)))
print('binary, c = const:  sum_i (c_i\')^2/c_i =', sp.factor(integrand))
th = sp.Function('theta')(x)
sub = sp.simplify(integrand.subs(xm, sp.sin(th)**2).doit())
print('with x_1 = sin^2(theta):  =', sp.simplify(sub))
print('=> sigma = 4 kappa c D^2 int (d theta/dn)^2 dn,  theta = arcsin(sqrt(x_1))')
print('   for a profile of width l (theta linear in n):  sigma ~ 4 kappa c D^2 (dtheta)^2 / l')

print('\n' + '=' * 78)
print('5. Diffusion velocity of the code, Gibbs-Duhem')
print('=' * 78)
# code: u'_i = -D grad ln gamma_i + B (B = D grad c, same for every species)
#       u'_i += -sum_j x_j u'_j   (uxcorr)
g1, g2, Bs, x1 = sp.symbols('g_1 g_2 B x_1', real=True)
A1, A2 = -D * g1 + Bs, -D * g2 + Bs
corr = -(x1 * A1 + (1 - x1) * A2)
print("u'_1 =", sp.factor(sp.expand(A1 + corr)), "   (the grad c term B cancels)")
a12, a21 = sp.symbols('A_12 A_21', real=True)
xs = sp.Symbol('x', positive=True)
lng1 = (a12 + 2 * (a21 - a12) * xs) * (1 - xs)**2
lng2 = (a21 + 2 * (a12 - a21) * (1 - xs)) * xs**2
print('Margules Gibbs-Duhem x1 dlng1 + x2 dlng2 =',
      sp.simplify(xs * lng1.diff(xs) + (1 - xs) * lng2.diff(xs)))
print("=> u'_i = -D grad ln gamma_i exactly; at equilibrium grad ln(x_i gamma_i) = 0 (to O(kappa u'^2/cs^2)),")
print("   so u'_i = D grad ln x_i: grad ln gamma_i is the interface indicator, grad mu_i = grad ln a_i is ~0")

print('\n' + '=' * 78)
print('6. Galilean invariance: third-order moments, central Hermite equilibrium')
print('=' * 78)
# Target: each species a Gaussian with mean v_i = u + u'_i and covariance
# cs^2 I + (kappa - 1) u'_i u'_i  (second moments = scaled moments above).
# D2Q9 holds m21, m12, m22 but aliases m30 -> m10, m03 -> m01.
def gauss_raw(ci, vx_, vy_, Sxx, Syy, Sxy):
    return {(3, 0): ci * (vx_**3 + 3 * vx_ * Sxx), (0, 3): ci * (vy_**3 + 3 * vy_ * Syy),
            (2, 1): ci * (vx_**2 * vy_ + Sxx * vy_ + 2 * vx_ * Sxy),
            (1, 2): ci * (vy_**2 * vx_ + Syy * vx_ + 2 * vy_ * Sxy)}

def chm_raw(ci, wx, wy, kap2, kap3):
    """raw moments of the central-Hermite-moment equilibrium about u:
       a0 = c, a1 = c u', a2 = c kap2 u'u', a_xxy = c kap3 u'x^2 u'y, a_xyy = c kap3 u'x u'y^2"""
    a10, a01 = ci * wx, ci * wy
    a20, a02, a11 = ci * kap2 * wx**2, ci * kap2 * wy**2, ci * kap2 * wx * wy
    a21, a12 = ci * kap3 * wx**2 * wy, ci * kap3 * wx * wy**2
    k = {(0, 0): ci, (1, 0): a10, (0, 1): a01, (2, 0): a20 + cs2 * ci, (0, 2): a02 + cs2 * ci,
         (1, 1): a11, (2, 1): a21 + cs2 * a01, (1, 2): a12 + cs2 * a10}
    kk = lambda a, b: k.get((a, b), 0)
    raw = lambda a, b: sp.expand(sum(sp.binomial(a, i) * sp.binomial(b, j) * ux**(a - i) * uy**(b - j) * kk(i, j)
                                     for i in range(a + 1) for j in range(b + 1)))
    return {(2, 1): raw(2, 1), (1, 2): raw(1, 2), (3, 0): raw(1, 0), (0, 3): raw(0, 1)}   # D2Q9 aliasing

k3 = sp.Symbol('kappa_3', real=True)
w1x, w1y = p1x, p1y
w2x, w2y = -c1 * p1x / c2, -c1 * p1y / c2
def mix(fun, *extra):
    s1 = fun(c1, w1x, w1y, *extra); s2 = fun(c2, w2x, w2y, *extra)
    return {key: sp.expand(s1[key] + s2[key]) for key in s1}
tgt = lambda ci, wx, wy: gauss_raw(ci, ux + wx, uy + wy, cs2 + (kap - 1) * wx**2, cs2 + (kap - 1) * wy**2, (kap - 1) * wx * wy)
G = mix(tgt)
Hm = mix(chm_raw, kap, k3)
Pm = {}
for ci, wx, wy in ((c1, w1x, w1y), (c2, w2x, w2y)):     # Parker's equilibrium, raw rows 6, 7 and aliased 1, 2
    e = eqR.subs({c: ci, px: wx, py: wy})
    for key, row in (((2, 1), 6), ((1, 2), 7), ((3, 0), 1), ((0, 3), 2)):
        Pm[key] = sp.expand(Pm.get(key, 0) + e[row])
Txx = kap * (c1 * w1x**2 + c2 * w2x**2)
for label, LB in (('Parker (2nd-order Lallemand-Luo) equilibrium', Pm), ('central Hermite equilibrium', Hm)):
    print(f'\n{label}: target - LB third moments of the mixture, u-dependent part')
    for key in [(3, 0), (2, 1), (1, 2), (0, 3)]:
        d = sp.expand(G[key] - LB[key])
        d_u = sp.expand(d - d.subs({ux: 0, uy: 0}))
        print(f'  dQ_{"x" * key[0]}{"y" * key[1]} = {sp.factor(d_u) if d_u != 0 else 0}')
dxxx = sp.expand(G[(3, 0)] - Hm[(3, 0)])
dxxx_u = sp.expand(dxxx - dxxx.subs({ux: 0, uy: 0}))
print('\nwith the central Hermite equilibrium only the aliased xxx (and yyy) moment is wrong:')
print('  dQ_xxx(u) - (c u_x^3 + 3 u_x T_xx) =', sp.simplify(dxxx_u - ((c1 + c2) * ux**3 + 3 * ux * Txx)),
      '   T = kappa sum_i c_i u\'_i u\'_i (capillary stress)')
print('  -> viscous stress error -(1/omega_nu - 1/2) d/dx[c u_x^3 + 3 u_x T_xx]: the standard cubic defect')
print('     plus an interface term (capillary stress advected by u). Correction: a source in the xx (yy)')
print('     stress moment, S_xx = -d/dx dQ_xxx, entered as (1 - omega_nu/2) S_xx in the collision.')
