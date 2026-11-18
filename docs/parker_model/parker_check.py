"""
Numerical check of the surface-tension relation of parker_surface_tension.py
with a small NumPy implementation of Parker's 2D model (MRT, periodic),
Margules activity coefficients (parameters of margules.py / periodic.f90).

  flat   : planar interfaces; sigma_mech = int (Pi_xx - Pi_yy) dx against
           kappa sum_i c_i u'_i^2 (model) and kappa D^2 sum_i (c_i')^2/c_i
  laplace: droplets of several radii; dp R against sigma_mech

usage: python3 parker_check.py flat|laplace|osc [D] [kappa] [drift] [R (laplace) | key=value (osc: nu R eps n nsteps coll)] [coll (flat)]
  coll: mrt (Parker) | chm (central Hermite moments) | chmg (chm + xxx/yyy correction)
  drift A: Parker, u'_i = -D G[ln gamma_i]
        B: u'_i = -D G[ln a_i] + D G[c_i]/c_i
        C: u'_i = -D G[ln a_i] + D/(cs2 c_i) G.Pi_i^eq (lagged; compensates the kappa stress)
  (all followed by the mass-average correction u'_i -= sum_j x_j u'_j)
"""
import sys
import numpy as np
from scipy import optimize

cx = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1])
cy = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1])
cs2 = 1 / 3
c2 = cx**2 + cy**2
Mp = np.array([np.ones(9), -4 + 3 * c2, 4 - 10.5 * c2 + 4.5 * c2**2, cx, (-5 + 3 * c2) * cx,
               cy, (-5 + 3 * c2) * cy, cx**2 - cy**2, cx * cy])
Mpi = np.linalg.inv(Mp)
RAW = [(0, 0), (1, 0), (0, 1), (2, 0), (0, 2), (1, 1), (2, 1), (1, 2), (2, 2)]
Mri = np.linalg.inv(np.array([cx**a * cy**b for (a, b) in RAW], float))

R_, T_ = 8.314, 293.0
Am, Bm = 6359.0 / (R_ * T_), -384.0 / (R_ * T_)
a12, a21 = Am - Bm, Am + Bm


def lngamma(x1):
    x2 = 1 - x1
    return np.array([(a12 + 2 * (a21 - a12) * x1) * x2**2, (a21 + 2 * (a12 - a21) * x2) * x1**2])


def coexistence():
    def fun(v):
        xa, xb = v
        ga, gb = lngamma(xa), lngamma(xb)
        return [np.log(xa) + ga[0] - np.log(xb) - gb[0], np.log(1 - xa) + ga[1] - np.log(1 - xb) - gb[1]]
    return optimize.root(fun, [0.05, 0.95]).x


def grad(phi):
    """isotropic gradient (Parker's s = 4): weights 1/3 straight, 1/12 diagonal"""
    gx = np.zeros_like(phi); gy = np.zeros_like(phi)
    for k in range(1, 9):
        w = 1 / 3 if k <= 4 else 1 / 12
        sh = np.roll(phi, (-cy[k], -cx[k]), axis=(-2, -1))
        gx += w * cx[k] * sh
        gy += w * cy[k] * sh
    return gx, gy


class Parker:
    def __init__(self, x1, D, nu, kappa, c0=1.0, wE=1.0, drift='A', coll='mrt'):
        """coll: 'mrt' (Parker), 'chm' (central Hermite moments about u),
                 'chmg' (chm + source for the aliased xxx/yyy moments)"""
        self.D, self.kappa, self.drift, self.coll = D, kappa, drift, coll
        wD, wV = 1 / (D / cs2 + 0.5), 1 / (nu / cs2 + 0.5)
        self.wD, self.wV, self.wE = wD, wV, wE
        self.K = np.array([0, wE, wE, wD, wD, wD, wD, wV, wV])[:, None, None]
        ci = np.array([c0 * x1, c0 * (1 - x1)])
        w = np.array([4 / 9] + [1 / 9] * 4 + [1 / 36] * 4)
        self.f = ci[:, None] * w[None, :, None, None]
        self.Pi = np.array([ci * cs2, ci * cs2, 0 * ci])     # lagged species stress (xx, yy, xy) for drift C

    def macro(self):
        f = self.f
        ci = f.sum(axis=1)
        c = ci.sum(axis=0)
        jx = np.tensordot(cx, f, axes=(0, 1)).sum(axis=0)
        jy = np.tensordot(cy, f, axes=(0, 1)).sum(axis=0)
        x1 = ci[0] / c
        lg = lngamma(x1)
        xi = np.array([x1, 1 - x1])
        px = np.empty_like(ci); py = np.empty_like(ci)
        for i in range(2):
            if self.drift == 'A':          # Parker: -D grad ln gamma_i
                gx, gy = grad(lg[i])
            else:                          # chemical-potential form: -D grad ln a_i + LB Fickian part
                gx, gy = grad(np.log(xi[i]) + lg[i])
                if self.drift == 'B':      # + D grad c_i / c_i
                    fx, fy = grad(ci[i])
                else:                      # + D/(cs2 c_i) div Pi_i^eq (lagged, includes kappa stress)
                    fx = grad(self.Pi[0, i])[0] + grad(self.Pi[2, i])[1]
                    fy = grad(self.Pi[2, i])[0] + grad(self.Pi[1, i])[1]
                    fx, fy = fx / cs2, fy / cs2
                gx, gy = gx - fx / ci[i], gy - fy / ci[i]
            px[i], py[i] = -self.D * gx, -self.D * gy
        corr_x = -(ci * px).sum(axis=0) / c
        corr_y = -(ci * py).sum(axis=0) / c
        return ci, c, jx / c, jy / c, px + corr_x, py + corr_y, x1, lg

    def step(self):
        if self.coll != 'mrt':
            return self.step_chm()
        ci, c, ux, uy, px, py, x1, lg = self.macro()
        kap = self.kappa
        m = np.einsum('ak,ikyx->iayx', Mp, self.f)
        U2 = ux**2 + uy**2 + 2 * (ux * px + uy * py) + kap * (px**2 + py**2)
        eq = np.array([ci, ci * (-2 + 3 * U2), ci * (1 - 3 * U2),
                       ci * (ux + px), -ci * (ux + px), ci * (uy + py), -ci * (uy + py),
                       ci * (ux**2 - uy**2 + 2 * px * ux - 2 * py * uy + kap * (px**2 - py**2)),
                       ci * (ux * uy + px * uy + py * ux + kap * px * py)]).transpose(1, 0, 2, 3)
        self.mneq = m - eq
        self.Pi = np.array([ci * (cs2 + 2 * ux * px + kap * px**2), ci * (cs2 + 2 * uy * py + kap * py**2),
                            ci * (ux * py + uy * px + kap * px * py)])
        m = m - self.K[None] * (m - eq)
        fp = np.einsum('ka,iayx->ikyx', Mpi, m)
        self.fhalf = 0.5 * (self.f + fp).sum(axis=0)
        for k in range(1, 9):
            fp[:, k] = np.roll(fp[:, k], (cy[k], cx[k]), axis=(-2, -1))
        self.f = fp

    def step_chm(self):
        """collision in central Hermite moments about the mixture velocity u.
        Equilibrium a0 = c_i, a1 = c_i u'_i, a2 = kappa c_i u'_i u'_i, a3 = a4 = 0 (Parker's
        equilibrium at rest); rates: a1 and a3 with omega_D (as Parker's j and q), deviatoric a2
        with omega_nu, trace of a2 and a4 with omega_E. 'chmg' adds the source
        S = -d/dx dQ_xxx (and yy), dQ_xxx = c u_x^3 + 3 u_x T_xx, to the xx (yy) stress."""
        ci, c, ux, uy, px, py, x1, lg = self.macro()
        kap = self.kappa
        f = self.f
        dx = cx[:, None, None] - ux[None]           # (9, ny, nx)
        dy = cy[:, None, None] - uy[None]
        k = {(a, b): np.einsum('kyx,ikyx->iyx', dx**a * dy**b, f) for (a, b) in
             [(0, 0), (1, 0), (0, 1), (2, 0), (0, 2), (1, 1), (2, 1), (1, 2), (2, 2)]}
        # central Hermite moments
        h = {(0, 0): k[0, 0], (1, 0): k[1, 0], (0, 1): k[0, 1],
             (2, 0): k[2, 0] - cs2 * k[0, 0], (0, 2): k[0, 2] - cs2 * k[0, 0], (1, 1): k[1, 1],
             (2, 1): k[2, 1] - cs2 * k[0, 1], (1, 2): k[1, 2] - cs2 * k[1, 0],
             (2, 2): k[2, 2] - cs2 * (k[2, 0] + k[0, 2]) + cs2**2 * k[0, 0]}
        heq = {(0, 0): ci, (1, 0): ci * px, (0, 1): ci * py,
               (2, 0): kap * ci * px**2, (0, 2): kap * ci * py**2, (1, 1): kap * ci * px * py,
               (2, 1): 0 * ci, (1, 2): 0 * ci, (2, 2): 0 * ci}
        wD, wV, wE = self.wD, self.wV, self.wE
        hp = {}
        hp[0, 0] = h[0, 0]
        for key in [(1, 0), (0, 1), (2, 1), (1, 2)]:
            hp[key] = h[key] - wD * (h[key] - heq[key])
        hp[1, 1] = h[1, 1] - wV * (h[1, 1] - heq[1, 1])
        hp[2, 2] = h[2, 2] - wE * (h[2, 2] - heq[2, 2])
        tr, trq = h[2, 0] + h[0, 2], heq[2, 0] + heq[0, 2]
        dv, dvq = h[2, 0] - h[0, 2], heq[2, 0] - heq[0, 2]
        tr = tr - wE * (tr - trq)
        dv = dv - wV * (dv - dvq)
        if self.coll == 'chmg':
            T = kap * (ci * px**2).sum(axis=0), kap * (ci * py**2).sum(axis=0)
            dQx = c * ux**3 + 3 * ux * T[0]
            dQy = c * uy**3 + 3 * uy * T[1]
            Sxx, Syy = -grad(dQx)[0], -grad(dQy)[1]
            xi = ci / c                               # shared by the species in proportion to x_i
            tr = tr + (1 - 0.5 * wE) * xi * (Sxx + Syy)
            dv = dv + (1 - 0.5 * wV) * xi * (Sxx - Syy)
        hp[2, 0], hp[0, 2] = 0.5 * (tr + dv), 0.5 * (tr - dv)
        # back to central moments, raw moments (binomial shift), populations
        kp = {(0, 0): hp[0, 0], (1, 0): hp[1, 0], (0, 1): hp[0, 1],
              (2, 0): hp[2, 0] + cs2 * hp[0, 0], (0, 2): hp[0, 2] + cs2 * hp[0, 0], (1, 1): hp[1, 1],
              (2, 1): hp[2, 1] + cs2 * hp[0, 1], (1, 2): hp[1, 2] + cs2 * hp[1, 0]}
        kp[2, 2] = hp[2, 2] + cs2 * (kp[2, 0] + kp[0, 2]) - cs2**2 * hp[0, 0]
        from math import comb
        mr = []
        for (a, b) in RAW:
            mr.append(sum(comb(a, i) * comb(b, j) * ux**(a - i) * uy**(b - j) * kp[i, j]
                          for i in range(a + 1) for j in range(b + 1)))
        fp = np.einsum('ka,aiyx->ikyx', Mri, np.array(mr))
        self.fhalf = 0.5 * (f + fp).sum(axis=0)
        # bookkeeping for drift C (lagged species stress)
        self.Pi = np.array([ci * (cs2 + 2 * ux * px + kap * px**2), ci * (cs2 + 2 * uy * py + kap * py**2),
                            ci * (ux * py + uy * px + kap * px * py)])
        for q in range(1, 9):
            fp[:, q] = np.roll(fp[:, q], (cy[q], cx[q]), axis=(-2, -1))
        self.f = fp

    def stress(self):
        """mixture momentum flux to second order: average of pre- and post-collision second moments"""
        fh = self.fhalf
        Pxx = np.tensordot(cx**2, fh, axes=(0, 0))
        Pyy = np.tensordot(cy**2, fh, axes=(0, 0))
        return Pxx, Pyy


def run_flat(D, kappa, drift='A', nx=80, ny=2, nsteps=40000, coll='mrt'):
    xa, xb = coexistence()
    X = np.arange(nx) + 0.5
    prof = 0.5 * (np.tanh((X - nx / 4) / 2) - np.tanh((X - 3 * nx / 4) / 2))
    x1 = (xa + (xb - xa) * prof)[None, :].repeat(ny, axis=0)
    s = Parker(x1, D, nu=0.1, kappa=kappa, drift=drift, coll=coll)
    for t in range(nsteps):
        s.step()
    s.step()
    ci, c, ux, uy, px, py, x1, lg = s.macro()
    Pxx, Pyy = s.stress()
    sig_mech = 0.5 * (Pxx - Pyy)[0].sum()                       # two interfaces
    sig_model = 0.5 * kappa * (ci * px**2).sum(axis=0)[0].sum()
    dci = np.array([grad(ci[i])[0] for i in range(2)])
    sig_formula = 0.5 * kappa * D**2 * (dci**2 / ci).sum(axis=0)[0].sum()
    th = np.arcsin(np.sqrt(np.clip(x1[0], 0, 1)))
    width = np.abs(th.max() - th.min()) / np.abs(np.diff(th)).max()
    lna = np.log(x1) + lg[0]
    print(f'[{drift},{coll}] D={D:5.3f} kappa={kappa:5.2f}  x_coex=({xa:.4f},{xb:.4f}) sim=({x1.min():.4f},{x1.max():.4f})'
          f'  width(theta)={width:4.2f} cells  max|u|={np.abs(ux).max():.1e}')
    print(f'   sigma_mech={sig_mech:.5e}  kappa*sum c u\'^2={sig_model:.5e}  kappa D^2 sum c\'^2/c={sig_formula:.5e}'
          f'   I=sigma/(kappa D^2 c)={sig_mech / (kappa * D**2 * c.mean()):.4f}')
    print(f'   max|grad ln gamma_1|={np.abs(grad(lg[0])[0]).max():.3f}   max|grad ln a_1|={np.abs(grad(lna)[0]).max():.2e}'
          f'   Pxx range={np.ptp(Pxx):.1e}')
    return sig_mech


def run_laplace(D, kappa, drift='A', R=16, n=80, nsteps=15000):
    """static droplet: sigma_Laplace = dp R (2D)"""
    xa, xb = coexistence()
    Y, X = np.mgrid[0:n, 0:n] + 0.5
    r = np.hypot(X - n / 2, Y - n / 2)
    x1 = xa + (xb - xa) * 0.5 * (1 - np.tanh((r - R) / 2))
    s = Parker(x1, D, nu=0.1, kappa=kappa, drift=drift)
    for t in range(nsteps):
        s.step()
    s.step()
    ci, c, ux, uy, px, py, x1, lg = s.macro()
    Pxx, Pyy = s.stress()
    p = 0.5 * (Pxx + Pyy)
    lo, hi = x1.min(), x1.max()
    Req = np.sqrt(((x1 - lo) / (hi - lo)).sum() / np.pi)
    pin, pout = p[r < Req - 6].mean(), p[r > Req + 6].mean()
    print(f'[{drift}] D={D} kappa={kappa} R0={R:3d} R={Req:6.2f}  dp={pin - pout:.4e}  sigma=dp*R={(pin - pout) * Req:.4e}'
          f'  max|u|={np.hypot(ux, uy).max():.1e}')


def run_oscillation(D, kappa, drift='A', nu=0.02, R=20, eps=0.08, n=112, nsteps=12000, every=20, coll='mrt'):
    """n = 2 shape oscillation of a droplet. Prints t, deformation (Ixx - Iyy)/(Ixx + Iyy), R_eq and the
    dynamic surface tension sigma_eff = kappa sum_i int c_i |u'_i|^2 dA / perimeter (= sigma_mech when static).
    Inviscid 2D theory (equal densities rho = c): omega^2 = n (n^2 - 1) sigma / ((rho_in + rho_out) R^3)"""
    xa, xb = coexistence()
    Y, X = np.mgrid[0:n, 0:n] + 0.5
    a, b = R * (1 + eps), R / (1 + eps)
    r = np.hypot((X - n / 2) / a, (Y - n / 2) / b)
    dist = (r - 1) * np.sqrt(a * b)
    x1 = xa + (xb - xa) * 0.5 * (1 - np.tanh(dist / 2))
    s = Parker(x1, D, nu=nu, kappa=kappa, drift=drift, coll=coll)
    out = []
    for t in range(nsteps + 1):
        if t % every == 0:
            ci, c, ux, uy, px, py, x1, lg = s.macro()
            w = np.clip((x1 - xa) / (xb - xa), 0, 1)
            Ixx = (w * (X - n / 2)**2).sum(); Iyy = (w * (Y - n / 2)**2).sum()
            gx, gy = grad(w)
            E = kappa * (ci * (px**2 + py**2)).sum()
            out.append((t, (Ixx - Iyy) / (Ixx + Iyy), np.sqrt(w.sum() / np.pi), E / np.hypot(gx, gy).sum(),
                        np.hypot(ux, uy).max()))
        s.step()
    np.savetxt(sys.stdout, np.array(out), fmt='%.6e', header='t deformation R_eq sigma_eff max|u|')


if __name__ == '__main__':
    mode = sys.argv[1]
    D = float(sys.argv[2]) if len(sys.argv) > 2 else 0.1
    kappa = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    drift = sys.argv[4] if len(sys.argv) > 4 else 'A'
    if mode == 'flat':
        run_flat(D, kappa, drift, coll=sys.argv[5] if len(sys.argv) > 5 else 'mrt')
    elif mode == 'laplace':
        run_laplace(D, kappa, drift, R=int(sys.argv[5]) if len(sys.argv) > 5 else 16)
    else:   # osc D kappa drift key=value ... (nu, R, eps, n, nsteps)
        kw = {k: (int(v) if k in ('R', 'n', 'nsteps') else v if k == 'coll' else float(v)) for k, v in (a.split('=') for a in sys.argv[5:])}
        run_oscillation(D, kappa, drift, **kw)
