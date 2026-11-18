"""
Exact linear theory of the shape oscillation (mode n) of a 2D droplet of
radius R in an unbounded fluid, both incompressible with the same density
rho and kinematic viscosity nu, surface tension sigma.

Normal modes ~ exp(i n theta + s t), u = grad phi + curl(psi e_z):
  inside   phi = A r^n,     psi = B I_n(k r)
  outside  phi = C r^-n,    psi = E K_n(k r),    k^2 = s/nu,  p = -rho s phi
Interface r = R + eta: continuity of u_r, u_theta and the shear stress,
normal-stress jump tau_rr(in) - tau_rr(out) = -sigma (n^2-1) eta / R^2,
kinematic condition s eta = u_r.  s = -gamma + i omega.

usage: python3 droplet_dispersion.py sigma R nu [n]
"""
import sys
import numpy as np
from scipy import special, optimize


def _potential(m, n, s, R, mu, rho):
    """phi = r^m: u_r, u_theta, tau_r_theta, tau_rr at R"""
    ur = m * R**(m - 1)
    ut = 1j * n * R**(m - 1)
    trt = 2 * mu * 1j * n * (m - 1) * R**(m - 2)
    trr = rho * s * R**m + 2 * mu * m * (m - 1) * R**(m - 2)
    return ur, ut, trt, trr


def _vortical(Z, Zp, n, k, R, mu):
    """psi = Z(k r) (modified Bessel I_n or K_n): u_r, u_theta, tau_r_theta, tau_rr at R"""
    z = k * R
    Z0, Z1 = Z(n, z), Zp(n, z)
    Z2 = -Z1 / z + (1 + n * n / z**2) * Z0
    ur = 1j * n / R * Z0
    ut = -k * Z1
    trt = mu * (-k * k * Z2 + k * Z1 / R - n * n * Z0 / R**2)
    trr = 2 * mu * 1j * n * (k * Z1 / R - Z0 / R**2)
    return ur, ut, trt, trr


def det(s, sigma, R, nu, n=2, rho=1.0):
    mu = rho * nu
    k = np.sqrt(s / nu)                 # principal root: Re k > 0, K_n decays outside
    cols = [_potential(n, n, s, R, mu, rho), _vortical(special.iv, special.ivp, n, k, R, mu),
            _potential(-n, n, s, R, mu, rho), _vortical(special.kv, special.kvp, n, k, R, mu)]
    M = np.empty((4, 4), complex)
    for j, (ur, ut, trt, trr) in enumerate(cols):
        sign = 1 if j < 2 else -1
        eta = ur / s if j < 2 else 0
        M[:, j] = sign * np.array([ur, ut, trt, trr + sigma * (n * n - 1) / R**2 * eta])
        M[:, j] /= np.abs(M[:, j]).max()
    return np.linalg.det(M)


def mode(sigma, R, nu, n=2, rho=1.0):
    w0 = np.sqrt(n * (n * n - 1) * sigma / (2 * rho * R**3))
    g0 = 2 * n * (n - 1) * nu / R**2
    f = lambda v: [det(complex(v[0], v[1]), sigma, R, nu, n, rho).real, det(complex(v[0], v[1]), sigma, R, nu, n, rho).imag]
    sol = optimize.root(f, [-g0, w0], tol=1e-14)
    return complex(*sol.x), w0, g0


def sigma_from(omega, R, nu, n=2, rho=1.0):
    """surface tension for which the exact mode has frequency omega; returns (sigma, gamma)"""
    g = lambda sg: mode(sg, R, nu, n, rho)[0].imag - omega
    sg0 = omega**2 * 2 * rho * R**3 / (n * (n * n - 1))
    sg = optimize.brentq(g, 0.5 * sg0, 3 * sg0, xtol=1e-12)
    return sg, -mode(sg, R, nu, n, rho)[0].real


if __name__ == '__main__':
    sigma, R, nu = map(float, sys.argv[1:4])
    n = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    s, w0, g0 = mode(sigma, R, nu, n)
    print(f'exact: omega={s.imag:.5e} gamma={-s.real:.5e}   inviscid omega0={w0:.5e}, '
          f'irrotational gamma={g0:.5e}   omega/omega0={s.imag / w0:.4f}  gamma/gamma0={-s.real / g0:.3f}')
