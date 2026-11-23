"""
Fit the n = 2 oscillation histories of parker_check.py osc and compare with the
exact viscous linear theory (droplet_dispersion.py):
  sigma_osc  : surface tension for which the exact mode has the fitted frequency
  gamma_th   : damping of that exact mode, against the fitted damping
  sigma_eff  : dynamic surface tension from the capillary stress (time average and
               relative oscillation amplitude), against the static sigma_flat

usage: python3 fit_oscillation.py sigma_flat file:nu [file:nu ...]
"""
import os
import sys
import numpy as np
from scipy.optimize import curve_fit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from droplet_dispersion import mode, sigma_from

sig_flat = float(sys.argv[1])
for arg in sys.argv[2:]:
    fn, nu = arg.rsplit(':', 1)
    nu = float(nu)
    d = np.loadtxt(fn)
    t, e, R, se = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
    m = t >= 300
    f = lambda t, A, g, w, ph, b: A * np.exp(-g * t) * np.cos(w * t + ph) + b
    p, _ = curve_fit(f, t[m], e[m], p0=[0.1, 3e-4, 3e-3, 0, 0], maxfev=20000)
    A, g, w, ph, b = p
    Rm = R[m].mean()
    sig_osc, g_th = sigma_from(w, Rm, nu)
    s_flat = mode(sig_flat, Rm, nu)[0]
    late = t > t.max() / 2
    print(f'{os.path.basename(fn):12s} nu={nu:<6g} R={Rm:5.2f}  omega={w:.4e} (theory {s_flat.imag:.4e}, {100 * (w / s_flat.imag - 1):+5.1f}%)'
          f'  gamma={g:.3e} (theory {-s_flat.real:.3e}, x{g / -s_flat.real:.2f})'
          f'  sigma_osc={sig_osc:.4e} ({100 * (sig_osc / sig_flat - 1):+5.1f}%)'
          f'  sigma_eff={se[late].mean():.4e} ({100 * (se[late].mean() / sig_flat - 1):+5.1f}%, swing {100 * np.ptp(se[m]) / se[m].mean():.1f}%)')
