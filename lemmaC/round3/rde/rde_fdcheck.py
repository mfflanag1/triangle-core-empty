"""Finite-difference check of DT = -2P (displacement coordinates) and of the shift law."""
import numpy as np
from scipy.interpolate import CubicSpline
from rde_solve import solve, T
h = 0.02; o, x, f, G, G1 = solve(L=10.0, h=h, tol=1e-13); n = len(x)
cs = CubicSpline(x, f)
def pushed(eps, phi, dphi):
    # density of X + eps*phi(X): solve u + eps*phi(u) = x by fixed-point iteration
    u = x.copy()
    for _ in range(60): u = x - eps * phi(u)
    return np.clip(cs(np.clip(u, x[0], x[-1])), 0, None) / (1 + eps * dphi(u))
t = 2 * x[0] + h * np.arange(2 * n - 1)
g1 = np.array([(np.clip(tt + x, 0, None) ** 2 * f).sum() * h / 2 for tt in t])
K = g1[np.add.outer(np.arange(n), np.arange(n))]
m = (x > -3) & (x < 4)
for name, phi, dphi in [('sin', np.sin, np.cos), ('const', lambda u: 1 + 0 * u, lambda u: 0 * u),
                        ('tanh', np.tanh, lambda u: 1 / np.cosh(u) ** 2)]:
    eps = 1e-3
    Gp = T(pushed(eps, phi, dphi), x, h)[0]; Gm = T(pushed(-eps, phi, dphi), x, h)[0]
    psi = -(Gp - Gm) / (2 * eps * G1)
    Pphi = (h * f[None, :] * K * phi(x)[None, :]).sum(1) / G1
    print(name, 'max|psi_fd + 2 P phi| on [-3,4] = %.2e' % np.max(np.abs(psi[m] + 2 * Pphi[m])),
          ' max|2 P phi| = %.3f' % np.max(np.abs(2 * Pphi[m])))
# convergence rate of the re-centred iteration (mixture damping 1/3) from a perturbed start
f2 = np.exp(-(x - 1.0) ** 2 / 0.6); f2 /= f2.sum() * h; d_old = None
for it in range(40):
    Gn, G1n, _, _ = T(f2, x, h); fn = G1n * np.exp(-Gn); d = np.abs(fn - f2).sum() * h
    if it % 5 == 4: print(' it', it + 1, 'L1 defect %.3e' % d, 'ratio/iter %.3f' % ((d / d_old) ** 0.2 if d_old else float('nan'))); d_old = d
    f2 = (2 * f2 + fn) / 3
