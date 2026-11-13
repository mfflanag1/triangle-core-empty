"""Linearisation of T at the fixed point in displacement coordinates: DT = -2 P,
P(x,da) = f(a) g1(x+a) da / G'(x),  g1(t) = E (t+X)_+^2 / 2.  P is a Markov kernel, reversible
w.r.t. pi ~ f G'.  Spectrum via the symmetrised matrix.  Also a finite-difference check of DT."""
import numpy as np
from rde_solve import solve, T
for h, L in [(0.04, 10.0), (0.02, 10.0)]:
    o, x, f, G, G1 = solve(L=L, h=h, tol=1e-13)
    n = len(x)
    # g1 on grid t = 2*x[0] + h*k
    t = 2 * x[0] + h * np.arange(2 * n - 1)
    g1 = np.array([(np.clip(tt + x, 0, None) ** 2 * f).sum() * h / 2 for tt in t])
    K = g1[np.add.outer(np.arange(n), np.arange(n))]          # g1(x_i + x_j)
    keep = (f > 1e-13) & (G1 > 1e-13)
    w = np.sqrt(h * f[keep] / G1[keep])
    A = w[:, None] * K[np.ix_(keep, keep)] * w[None, :]
    ev = np.linalg.eigvalsh(A)
    rows = (h * f[None, keep] * K[np.ix_(keep, keep)] / G1[keep, None]).sum(1)
    print('h', h, 'n', keep.sum(), 'row sums min/max', rows.min(), rows.max())
    print(' top eigenvalues of P:', np.round(ev[::-1][:8], 6))
    print(' bottom eigenvalues of P:', np.round(ev[:8], 6))
    print(' => eigenvalues of DT=-2P: shift', -2 * ev[-1], ' next largest |.|:', 2 * max(abs(ev[0]), abs(ev[-2])))
    print(' min |1+2p| =', np.min(np.abs(1 + 2 * ev)), '  => ||(I-DT)^-1||_{L2(pi)} ~', 1 / np.min(np.abs(1 + 2 * ev)))
# finite-difference check of DT = -2P on a test displacement phi(x) = sin(x)
o, x, f, G, G1 = solve(L=10.0, h=0.02, tol=1e-13); h = 0.02; n = len(x)
eps = 1e-4; phi = np.sin(x)
xi = x + eps * phi                                           # displaced atoms, weights h*f
# density of displaced law on grid: push forward f by x -> x+eps*phi
fe = np.interp(x, xi, f / (1 + eps * np.cos(x)))
Ge, G1e, _, _ = T(fe, x, h)
psi_fd = -(Ge - G) / (eps * np.where(G1 > 1e-12, G1, np.nan))    # output displacement
t = 2 * x[0] + h * np.arange(2 * n - 1)
g1 = np.array([(np.clip(tt + x, 0, None) ** 2 * f).sum() * h / 2 for tt in t])
K = g1[np.add.outer(np.arange(n), np.arange(n))]
Pphi = (h * f[None, :] * K * phi[None, :]).sum(1) / np.where(G1 > 1e-12, G1, np.nan)
m = (x > -3) & (x < 4)
print('FD check: max |psi_fd - (-2 P phi)| on [-3,4] =', np.nanmax(np.abs(psi_fd[m] + 2 * Pphi[m])), ' scale', np.nanmax(np.abs(2 * Pphi[m])))
