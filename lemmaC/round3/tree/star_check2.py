"""Re-test the one perturbation that gave L > L* in star_check.py, with ln w renormalised
(max = 0) so that the Laplace grid in t covers the support of W."""
import numpy as np
from star_first_moment import xi, lam, stats
from star_check import L_of
beta = 0.55729
lnw = -3 * beta * xi
for it in range(300):
    rho, _ = stats(lnw); rho = np.clip(rho, 1e-300, 1 - 1e-12)
    new = -beta * xi + (2 / 3) * (np.log(rho) - np.log1p(-rho)); new -= new.max()
    lnw = 0.5 * new + 0.5 * lnw
Ls, _, _ = L_of(lnw, beta)
for eps in (0.3, 0.6, 1.0, 2.0):
    p = lnw + eps * xi**2 / 4; raw = L_of(p, beta); p = p - p.max(); nrm = L_of(p, beta)
    rho, _ = stats(p)
    print("xi2 eps=%.1f raw L-L*=%+.3e | normalised L-L*=%+.3e Sigma=%.3f e=%.3f norm=%.4f" % (eps, raw[0] - Ls, nrm[0] - Ls, nrm[1], nrm[2], (rho * lam).sum()))
