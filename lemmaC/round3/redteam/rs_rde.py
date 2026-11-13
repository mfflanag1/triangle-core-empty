"""Replica-symmetric (min-sum / tree) fixed point for d-index matching on the Poisson hypertree.
Teams at an agent: Poisson process of costs with intensity t^r/r!  (our model: d=3, r=2, cost unit q^(-2/3)).
RDE: X = min_k (xi_k - sum of d-1 iid copies of X).  P(X>x) = exp(-E[(x+S)_+^(r+1)]/(r+1)!), S = sum of d-1 copies.
Energy per sector-1 agent: a_RS = int t * t^r/r! * E[G(t-S)] dt.  Deterministic grid iteration; evidence only."""
import numpy as np, sys
from math import factorial
from scipy.signal import fftconvolve
def solve(d=3, r=2, L=14.0, h=0.01, iters=400, damp=0.5, verbose=False):
    x = np.arange(-L, L + h / 2, h); n = len(x)
    f = np.exp(-(x - 1.0) ** 2 / 0.5); f /= f.sum() * h
    hist = []
    for it in range(iters):
        g = f.copy()
        for _ in range(d - 2):
            g = fftconvolve(g, f) * h          # density of S on grid from -(d-1)L
        s = np.arange(len(g)) * h - (d - 1) * L * (1 if d > 2 else 1)
        if d == 2: s = x
        # Phi(x) = E[(x+S)_+^(r+1)]/(r+1)!
        # compute via direct sum on a coarser evaluation: Phi(x_i) = sum_j g_j h (x_i+s_j)_+^(r+1)
        # use correlation: kernel k(u)=u_+^(r+1) on grid
        u = np.arange(-(len(s) + n), (len(s) + n)) * h
        Phi = np.empty(n)
        # vectorised in chunks
        for a in range(0, n, 512):
            xx = x[a:a + 512, None] + s[None, :]
            Phi[a:a + 512] = (np.clip(xx, 0, None) ** (r + 1) * g[None, :]).sum(1) * h / factorial(r + 1)
        G = np.exp(-Phi)                       # P(X > x)
        fn = -np.gradient(G, h); fn = np.clip(fn, 0, None); fn /= fn.sum() * h
        diff = np.abs(fn - f).sum() * h
        f = damp * fn + (1 - damp) * f
        hist.append(diff)
        if diff < 1e-10: break
    # energy
    t = np.arange(0, 3 * L, h)
    # E[G(t-S)]: G extended: G=1 below -L, 0 above L
    def Gext(y):
        return np.where(y < -L, 1.0, np.where(y > L, 0.0, np.interp(y, x, G)))
    EG = np.empty(len(t))
    for a in range(0, len(t), 512):
        EG[a:a + 512] = (Gext(t[a:a + 512, None] - s[None, :]) * g[None, :]).sum(1) * h
    w = t ** r / factorial(r) * EG
    norm = w.sum() * h; en = (t * w).sum() * h
    mean = (x * f).sum() * h; sd = np.sqrt(((x - mean) ** 2 * f).sum() * h)
    return dict(d=d, r=r, iters=it + 1, last_diff=hist[-1], norm=norm, energy=en, meanX=mean, sdX=sd,
                q01=float(x[np.searchsorted(np.cumsum(f) * h, 0.01)]), q99=float(x[np.searchsorted(np.cumsum(f) * h, 0.99)]))
if __name__ == '__main__':
    print(solve(d=2, r=0, L=20, h=0.02))      # check: pi^2/6 = 1.6449
    print(solve(d=3, r=0, L=14, h=0.02))      # MMR uniform-cost 3-index model
    print(solve(d=3, r=2, L=10, h=0.02))      # our model
    print(solve(d=3, r=2, L=10, h=0.01))
