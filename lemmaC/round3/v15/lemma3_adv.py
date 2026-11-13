"""v15: Lemma 3 (Gamma antitone in the dispersive order) adversarially, general nondecreasing Q (no concavity),
slopes bounded, maximiser re-evaluated on finer quadrature.  Functions copied from cost_cmp.py."""
import numpy as np, time
from scipy.optimize import minimize
from adv_M import Xnodes
def nodes(K, m, th_lo=1e-4, th_hi=16.0):
    bp = np.concatenate([[0.0], np.geomspace(th_lo, th_hi, K)])
    gl, gw = np.polynomial.legendre.leggauss(m); th = []; w = []; cell = []
    for kk in range(len(bp) - 1):
        lo_, hi_ = bp[kk], bp[kk + 1]; t = (lo_ + hi_) / 2 + (hi_ - lo_) / 2 * gl
        th += list(t); w += list((hi_ - lo_) / 2 * gw * np.exp(-t)); cell += [kk] * m
    ll, lw = np.polynomial.laguerre.laggauss(8)
    th += list(bp[-1] + ll); w += list(lw * np.exp(-bp[-1])); cell += [len(bp) - 1] * 8
    return bp, np.array(th), np.array(w), np.array(cell)
glt, gwt = np.polynomial.legendre.leggauss(8)
def Gamma(q, env, npan=120):
    bp, th, w, cell = env; X = Xnodes(q, bp, th, cell)
    Y = (X[:, None] + X[None, :]).ravel(); W = (w[:, None] * w[None, :]).ravel()
    tmin = -Y.max(); EY = np.sum(W * Y); edges = np.linspace(tmin, -EY + 12.0, npan + 1)
    tl, tr = edges[:-1], edges[1:]
    t = (((tl + tr) / 2)[:, None] + ((tr - tl) / 2)[:, None] * glt[None, :]).ravel(); wt = (((tr - tl) / 2)[:, None] * gwt[None, :]).ravel()
    g = (np.maximum(t[:, None] + Y[None, :], 0) ** 3 @ W) / 6
    return tmin + np.sum(wt * np.exp(-g)) + 2 * np.sum(w * X)
env = nodes(12, 6); envf = nodes(12, 20); n = len(env[0])
def pair(p):
    q1 = np.exp(np.clip(p[:n], -8, 1.5)); return q1 / (1 + np.exp(-p[n:])), q1
def gap(p, e=env, npan=120):
    q0, q1 = pair(p); return -(Gamma(q1, e, npan) - Gamma(q0, e, npan))
rng = np.random.default_rng(21); t0 = time.time()
for s in range(5):
    res = minimize(gap, np.concatenate([rng.normal(-1, 1.5, n), rng.normal(0, 3, n)]), method='Nelder-Mead', options=dict(maxfev=700, xatol=1e-4, fatol=1e-13))
    q0, q1 = pair(res.x)
    print("start %d: max Gamma(Q1)-Gamma(Q0): coarse %+.3e, fine quadrature %+.3e  (Gamma(Q1) = %.4f, max q1/q0 - 1 = %.2e)  %.0fs"
          % (s, -res.fun, -gap(res.x, envf, 400), Gamma(q1, envf, 400), np.max(q1 / q0 - 1), time.time() - t0), flush=True)
