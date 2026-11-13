"""v15: exact re-evaluation (root solve, no interpolation in t) of the 'violations' the Powell search in adv_M.py reported."""
import numpy as np
from scipy.optimize import minimize, brentq
import adv_M as A

def nodes(bp, m, mlag, split_first=1):
    gl, gw = np.polynomial.legendre.leggauss(m)
    th = []; w = []; cell = []
    for k in range(len(bp) - 1):
        a, b = bp[k], bp[k + 1]
        sub = np.linspace(a, b, (split_first if k == 0 else 1) + 1)
        for u, v in zip(sub[:-1], sub[1:]):
            t = (u + v) / 2 + (v - u) / 2 * gl
            th += list(t); w += list((v - u) / 2 * gw * np.exp(-t)); cell += [k] * m
    ll, lw = np.polynomial.laguerre.laggauss(mlag)
    th += list(bp[-1] + ll); w += list(lw * np.exp(-bp[-1])); cell += [len(bp) - 1] * mlag
    return np.array(th), np.array(w), np.array(cell)

def s_exact(q, bp, th, w, cell, levels):
    X = A.Xnodes(q, bp, th, cell)
    Y = (X[:, None] + X[None, :]).ravel(); W = (w[:, None] * w[None, :]).ravel()
    g = lambda t: np.sum(W * np.maximum(t + Y, 0) ** 3) / 6
    lo = -Y.max(); out = []
    for L in levels:
        hi = -Y.min() + (6 * L / W.sum()) ** (1 / 3) + 1.0
        t = brentq(lambda t: g(t) - L, lo, hi, xtol=1e-13, rtol=1e-14)
        out.append(np.sum(W * np.maximum(t + Y, 0) ** 2) / 2)
    return np.array(out)

rng = np.random.default_rng(1)
bp, th, w, cell = A.make_nodes()
levels = np.geomspace(1e-3, 20, 80); env = (bp, th, w, cell, levels); n = len(bp)
for _ in range(40): rng.normal(0, 2, n), rng.normal(0, 2, n)
for _ in range(300):
    rng.normal(-1, 2, n), rng.normal(0, 3, n); rng.random()
lev2 = np.geomspace(1e-3, 20, 40)
for start in range(8):
    p0 = np.concatenate([rng.normal(-1, 2, n), rng.normal(0, 3, n)])
    if start % 3 == 0: p0[n - 1] = -20
    if start % 5 == 4: p0[n:] = rng.normal(4, 1, n)
    if start not in (1, 3, 5, 6): continue
    res = minimize(A.violation, p0, args=(env,), method='Powell', options=dict(maxfev=1200, xtol=1e-3, ftol=1e-9))
    q0, q1 = A.pair_from_params(res.x, n)
    print("start %d: interpolated objective %+.3e" % (start, -res.fun))
    print("   q0 =", np.array2string(q0, precision=3)); print("   q1 =", np.array2string(q1, precision=3))
    for (m, ml, sf) in ((5, 4, 1), (20, 12, 1), (20, 24, 8)):
        th2, w2, c2 = nodes(bp, m, ml, sf)
        s0 = s_exact(q0, bp, th2, w2, c2, lev2); s1 = s_exact(q1, bp, th2, w2, c2, lev2)
        r = (s1 - s0) / s0
        print("   exact solve, %2d GL nodes/cell, %2d Laguerre, first cell split %d: max (s1-s0)/s0 = %+.3e at theta=%.3g" % (m, ml, sf, r.max(), lev2[r.argmax()]), flush=True)
