"""v15: (i) identities at the numerical fixed point: E S_+^3 = 6, E S_+^4/8 = 9 E X = c*.
(ii) Lemma 3 adversarially: maximise Gamma(Q1) - Gamma(Q0) over dispersive-ordered pairs q0 <= q1.
(iii) cost comparison: minimise c(Q) = 3 Gamma(Q) over IFR laws Q below the fixed point (q <= q*, q nonincreasing)."""
import numpy as np
from scipy.optimize import minimize
from iter_cn import step
from adv_M import Xnodes

d = np.load('fixedpoint_h8.npz'); x = d['x']; f = d['f']; G = d['G']; Gp = d['Gp']; cstar = float(d['c'])
h = x[1] - x[0]; J = (len(x) - 1) // 2
# (i)
G2, Gp2 = step(f, h, J); x2 = np.arange(-2 * J, 2 * J + 1) * h
fn = Gp2 * np.exp(-G2); ET = h * np.sum(x2 * fn); EXc = h * np.sum(x * f); k = ET - EXc; a = k / 3
fY = h * np.convolve(f, f); fS = h * np.convolve(fY, f); x3 = np.arange(-3 * J, 3 * J + 1) * h
Sp = np.maximum(x3 + 3 * a, 0)
print("(i) fixed point: shift k = %.8f, E X* = %.8f, 9 E X* = %.8f, E S_+^3 = %.7f (should be 6), E S_+^4/8 = %.7f, c* (iteration) = %.8f"
      % (k, EXc + a, 9 * (EXc + a), h * np.sum(fS * Sp ** 3), h * np.sum(fS * Sp ** 4) / 8, cstar))
# check T(X_c) = X_c + k in law: compare survival functions
i0 = np.searchsorted(x2, x[0] + 1e-12)
Gshift = np.interp(x + k, x2, G2)
ok = np.isfinite(G) & (G > 1e-12) & (G < 30)
print("    sup |G_T(x+k) - G(x)|/G(x) over 1e-12<G<30: %.2e" % np.max(np.abs(Gshift[ok] - G[ok]) / G[ok]))

# theta-mesh machinery
def nodes(K, m, th_lo=1e-4, th_hi=16.0):
    bp = np.concatenate([[0.0], np.geomspace(th_lo, th_hi, K)])
    gl, gw = np.polynomial.legendre.leggauss(m)
    th = []; w = []; cell = []
    for kk in range(len(bp) - 1):
        lo_, hi_ = bp[kk], bp[kk + 1]; t = (lo_ + hi_) / 2 + (hi_ - lo_) / 2 * gl
        th += list(t); w += list((hi_ - lo_) / 2 * gw * np.exp(-t)); cell += [kk] * m
    ll, lw = np.polynomial.laguerre.laggauss(4)
    th += list(bp[-1] + ll); w += list(lw * np.exp(-bp[-1])); cell += [len(bp) - 1] * 4
    return bp, np.array(th), np.array(w), np.array(cell)

glt, gwt = np.polynomial.legendre.leggauss(6)
def Gamma(q, env):
    bp, th, w, cell = env
    X = Xnodes(q, bp, th, cell)
    Y = (X[:, None] + X[None, :]).ravel(); W = (w[:, None] * w[None, :]).ravel()
    tmin = -Y.max(); EY = np.sum(W * Y)
    edges = np.concatenate([tmin + np.geomspace(1e-6, 1.0, 20) * 0 + np.linspace(0, 1, 21)[:-1] * 0, np.linspace(tmin, -EY + 12.0, 90)])
    edges = np.unique(edges)
    tl, tr = edges[:-1], edges[1:]
    t = ((tl + tr) / 2)[:, None] + ((tr - tl) / 2)[:, None] * glt[None, :]
    wt = ((tr - tl) / 2)[:, None] * gwt[None, :]
    t = t.ravel(); wt = wt.ravel()
    U = np.maximum(t[:, None] + Y[None, :], 0)
    g = (U ** 3 @ W) / 6
    ET_ = tmin + np.sum(wt * np.exp(-g))
    return ET_ + 2 * np.sum(w * X)

rng = np.random.default_rng(3)
env = nodes(24, 5); bp = env[0]; n = len(bp)
# q* at the right end of each theta-cell (last cell: 0 -> atom, the certificate's rule)
okk = np.isfinite(G) & (G > 0)
qstar = lambda th_: np.exp(np.interp(np.log(th_), np.log(G[okk]), -np.log(Gp[okk])))
cap = np.concatenate([qstar(bp[1:]), [0.0]])
c_cap = 3 * Gamma(cap, env)
print("(iii) c* = %.6f; one-sided rounding of Q* on %d cells (q = q* at right ends, atom on top): c = %.6f" % (cstar, n, c_cap))
def q_from(p):
    q = cap / (1 + np.exp(-p)); return np.minimum.accumulate(q)
vals = [3 * Gamma(q_from(rng.normal(1, 3, n)), env) for _ in range(200)]
print("     200 random IFR laws below Q*: min c = %.6f, max c = %.6f  (claim: all >= c* = %.6f)" % (min(vals), max(vals), cstar))
best = 9.0
for s in range(6):
    res = minimize(lambda p: 3 * Gamma(q_from(p), env), rng.normal(2, 2, n), method='Powell', options=dict(maxfev=3000, xtol=1e-3, ftol=1e-10))
    best = min(best, res.fun)
    print("     adversarial start %d: min c over IFR laws below Q* = %.6f  (gap to cap law %.2e, to c* %+.2e)" % (s, res.fun, res.fun - c_cap, res.fun - cstar), flush=True)
# scalings tau*Q* (tau<1) and T^n(0) are also below Q*: tau*cap
for tau in (0.0, 0.25, 0.5, 0.75, 0.9, 1.0):
    print("     tau = %.2f: c(tau * rounded Q*) = %.6f" % (tau, 3 * Gamma(tau * cap, env)))
# (ii) Lemma 3: any pair q0 <= q1 (not nec. IFR): Gamma(Q0) >= Gamma(Q1)
env2 = nodes(12, 6); n2 = len(env2[0])
def pair(p, ifr):
    e = np.exp(np.clip(p[:n2], -30, 6)); q1 = np.cumsum(e[::-1])[::-1] if ifr else e
    q0 = q1 / (1 + np.exp(-p[n2:]))
    if ifr: q0 = np.minimum.accumulate(q0)
    return q0, q1
def gap(p, ifr=True):
    q0, q1 = pair(p, ifr); return -(Gamma(q1, env2) - Gamma(q0, env2))
for ifr in (True, False):
    w_ = max(-gap(np.concatenate([rng.normal(-1, 2, n2), rng.normal(0, 3, n2)]), ifr) for _ in range(200))
    best2 = -9
    for s in range(4):
        res = minimize(gap, np.concatenate([rng.normal(-1, 2, n2), rng.normal(0, 3, n2)]), args=(ifr,), method='Powell', options=dict(maxfev=2500, xtol=1e-3, ftol=1e-10))
        best2 = max(best2, -res.fun)
    print("(ii) Lemma 3, %s pairs q0<=q1: max over 200 random of Gamma(Q1)-Gamma(Q0) = %+.3e; adversarial max = %+.3e  (claim: <= 0)" % ("IFR" if ifr else "general monotone", w_, best2), flush=True)
