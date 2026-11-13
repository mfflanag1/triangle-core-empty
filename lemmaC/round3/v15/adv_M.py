"""v15: adversarial test of Theorem M (dispersive monotonicity of T on IFR laws) and of Lemma 3 / cost comparison.
Laws are given by slopes q_k >= 0 (nonincreasing = IFR) of the quantile function Q on a fixed theta-mesh
(X = Q(E), E ~ Exp(1)); q_last = 0 gives an atom at the top.  Q0 <= Q1 in dispersive order iff q0 <= q1.
g_Q(t) = E(t+X1+X2)_+^3/6 by tensor Gauss quadrature in theta; s_Q(theta) = g'(g^{-1}(theta)).
Violation of Theorem M  <=>  s_{Q1}(theta) > s_{Q0}(theta) for some theta."""
import sys, numpy as np
from scipy.optimize import minimize

def make_nodes(K=10, m=5, th_lo=1e-3, th_hi=14.0):
    bp = np.concatenate([[0.0], np.geomspace(th_lo, th_hi, K)])
    gl, gw = np.polynomial.legendre.leggauss(m)
    th = []; w = []; cell = []
    for k in range(len(bp) - 1):
        a, b = bp[k], bp[k + 1]
        t = (a + b) / 2 + (b - a) / 2 * gl
        th += list(t); w += list((b - a) / 2 * gw * np.exp(-t)); cell += [k] * m
    ll, lw = np.polynomial.laguerre.laggauss(4)
    th += list(bp[-1] + ll); w += list(lw * np.exp(-bp[-1])); cell += [len(bp) - 1] * 4
    return bp, np.array(th), np.array(w), np.array(cell)

def Xnodes(q, bp, th, cell):
    # Q piecewise linear, Q(bp[1])=0 reference (shift irrelevant)
    Qbp = np.concatenate([[0.0], np.cumsum(q[:-1] * np.diff(bp))])
    return Qbp[cell] + q[cell] * (th - bp[cell])

def curve(q, bp, th, w, cell, nt=300):
    X = Xnodes(q, bp, th, cell)
    Y = (X[:, None] + X[None, :]).ravel(); W = (w[:, None] * w[None, :]).ravel()
    o = np.argsort(Y); Y = Y[o]; W = W[o]
    sd = max(np.sqrt(np.sum(W * (Y - np.sum(W * Y)) ** 2)), 1e-3)
    t = np.linspace(-Y[-1] + 1e-9, -np.sum(W * Y) + 6 * sd + 6.0, nt)
    U = np.maximum(t[:, None] + Y[None, :], 0)
    return t, (U ** 3 @ W) / 6, (U ** 2 @ W) / 2, X

def s_at(levels, g, gp):
    ok = g > 1e-300
    return np.exp(np.interp(np.log(levels), np.log(g[ok]), np.log(gp[ok]), left=np.nan, right=np.nan))

def pair_from_params(p, n, ifr=True):
    # q1 = cumulative sums from the top (nonincreasing), q0 = running min of q1*sigmoid(c) (nonincreasing, <= q1)
    b = p[:n]; c = p[n:]
    e = np.exp(np.clip(b, -30, 8)); e[-1] = e[-1] if b[-1] > -8 else 0.0     # allow exact atom on top
    q1 = np.cumsum(e[::-1])[::-1] if ifr else e
    q0 = q1 / (1 + np.exp(-c))
    if ifr: q0 = np.minimum.accumulate(q0)
    return q0, q1

def violation(p, env, ifr=True, full=False):
    bp, th, w, cell, levels = env
    q0, q1 = pair_from_params(p, len(bp), ifr)
    t0, g0, gp0, _ = curve(q0, bp, th, w, cell); t1, g1, gp1, _ = curve(q1, bp, th, w, cell)
    s0 = s_at(levels, g0, gp0); s1 = s_at(levels, g1, gp1)
    r = (s1 - s0) / s0
    r = r[np.isfinite(r)]
    v = r.max() if len(r) else -1.0
    return (v, q0, q1) if full else -v

if __name__ == '__main__':
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    rng = np.random.default_rng(seed)
    bp, th, w, cell = make_nodes()
    levels = np.geomspace(1e-3, 20, 80)
    env = (bp, th, w, cell, levels); n = len(bp)
    print("nodes %d, slopes per law %d, quadrature mass %.12f" % (len(th), n, w.sum()))
    # control: non-IFR pairs must show violations (power of the test)
    best = -1
    for _ in range(40):
        p = np.concatenate([rng.normal(0, 2, n), rng.normal(0, 2, n)])
        best = max(best, violation(p, env, ifr=False, full=True)[0])
    print("control (non-IFR random pairs, 40): max (s1-s0)/s0 = %+.4f   (positive = violation, expected)" % best, flush=True)
    # random IFR pairs
    vals = []
    for _ in range(300):
        p = np.concatenate([rng.normal(-1, 2, n), rng.normal(0, 3, n)])
        if rng.random() < 0.3: p[n - 1] = -20     # atom on top of both
        vals.append(violation(p, env, full=True)[0])
    print("random IFR pairs (300): max (s1-s0)/s0 = %+.3e, median %+.3e" % (max(vals), np.median(vals)), flush=True)
    # adversarial: maximise the violation by Nelder-Mead / Powell from many starts
    worst = -1; worstp = None
    for start in range(8):
        p0 = np.concatenate([rng.normal(-1, 2, n), rng.normal(0, 3, n)])
        if start % 3 == 0: p0[n - 1] = -20
        if start % 5 == 4: p0[n:] = rng.normal(4, 1, n)       # q0 close to q1
        res = minimize(violation, p0, args=(env,), method='Powell', options=dict(maxfev=1200, xtol=1e-3, ftol=1e-9))
        v = -res.fun
        if v > worst: worst, worstp = v, res.x
        print("  start %2d: sup over levels of (s1-s0)/s0 after maximisation = %+.3e (nfev %d)" % (start, v, res.nfev), flush=True)
    print("ADVERSARIAL MAX over IFR pairs: %+.3e" % worst)
    # re-evaluate the worst pair on a finer quadrature
    q0, q1 = pair_from_params(worstp, n)
    for (K2, m2) in ((10, 5), (10, 30)):
        gl, gw = np.polynomial.legendre.leggauss(m2)
        th2 = []; w2 = []; cell2 = []
        for k in range(len(bp) - 1):
            a, b = bp[k], bp[k + 1]; tt = (a + b) / 2 + (b - a) / 2 * gl
            th2 += list(tt); w2 += list((b - a) / 2 * gw * np.exp(-tt)); cell2 += [k] * m2
        ll, lw = np.polynomial.laguerre.laggauss(10)
        th2 += list(bp[-1] + ll); w2 += list(lw * np.exp(-bp[-1])); cell2 += [len(bp) - 1] * 10
        th2 = np.array(th2); w2 = np.array(w2); cell2 = np.array(cell2)
        _, g0, gp0, _ = curve(q0, bp, th2, w2, cell2, nt=1500); _, g1, gp1, _ = curve(q1, bp, th2, w2, cell2, nt=1500)
        r = (s_at(levels, g1, gp1) - s_at(levels, g0, gp0)) / s_at(levels, g0, gp0)
        print("  worst pair re-evaluated with %d nodes/cell: max (s1-s0)/s0 = %+.3e at theta = %.3g" % (m2, np.nanmax(r), levels[np.nanargmax(r)]))
