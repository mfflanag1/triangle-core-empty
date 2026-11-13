"""v15: (a) c(tau Q*) on the x-lattice for tau near 1 (laws below Q* for tau<1): Gamma antitone => c(tau Q*) >= c* iff tau <= 1.
(b) infinitesimal adversarial test of Theorem M with exact root solves: q1 = running-min(q0 (1 + eps d)), d in [0,1]."""
import numpy as np, time
from scipy.optimize import minimize
from iter_cn import step
import adv_M as A
from adv_M_diag import nodes, s_exact
d = np.load('fixedpoint_h8.npz'); x = d['x']; f = d['f']; cstar = float(d['c']); h = x[1] - x[0]; J = (len(x) - 1) // 2
for tau in (0.5, 0.9, 0.99, 0.999, 1.0, 1.001, 1.01, 1.1):
    ht = tau * h; ft = f / tau
    G2, Gp2 = step(ft, ht, J); x2 = np.arange(-2 * J, 2 * J + 1) * ht
    ET = ht * np.sum(x2 * Gp2 * np.exp(-G2)); EQ = ht * np.sum(x * tau * ft)
    print("(a) tau = %.3f: c(tau Q*) = %.8f   c - c* = %+.3e" % (tau, 3 * (ET + 2 * EQ), 3 * (ET + 2 * EQ) - cstar), flush=True)
bp, th, w, cell = A.make_nodes(); n = len(bp); lev = np.geomspace(2e-3, 15, 24); eps = 1e-3
th2, w2, c2 = nodes(bp, 20, 12, 4)
def pair(p):
    e = np.exp(np.clip(p[:n], -12, 2.5))
    if p[n - 1] < -8: e[-1] = 0.0
    q0 = np.cumsum(e[::-1])[::-1]; dd = 1 / (1 + np.exp(-p[n:])); dd = dd / dd.max()
    return q0, np.minimum.accumulate(q0 * (1 + eps * dd))
def obj(p, fine=False):
    q0, q1 = pair(p); a = (th2, w2, c2) if fine else (th, w, cell)
    s0 = s_exact(q0, bp, *a, lev); s1 = s_exact(q1, bp, *a, lev)
    return -np.max((s1 - s0) / s0) / eps
rng = np.random.default_rng(11); t0 = time.time(); worst = -9
r = [-obj(np.concatenate([rng.normal(-1, 1.5, n), rng.normal(0, 3, n)])) for _ in range(150)]
print("(b) 150 random IFR laws and directions: max of sup_theta d log s / d eps = %+.3e (claim <= 0)" % max(r), flush=True)
for start in range(8):
    p0 = np.concatenate([rng.normal(-1, 1.5, n), rng.normal(0, 3, n)])
    if start % 2 == 0: p0[n - 1] = -20
    res = minimize(obj, p0, method='Nelder-Mead', options=dict(maxfev=500, xatol=1e-4, fatol=1e-12))
    vf = -obj(res.x, fine=True); worst = max(worst, vf)
    print("    start %d: maximised sup_theta d log s/d eps: coarse nodes %+.3e, fine nodes %+.3e (%.0fs)" % (start, -res.fun, vf, time.time() - t0), flush=True)
    if time.time() - t0 > 200: break
print("(b) ADVERSARIAL MAX (fine nodes): %+.3e" % worst)
