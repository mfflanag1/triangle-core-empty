"""v15: adversarial search for a violation of Theorem M with an EXACT objective (root solve in t, no interpolation)."""
import numpy as np, time
from scipy.optimize import minimize
import adv_M as A
from adv_M_diag import nodes, s_exact
bp, th, w, cell = A.make_nodes(); n = len(bp); lev = np.geomspace(2e-3, 15, 24)
th2, w2, c2 = nodes(bp, 20, 12, 4)
def pair(p):
    e = np.exp(np.clip(p[:n], -12, 2.5))
    if p[n - 1] < -8: e[-1] = 0.0
    q1 = np.cumsum(e[::-1])[::-1]; q0 = np.minimum.accumulate(q1 / (1 + np.exp(-p[n:]))); return q0, q1
def obj(p, fine=False):
    q0, q1 = pair(p)
    a = (th2, w2, c2) if fine else (th, w, cell)
    s0 = s_exact(q0, bp, *a, lev); s1 = s_exact(q1, bp, *a, lev)
    return -np.max((s1 - s0) / s0)
rng = np.random.default_rng(5); t0 = time.time(); worst = -9
for start in range(6):
    p0 = np.concatenate([rng.normal(-1, 1.5, n), rng.normal(0, 2, n)])
    if start % 2 == 0: p0[n - 1] = -20
    res = minimize(obj, p0, method='Nelder-Mead', options=dict(maxfev=450, xatol=1e-4, fatol=1e-12))
    vf = -obj(res.x, fine=True); worst = max(worst, vf)
    print("start %d: max (s1-s0)/s0 coarse nodes %+.3e, same pair on fine nodes %+.3e  (%.0fs)" % (start, -res.fun, vf, time.time() - t0), flush=True)
    if time.time() - t0 > 170: break
print("EXACT-OBJECTIVE ADVERSARIAL MAX (fine nodes): %+.3e" % worst)
