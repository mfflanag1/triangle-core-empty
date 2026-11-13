"""How many kinds of pair-keeping teams does the gap need?  Floating point; evidence only.  Same seeds as clone.py.
Restricted LPs over M + teams sharing no pair with M + pair-keeping teams of the listed kinds
(kind AB = teams (a_i, b_i, c_k) keeping the AB pair of M-team i, etc.)."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import sys, json
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from clone import amat
HERE = os.path.dirname(os.path.abspath(__file__))
def run(q, idx):
    rng = np.random.default_rng(np.random.SeedSequence([20261111, 3, 1, q, idx]))
    A, B, C = rng.random((q, q)), rng.random((q, q)), rng.random((q, q))
    vf = (A[:, :, None] + B[None, :, :] + C[:, None, :]).reshape(-1); Ae = amat(q); n = 3 * q
    lp = linprog(-vf, A_eq=Ae, b_eq=np.ones(n), bounds=(0, None), method='highs-ds'); LPs = -lp.fun
    if np.max(np.minimum(lp.x, 1 - lp.x)) < 1e-6: x = np.round(lp.x); T = LPs
    else:
        mi = milp(-vf, constraints=LinearConstraint(Ae, np.ones(n), np.ones(n)), integrality=np.ones(q ** 3), bounds=Bounds(0, 1),
                  options=dict(mip_rel_gap=0.0)); x = np.round(mi.x); T = -mi.fun
    Mi = np.nonzero(x > 0.5)[0]; a, b, c = np.unravel_index(np.arange(q ** 3), (q, q, q))
    ta = np.empty(q, int); tb = np.empty(q, int); tc = np.empty(q, int)
    ma, mb, mc = np.unravel_index(Mi, (q, q, q)); ta[ma] = np.arange(q); tb[mb] = np.arange(q); tc[mc] = np.arange(q)
    i, j, k = ta[a], tb[b], tc[c]
    base = ((i == j) & (j == k)) | ((i != j) & (j != k) & (i != k))
    kAB = (i == j) & (k != i); kBC = (j == k) & (i != j); kAC = (i == k) & (j != i)
    out = dict(q=q, idx=idx, D=LPs - T)
    for nm, mask in (('AB', base | kAB), ('AB+BC', base | kAB | kBC), ('AB+AC', base | kAB | kAC), ('BC+AC', base | kBC | kAC)):
        cols = np.nonzero(mask)[0]
        r = linprog(-vf[cols], A_eq=Ae[:, cols], b_eq=np.ones(n), bounds=(0, None), method='highs-ds'); out[nm] = -r.fun - T
    return out
if __name__ == '__main__':
    for q, nn in ((8, 300), (12, 200), (16, 200)):
        R = [run(q, i) for i in range(nn)]
        json.dump(R, open(os.path.join(HERE, f'kinds_q{q:02d}.json'), 'w'))
        D = np.array([r['D'] for r in R])
        print(f"q={q} n={nn} P(D>0)={(D > 1e-7).mean():.3f} ED={D.mean():.4f} | " + ' '.join(
            f"{nm}: P={(np.array([r[nm] for r in R]) > 1e-7).mean():.3f} E={np.mean([r[nm] for r in R]):.4f}" for nm in ('AB', 'AB+BC', 'AB+AC', 'BC+AC')), flush=True)
