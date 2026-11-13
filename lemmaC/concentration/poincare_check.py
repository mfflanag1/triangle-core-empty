"""Numerical sanity check of the exponential-Poincare bounds in CONCENTRATION.md (evidence only).
Per instance: T, LP*, mu* (signed, by bisection), the maximizing vertex w* just below mu*, its
d(w*), the Poincare integrand Q = sum_i c_i^2 (1-X_i)^2 with c = (pairload(w*) - pairload(M))/d,
and R_T = sum_{M-pairs}(1-X_i)^2. Compare Var(mu*) with 4 E Q and Var(T) with 4 E R_T."""
import json, sys, os, argparse
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'gate1'))
from gate1 import instance, structure, ok

def lpsolve(c, Mm, n):
    r = ok(linprog(-c, A_eq=Mm.T.tocsr(), b_eq=np.ones(n), bounds=(0, None), method='highs-ds'), 'LP')
    return r.x, -r.fun

def pairloads(x, q):
    X = x.reshape(q, q, q)  # v[a,b,c] = A[a,b] + B[a,c] + C[b,c]
    return np.concatenate([X.sum(2).ravel(), X.sum(1).ravel(), X.sum(0).ravel()])

def run(args):
    q, idx = args
    A, B, C = instance('triangle', 'unif', q, idx)
    v, Mm = structure(A, B, C); n = 3 * q
    Xp = np.concatenate([A.ravel(), B.ravel(), C.ravel()])
    lam, V = lpsolve(v, Mm, n)
    mi = milp(-v, constraints=LinearConstraint(Mm.T.tocsr(), np.ones(n), np.ones(n)), integrality=np.ones(q**3),
              bounds=Bounds(0, 1), options=dict(time_limit=900, mip_rel_gap=1e-12))
    ind = (mi.x > 0.5).astype(float); T = float(v @ ind)
    lo, hi = -6.0, 1.0
    def isM(mu):
        x, _ = lpsolve(v + mu * ind, Mm, n); return np.max(np.abs(x - ind)) < 1e-6, x
    while not isM(hi)[0]: hi *= 2
    xl = isM(lo)[1]
    for _ in range(45):
        mid = (lo + hi) / 2
        good, x = isM(mid)
        if good: hi = mid
        else: lo, xl = mid, x
    d = float((1 - xl[ind > 0.5]).sum())
    c = (pairloads(xl, q) - pairloads(ind, q)) / d
    Q = float((c**2 * (1 - Xp)**2).sum())
    mpairs = pairloads(ind, q) > 0.5
    R_T = float(((1 - Xp[mpairs])**2).sum())
    teamdef = 3 - v[ind > 0.5]
    barD = float(((1 - xl[ind > 0.5]) * teamdef).sum() / d)
    return dict(q=q, idx=idx, T=T, LP=float(V), mu=lo, d=d, Q=Q, c1=float(np.abs(c).sum()), cinf=float(np.abs(c).max()),
                R_T=R_T, def_T=3*q - T, dmax=float(teamdef.max()), barD=barD, Qbound=(2 * barD - lo) / d)

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--n', type=int, nargs='+'); ap.add_argument('--tag', default='')
    a = ap.parse_args()
    for q, nn in zip(a.q, a.n):
        with ProcessPoolExecutor(4) as ex, open(os.path.join(HERE, 'out', f'poin{a.tag}_q{q:02d}.jsonl'), 'w') as f:
            for r in ex.map(run, [(q, i) for i in range(nn)]):
                f.write(json.dumps(r) + '\n'); f.flush()
        print('done', q, flush=True)
