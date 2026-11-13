"""Single-deletion regrets on the random triangle market. Floating point; evidence only.
For an instance I with integer optimum M, LP optimum x: for sampled pairs i in M set the pair value to 0
(deficit 1) and recompute T, LP*.  rho_IP = T - T^{-i}, rho_LP = LP - LP^{-i}, D^{-i} = LP^{-i} - T^{-i}.
Identity under test (REPORT Lemma S2): E sum_{i in M} rho_IP(i) = E(3q - T) up to exp-small terms.
Seeds: SeedSequence([20261111, 3, 1, q, idx])."""
import os, sys, json
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling'))
from cg import solve_lp, amat

def solve(A, B, C, q):
    V = A[:, :, None] + B[:, None, :] + C[None, :, :]
    vf = V.reshape(-1)
    LPs, p, r, cols, x, it = solve_lp(V, q)
    if np.max(np.minimum(x, 1 - x)) < 1e-6:
        sel = cols[x > 0.5]; return LPs, float(vf[sel].sum()), sel
    delta = 0.15; n = 3 * q
    while True:
        cand = np.nonzero(r >= -delta)[0]
        mi = milp(-vf[cand], constraints=LinearConstraint(amat(cand, q), np.ones(n), np.ones(n)),
                  integrality=np.ones(len(cand)), bounds=Bounds(0, 1), options=dict(time_limit=600, mip_rel_gap=1e-12))
        if mi.status == 0 and -mi.fun >= LPs - delta + 1e-9:
            return LPs, -mi.fun, cand[mi.x > 0.5]
        delta *= 2

def run(args):
    q, idx, ns = args
    rng = np.random.default_rng(np.random.SeedSequence([20261111, 3, 1, q, idx]))
    A, B, C = rng.random((q, q)), rng.random((q, q)), rng.random((q, q))
    LP, T, sel = solve(A, B, C, q)
    a, b, c = np.unravel_index(sel, (q, q, q))
    pairs = [(0, a[k], b[k]) for k in range(q)] + [(1, a[k], c[k]) for k in range(q)] + [(2, b[k], c[k]) for k in range(q)]  # cg.py: V[x,y,z]=A[x,y]+B[x,z]+C[y,z]
    mats = [A, B, C]
    out = []
    for j in rng.choice(3 * q, size=min(ns, 3 * q), replace=False):
        t, u, w = pairs[j]
        old = mats[t][u, w]; mats[t][u, w] = 0.0
        LP2, T2, _ = solve(A, B, C, q)
        mats[t][u, w] = old
        out.append(dict(d=1 - old, rIP=T - T2, rLP=LP - LP2, D2=LP2 - T2))
    return dict(q=q, idx=idx, LP=LP, T=T, D=LP - T, del_=out)

if __name__ == '__main__':
    q, n, ns = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    with ProcessPoolExecutor(2) as ex, open(os.path.join(HERE, f'regret_q{q:03d}.jsonl'), 'a') as f:
        for rec in ex.map(run, [(q, i, ns) for i in range(n)]):
            f.write(json.dumps(rec, default=float) + '\n'); f.flush()
