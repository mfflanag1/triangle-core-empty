"""Signed stability margin of the integer optimum M:  mu* = inf{mu : M is an optimal LP vertex for
v + mu * 1[team in M]}.  Core empty iff mu* > 0 (a.s.).  For nonempty cores mu* <= 0 measures how
far M is from losing LP optimality.  Bisection on mu in [-3, 3]. Floating point; evidence only."""
import json, sys, os, argparse
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'gate1'))
from gate1 import instance, structure, ok

def run(args):
    q, idx = args
    A, B, C = instance('triangle', 'unif', q, idx)
    v, Mm = structure(A, B, C); n = 3 * q; Ae = Mm.T.tocsr()
    sol = lambda c: ok(linprog(-c, A_eq=Ae, b_eq=np.ones(n), bounds=(0, None), method='highs-ds'), 'LP').x
    x = sol(v)
    if np.max(np.minimum(x, 1 - x)) < 1e-6:
        ind = (x > 0.5).astype(float)
    else:
        mi = milp(-v, constraints=LinearConstraint(Ae, np.ones(n), np.ones(n)), integrality=np.ones(q**3),
                  bounds=Bounds(0, 1), options=dict(time_limit=900, mip_rel_gap=1e-12))
        ind = (mi.x > 0.5).astype(float)
    isM = lambda mu: np.max(np.abs(sol(v + mu * ind) - ind)) < 1e-6
    lo, hi = -3.0, 3.0   # isM(hi) True; isM(lo) False (mu=-3 makes M worst)
    for _ in range(30):
        mid = (lo + hi) / 2
        if isM(mid): hi = mid
        else: lo = mid
    return dict(q=q, idx=idx, mu_star=(lo + hi) / 2, empty=bool((lo + hi) / 2 > 0))

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--n', type=int)
    a = ap.parse_args()
    for q in a.q:
        with ProcessPoolExecutor(6) as ex, open(os.path.join(HERE, 'out', f'margin_q{q:02d}.jsonl'), 'w') as f:
            for r in ex.map(run, [(q, i) for i in range(a.n)]):
                f.write(json.dumps(r) + '\n'); f.flush()
        print('done', q, flush=True)
