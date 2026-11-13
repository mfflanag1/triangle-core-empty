"""Spread sweep: pair values exp(K*(U-1)), U iid U[0,1]; K = 0 means plain U[0,1] values.
Large K approaches the ordinal model of Theorem O, where the core is always nonempty.  Reports the
empty share (fractional LP vertex) per (q, K).  Floating point; for K >= 20 values span ~1e-9 and
solver tolerance makes results unreliable, so K <= 16 only."""
import os, sys, argparse, json
import numpy as np
from scipy.optimize import linprog
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'gate1'))
from gate1 import instance, structure

def run(args):
    q, K, idx = args
    A, B, C = instance('triangle', 'unif', q, idx)
    f = (lambda X: X) if K == 0 else (lambda X: np.exp(K * (X - 1)))
    v, Mm = structure(f(A), f(B), f(C))
    v = v / v.max()
    x = linprog(-v, A_eq=Mm.T.tocsr(), b_eq=np.ones(3 * q), bounds=(0, None), method='highs-ds',
                options=dict(primal_feasibility_tolerance=1e-10, dual_feasibility_tolerance=1e-10)).x
    return bool(np.max(np.minimum(x, 1 - x)) > 1e-6)

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--n', type=int)
    a = ap.parse_args(); res = {}
    with ProcessPoolExecutor(4) as ex:
        for q in a.q:
            for K in [0, 1, 2, 4, 8, 16]:
                e = list(ex.map(run, [(q, K, i) for i in range(a.n)]))
                res[f'{q},{K}'] = float(np.mean(e)); print(q, K, 'P(empty)=%.3f' % np.mean(e), flush=True)
    json.dump(res, open(os.path.join(HERE, 'out', 'spread.json'), 'w'), indent=1)
