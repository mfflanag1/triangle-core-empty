"""Rotation-restricted LP. Floating point; evidence only.
D_rot = LP*(teams sharing at least one pair with M) - T.  Always 0 <= D_rot <= D (M is feasible and optimal).
Seeds: SeedSequence([20261111, 3, 4, q, idx])."""
import os, sys, json
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.optimize import linprog
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling'))
from cg import amat
from regret import solve

def run(args):
    q, idx = args
    rng = np.random.default_rng(np.random.SeedSequence([20261111, 3, 4, q, idx]))
    A, B, C = rng.random((q, q)), rng.random((q, q)), rng.random((q, q))
    LP, T, sel = solve(A, B, C, q)
    V = A[:, :, None] + B[:, None, :] + C[None, :, :]
    a, b, c = np.unravel_index(sel, (q, q, q))
    keep = np.zeros((q, q, q), bool)
    keep[a, b, :] = True; keep[a, :, c] = True; keep[:, b, c] = True
    cols = np.nonzero(keep.reshape(-1))[0]
    lp = linprog(-V.reshape(-1)[cols], A_eq=amat(cols, q), b_eq=np.ones(3 * q), bounds=(0, None), method='highs-ds')
    return dict(q=q, idx=idx, D=LP - T, Drot=-lp.fun - T, ncols=int(len(cols)))

if __name__ == '__main__':
    q, n = int(sys.argv[1]), int(sys.argv[2])
    with ProcessPoolExecutor(2) as ex:
        R = list(ex.map(run, [(q, i) for i in range(n)]))
    with open(os.path.join(HERE, f'drot_q{q:03d}.jsonl'), 'w') as f:
        for r in R: f.write(json.dumps(r, default=float) + '\n')
    D = np.array([r['D'] for r in R]); Dr = np.array([r['Drot'] for r in R])
    print(f'q={q} n={n} E D={D.mean():.4f} E D_rot={Dr.mean():.4f} ratio={Dr.mean()/max(D.mean(),1e-12):.3f}  P(D>0)={(D>1e-7).mean():.3f} P(D_rot>0)={(Dr>1e-7).mean():.3f}  min(D-D_rot)={np.min(D-Dr):.2e} min D_rot={Dr.min():.2e}')
