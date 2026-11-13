"""Check of the deletion flow balance (REPORT Prop. S3), floating point, evidence only:
   E[N_kill ; D>0] = E[N_create ; D=0],   E sum_i (D^{-i} - D) = E D,
N_create = #{i in M : D(I^{i->0}) > 0} (on D=0),  N_kill = #{pairs i : D(I^{i->0}) = 0} (on D>0; only pairs of M or
pairs loaded by the LP vertex can matter).  All such pairs are deleted one at a time (value set to 0).
Seeds: SeedSequence([20261111, 3, 2, q, idx])."""
import os, sys, json
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling'))
from cg import solve_lp
from regret import solve
EPS = 1e-7

def run(args):
    q, idx = args
    rng = np.random.default_rng(np.random.SeedSequence([20261111, 3, 2, q, idx]))
    mats = [rng.random((q, q)) for _ in range(3)]; A, B, C = mats
    LP, T, sel = solve(A, B, C, q); D = LP - T
    V = A[:, :, None] + B[:, None, :] + C[None, :, :]
    _, _, _, cols, x, _ = solve_lp(V, q)
    def prs(teams):
        a, b, c = np.unravel_index(teams, (q, q, q))
        return {(0, int(u), int(w)) for u, w in zip(a, b)} | {(1, int(u), int(w)) for u, w in zip(a, c)} | {(2, int(u), int(w)) for u, w in zip(b, c)}
    PM = prs(sel); PX = prs(cols[x > EPS])
    res = []
    for (t, u, w) in sorted(PM | PX):
        old = mats[t][u, w]; mats[t][u, w] = 0.0
        LP2, T2, _ = solve(A, B, C, q); mats[t][u, w] = old
        res.append((int((t, u, w) in PM), T - T2, LP - LP2))
    return dict(q=q, idx=idx, T=T, LP=LP, D=D, nM=len(PM), nX=len(PX), nU=len(PM | PX), res=res)

if __name__ == '__main__':
    q, n = int(sys.argv[1]), int(sys.argv[2]); start = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    with ProcessPoolExecutor(2) as ex, open(os.path.join(HERE, f'balance_q{q:03d}.jsonl'), 'a') as f:
        for rec in ex.map(run, [(q, i) for i in range(start, start + n)]):
            f.write(json.dumps(rec, default=float) + '\n'); f.flush()
