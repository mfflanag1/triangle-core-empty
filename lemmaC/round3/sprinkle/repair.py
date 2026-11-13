"""Size of the integer and fractional repairs after deleting one pair of M. Floating point; evidence only.
For sampled pairs i of M: kIP = number of M-teams not in M_2 (optimum of I^{-i});
kLP = number of M-teams with weight < 1 in the LP optimum of I^{-i}; nfrac = number of fractional LP teams.
Seeds: SeedSequence([20261111, 3, 3, q, idx])."""
import os, sys, json
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling'))
from cg import solve_lp
from regret import solve

def run(args):
    q, idx, ns = args
    rng = np.random.default_rng(np.random.SeedSequence([20261111, 3, 3, q, idx]))
    mats = [rng.random((q, q)) for _ in range(3)]; A, B, C = mats
    LP, T, sel = solve(A, B, C, q); D = LP - T
    a, b, c = np.unravel_index(sel, (q, q, q)); Mset = set(int(s) for s in sel)
    pairs = [(0, a[k], b[k]) for k in range(q)] + [(1, a[k], c[k]) for k in range(q)] + [(2, b[k], c[k]) for k in range(q)]
    out = []
    for j in rng.choice(3 * q, size=ns, replace=False):
        t, u, w = pairs[j]; old = mats[t][u, w]; mats[t][u, w] = 0.0
        LP2, T2, sel2 = solve(A, B, C, q)
        V = A[:, :, None] + B[:, None, :] + C[None, :, :]
        _, _, _, cols, x, _ = solve_lp(V, q)
        mats[t][u, w] = old
        w_on_M = {int(cc): float(xx) for cc, xx in zip(cols, x) if int(cc) in Mset}
        kLP = sum(1 for s in Mset if w_on_M.get(s, 0.0) < 1 - 1e-7)
        out.append(dict(D2=LP2 - T2, rIP=T - T2, kIP=int(q - len(Mset & set(int(s) for s in sel2))), kLP=int(kLP),
                        nfrac=int(np.sum((x > 1e-7) & (x < 1 - 1e-7)))))
    return dict(q=q, idx=idx, D=D, out=out)

if __name__ == '__main__':
    q, n, ns = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    with ProcessPoolExecutor(2) as ex, open(os.path.join(HERE, f'repair_q{q:03d}.jsonl'), 'a') as f:
        for rec in ex.map(run, [(q, i, ns) for i in range(n)]):
            f.write(json.dumps(rec, default=float) + '\n'); f.flush()
