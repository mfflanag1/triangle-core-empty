"""Sister model: iid team values. Floating point; evidence only.
Models: P  = pair-additive (the brief's model), v = A[x,y]+B[x,z]+C[y,z];
        IH = each team value an independent sum of three U[0,1] (same marginal as P, so P(deficit<t)=t^3/6 for t<=1);
        U  = each team value an independent U[0,1].
Per instance: LP* (column generation from ../../scaling/cg.py), integrality of the LP optimum, T and D by reduced-cost-fixing MILP."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import sys, json, time, argparse
import numpy as np
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling'))
import cg
from scout import integer_opt_sel
TAG = 770032

def run(args):
    model, q, idx, noD = args
    t0 = time.time()
    rng = np.random.default_rng(np.random.SeedSequence([20261111, TAG, {'P': 0, 'IH': 1, 'U': 2}[model], q, idx]))
    if model == 'P':
        A, B, C = rng.random((q, q)), rng.random((q, q)), rng.random((q, q))
        V = A[:, :, None] + B[:, None, :] + C[None, :, :]
    elif model == 'IH':
        V = rng.random((3, q, q, q)).sum(0)
    else:
        V = rng.random((q, q, q))
    vf = V.reshape(-1)
    LPs, p, r, cols, x, it = cg.solve_lp(V, q)
    integral = bool(np.max(np.minimum(x, 1 - x)) < 1e-6)
    if integral: T = LPs
    elif noD: T = None
    else:
        _, T = integer_opt_sel(vf, q, LPs, r, 0.15)
    return dict(model=model, q=q, idx=idx, LP=LPs, T=T, D=(None if T is None else LPs - T), integral=integral, secs=time.time() - t0)

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--models', nargs='+', default=['IH', 'P', 'U']); ap.add_argument('--q', type=int, nargs='+')
    ap.add_argument('--n', type=int); ap.add_argument('--start', type=int, default=0); ap.add_argument('--workers', type=int, default=1); ap.add_argument('--noD', action='store_true')
    a = ap.parse_args()
    os.makedirs(os.path.join(HERE, 'out'), exist_ok=True)
    for q in a.q:
        for m in a.models:
            fn = os.path.join(HERE, 'out', f'sister_{m}_q{q:03d}.jsonl')
            with ProcessPoolExecutor(a.workers) as ex, open(fn, 'a') as f:
                for rr in ex.map(run, [(m, q, i, a.noD) for i in range(a.start, a.start + a.n)]):
                    f.write(json.dumps(rr) + '\n'); f.flush()
            print('done', m, q, flush=True)
