"""Sister model (a3): iid team values with the SAME marginal law as the pair-additive model.
Each of the q^3 teams gets an independent deficit U1+U2+U3 (three fresh uniforms per team), value v = 3 - deficit.
Pair-sharing is removed; the one-team law, the first-moment bound (fact 4) and the local (tree) limit of the
near-tight hypergraph are unchanged. Floating point; evidence only. Compare with the brief's table.
Seeds: SeedSequence([20261111, 3, 2, q, idx])."""
import os, sys, json, time, argparse
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling')); sys.path.insert(0, HERE)
import cg
from iid_axial import integer_opt

def run(args):
    q, idx, need_ip = args
    t0 = time.time()
    rng = np.random.default_rng(np.random.SeedSequence([20261111, 3, 2, q, idx]))
    V = 3.0 - (rng.random((q, q, q)) + rng.random((q, q, q)) + rng.random((q, q, q)))
    vf = V.reshape(-1)
    LPs, p, r, cols, x, it = cg.solve_lp(V, q)
    integral = bool(np.max(np.minimum(x, 1 - x)) < 1e-6)
    rec = dict(q=q, idx=idx, LPdef=3 * q - LPs, lp_integral=integral, cg_iters=it, support=int((x > 1e-9).sum()))
    if integral:
        rec['IPdef'] = 3 * q - float(vf[cols[x > 0.5]].sum())
    elif need_ip:
        T, tries = integer_opt(vf, q, LPs, r, 0.4)
        rec['IPdef'] = 3 * q - T; rec['tries'] = tries
    if 'IPdef' in rec: rec['gap'] = rec['IPdef'] - rec['LPdef']
    rec['secs'] = time.time() - t0
    return rec

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--inst', type=int); ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--workers', type=int, default=2); ap.add_argument('--noip', action='store_true')
    a = ap.parse_args()
    for q in a.q:
        fn = os.path.join(HERE, 'out', f"sum3_q{q:03d}.jsonl")
        with ProcessPoolExecutor(a.workers) as ex, open(fn, 'a') as f:
            for rr in ex.map(run, [(q, i, not a.noip) for i in range(a.start, a.start + a.inst)]):
                f.write(json.dumps(rr) + '\n'); f.flush()
        print('done', q, flush=True)
