"""Sister model (a): random axial 3-index assignment with iid team costs. Floating point; evidence only.
Costs c(a,b,c) iid Exp(1), scaled by n^2 so each agent sees a rate-1 Poisson process of team costs.
Reuses solve_lp / integer_opt from ../../scaling/cg.py (they maximise V; here V = -n^2 * cost).
Records LP, IP (scaled totals), gap = IP_cost - LP_cost >= 0, and whether the LP optimum is integral.
Seeds: SeedSequence([20261111, 3, 1, n, idx])."""
import os, sys, json, time, argparse
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling'))
import cg

def integer_opt(vf, q, LPs, r, delta0):
    """Reduced-cost fixing as in cg.integer_opt, without its delta cap (the iid gap is Theta(n) in scaled units)."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    n = 3 * q; delta = delta0; tries = []
    while True:
        cand = np.nonzero(r >= -delta)[0]
        mi = milp(-vf[cand], constraints=LinearConstraint(cg.amat(cand, q), np.ones(n), np.ones(n)),
                  integrality=np.ones(len(cand)), bounds=Bounds(0, 1), options=dict(time_limit=1800, mip_rel_gap=1e-12))
        tries.append((delta, len(cand), int(mi.status)))
        if mi.status == 0 and -mi.fun >= LPs - delta + 1e-9:
            return -mi.fun, tries
        if delta > 40 * q: raise RuntimeError('delta blew up')
        delta *= 2

def run(args):
    n, idx, need_ip = args
    t0 = time.time()
    rng = np.random.default_rng(np.random.SeedSequence([20261111, 3, 1, n, idx]))
    V = -(n ** 2) * rng.exponential(size=(n, n, n))
    vf = V.reshape(-1)
    LPs, p, r, cols, x, it = cg.solve_lp(V, n)
    integral = bool(np.max(np.minimum(x, 1 - x)) < 1e-6)
    rec = dict(n=n, idx=idx, LP=-LPs, lp_integral=integral, cg_iters=it, lp_secs=time.time() - t0,
               support=int((x > 1e-9).sum()), minpos=float(x[x > 1e-9].min()))
    if integral:
        rec['IP'] = -float(vf[cols[x > 0.5]].sum())
    elif need_ip:
        T, tries = integer_opt(vf, n, LPs, r, 0.3 * n)
        rec['IP'] = -T; rec['tries'] = tries
    if 'IP' in rec: rec['gap'] = rec['IP'] - rec['LP']
    rec['secs'] = time.time() - t0
    return rec

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, nargs='+'); ap.add_argument('--inst', type=int); ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--workers', type=int, default=2); ap.add_argument('--noip', action='store_true')
    a = ap.parse_args()
    for n in a.n:
        fn = os.path.join(HERE, 'out', f"iid_n{n:03d}.jsonl")
        with ProcessPoolExecutor(a.workers) as ex, open(fn, 'a') as f:
            for rr in ex.map(run, [(n, i, not a.noip) for i in range(a.start, a.start + a.inst)]):
                f.write(json.dumps(rr) + '\n'); f.flush()
        print('done', n, flush=True)
