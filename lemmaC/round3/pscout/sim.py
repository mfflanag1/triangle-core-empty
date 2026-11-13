"""Poissonised 3-sector model: LP and IP optima. Usage: sim.py q K nseeds seed0 tlim nproc out"""
import numpy as np, sys, json, time
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import csc_matrix
from multiprocessing import Pool

KGEN = 8.0  # instances generated to cost KGEN, truncated to K (so K=6 and K=8 runs with equal seed are coupled)

def gen(q, K, seed):
    rng = np.random.default_rng([seed, q])
    n = rng.poisson(q * KGEN**3 / 6)
    a = rng.integers(q, size=n); b = rng.integers(q, size=n); c = rng.integers(q, size=n)
    xi = KGEN * rng.random(n) ** (1 / 3)          # density prop. to xi^2 on [0,KGEN]
    key = (a * q + b) * q + c
    order = np.lexsort((xi, key))
    ks = key[order]
    keep = np.ones(n, bool); keep[1:] = ks[1:] != ks[:-1]   # duplicates: keep the cheaper
    idx = order[keep]
    idx = idx[xi[idx] <= K]
    return a[idx], b[idx], c[idx], xi[idx]

def solve(args):
    q, K, seed, tlim = args
    a, b, c, xi = gen(q, K, seed)
    n = len(xi)
    rows = np.concatenate([a, q + b, 2 * q + c]); cols = np.tile(np.arange(n), 3)
    A = csc_matrix((np.ones(3 * n), (rows, cols)), shape=(3 * q, n))
    rec = dict(q=q, K=K, seed=seed, n=n)
    t0 = time.time()
    r = linprog(xi, A_eq=A, b_eq=np.ones(3 * q), bounds=(0, None), method='highs-ds')
    if r.status != 0:
        rec['feas'] = 0; return rec
    x = r.x
    rec.update(feas=1, lp=r.fun / q, lpint=int(np.all((x < 1e-7) | (x > 1 - 1e-7))),
               lpmax=float(xi[x > 1e-9].max()), tlp=round(time.time() - t0, 2))
    if tlim <= 0:
        return rec
    if rec['lpint']:
        rec.update(ub=rec['lp'], lb=rec['lp'], closed=1, ipmax=rec['lpmax'], tip=0.0)
        return rec
    t0 = time.time()
    m = milp(xi, constraints=LinearConstraint(A, 1, 1), integrality=np.ones(n), bounds=Bounds(0, 1),
             options=dict(time_limit=tlim, mip_rel_gap=1e-7))
    rec['tip'] = round(time.time() - t0, 1)
    rec['lb'] = m.mip_dual_bound / q if getattr(m, 'mip_dual_bound', None) is not None else rec['lp']
    if m.x is not None:
        rec.update(ub=m.fun / q, closed=int(m.status == 0), ipmax=float(xi[m.x > 0.5].max()))
    else:
        rec.update(ub=None, closed=0)
    return rec

if __name__ == '__main__':
    q, K, ns, s0, tlim, nproc, out = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6]), sys.argv[7]
    jobs = [(q, K, s, tlim) for s in range(s0, s0 + ns)]
    with open(out, 'a') as f:
        if nproc == 1:
            it = map(solve, jobs)
        else:
            pool = Pool(nproc); it = pool.imap_unordered(solve, jobs)
        for rec in it:
            f.write(json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in rec.items()}) + '\n'); f.flush()
