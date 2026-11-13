"""Cloning (hub) test. Floating point (HiGHS); evidence only.
T(z) = max v.x over integer x >= 0 with load z_t on the three agents of M-team t  (z = 1: T; M = optimal assignment).
Blow-up lemma: D = 0  =>  T(z) = sum_t z_t v(e_t) for every z.  Tested families:
  hub:  z = 1 + e_t  (one M-team cloned),  g_t = T(1+e_t) - T - v(e_t) >= 0,  D >= g_t / 2
  half: z = 2,                               D_half = T(2)/2 - T,             D >= D_half
Seeds: SeedSequence([20261111, 3, 1, q, idx])."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import sys, json, time, argparse
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); TOL = 1e-7

def amat(q):
    idx = np.arange(q ** 3); a, b, c = np.unravel_index(idx, (q, q, q))
    cols = np.stack([a, q + b, 2 * q + c], 1).reshape(-1); rows = np.repeat(idx, 3)
    return sp.csr_matrix((np.ones(3 * q ** 3), (cols, rows)), shape=(3 * q, q ** 3))

def Tz(vf, Ae, rhs, ub):
    mi = milp(-vf, constraints=LinearConstraint(Ae, rhs, rhs), integrality=np.ones(len(vf)), bounds=Bounds(0, ub),
              options=dict(time_limit=600, mip_rel_gap=0.0))
    assert mi.status == 0, mi.message
    return -mi.fun, np.round(mi.x)

def run(args):
    q, idx, half = args; t0 = time.time()
    rng = np.random.default_rng(np.random.SeedSequence([20261111, 3, 1, q, idx]))
    A, B, C = rng.random((q, q)), rng.random((q, q)), rng.random((q, q))
    vf = (A[:, :, None] + B[None, :, :] + C[:, None, :]).reshape(-1); Ae = amat(q); n = 3 * q
    lp = linprog(-vf, A_eq=Ae, b_eq=np.ones(n), bounds=(0, None), method='highs-ds'); LPs = -lp.fun
    T, x = Tz(vf, Ae, np.ones(n), 1)
    Mi = np.nonzero(x > 0.5)[0]; tm = np.array(np.unravel_index(Mi, (q, q, q))).T
    g = []; ks = []
    for t, (a, b, c) in enumerate(tm):
        rhs = np.ones(n); rhs[[a, q + b, 2 * q + c]] = 2
        Tt, xt = Tz(vf, Ae, rhs, 2)
        gt = Tt - T - vf[Mi[t]]; g.append(float(gt))
        if gt > TOL:  # size of the improving structure: M-teams whose multiplicity dropped below z_t
            zt = np.ones(q); zt[t] = 2
            ks.append(int(np.sum(xt[Mi] < zt - 0.5)))
    rec = dict(q=q, idx=idx, D=LPs - T, T=T, g=g, ghub=max(g), nhub=int(sum(x_ > TOL for x_ in g)), khub=ks)
    if half:
        T2, _ = Tz(vf, Ae, 2 * np.ones(n), 2); rec['Dhalf'] = T2 / 2 - T
    rec['secs'] = time.time() - t0
    return rec

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--n', type=int)
    ap.add_argument('--half', action='store_true'); ap.add_argument('--workers', type=int, default=2)
    a = ap.parse_args()
    for q in a.q:
        with ProcessPoolExecutor(a.workers) as ex, open(os.path.join(HERE, f'clone_q{q:02d}.jsonl'), 'a') as f:
            for r in ex.map(run, [(q, i, a.half) for i in range(a.n)]):
                f.write(json.dumps(r) + '\n'); f.flush()
        print('done', q, flush=True)
