"""For empty-core instances: the most efficient improving move at the integer optimum M.
Objective v_mu = v + mu * 1[team in M]; mu* = sup{mu : M not optimal}. Just below mu*, the LP
vertex differs from M by an irreducible structure (an edge of the axial polytope at M with the
best gain per unit of disrupted M-weight). Record its number of disrupted M-teams k_eff,
its fractional weights, and the ratio mu*. Floating point; evidence only."""
import json, sys, os, time, argparse
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'gate1'))
from gate1 import instance, structure, ok

def lpsolve(c, Mm, n):
    r = ok(linprog(-c, A_eq=Mm.T.tocsr(), b_eq=np.ones(n), bounds=(0, None), method='highs-ds'), 'LP')
    return r.x, -r.fun

def run(args):
    q, idx = args
    A, B, C = instance('triangle', 'unif', q, idx)
    v, Mm = structure(A, B, C); n = 3 * q
    lam, V = lpsolve(v, Mm, n)
    if np.max(np.minimum(lam, 1 - lam)) < 1e-6:
        return dict(q=q, idx=idx, empty=False)
    mi = milp(-v, constraints=LinearConstraint(Mm.T.tocsr(), np.ones(n), np.ones(n)), integrality=np.ones(q**3),
              bounds=Bounds(0, 1), options=dict(time_limit=900, mip_rel_gap=1e-12))
    sel = mi.x > 0.5; T = -mi.fun
    if V - T < 1e-7 * T:
        return dict(q=q, idx=idx, empty=False, flag='zero_gap')
    ind = sel.astype(float)
    lo, hi = 0.0, 1.0
    while True:  # find hi with M optimal
        x, _ = lpsolve(v + hi * ind, Mm, n)
        if np.max(np.abs(x - ind)) < 1e-6: break
        lo, hi = hi, 2 * hi
    for _ in range(40):
        mid = (lo + hi) / 2
        x, _ = lpsolve(v + mid * ind, Mm, n)
        if np.max(np.abs(x - ind)) < 1e-6: hi = mid
        else: lo, xl = mid, x
    x, _ = lpsolve(v + lo * ind, Mm, n)
    dis = 1 - x[sel]
    k = int((dis > 1e-6).sum())
    fr = x[(x > 1e-6) & (x < 1 - 1e-6)]
    inv = sorted(set(round(1 / w, 3) for w in fr))
    gain = float(v @ x - T)
    return dict(q=q, idx=idx, empty=True, gap=float(V - T), mu_star=lo, k_eff=k, disrupted=float(dis.sum()),
                nfrac=int(len(fr)), frac_inv=inv, gain_at_threshold=gain)

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--n', type=int)
    a = ap.parse_args()
    for q in a.q:
        with ProcessPoolExecutor(6) as ex, open(os.path.join(HERE, 'out', f'eff_q{q:02d}.jsonl'), 'w') as f:
            for r in ex.map(run, [(q, i) for i in range(a.n)]):
                f.write(json.dumps(r) + '\n'); f.flush()
        print('done', q, flush=True)
