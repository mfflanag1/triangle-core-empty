"""Local-route check. For each instance (gate1 seeds): integer optimum M (LP vertex if integral,
else scipy milp); N3 = number of 3-subsets S of M whose 9-agent sub-market has an empty core
(LP_S > T_S, exact vertex list of the 3x3x3 image from verts3.py); if the core is empty and
N3 = 0, N4 = number of 4-subsets with empty sub-core (HiGHS LP per subset).
Floating point; evidence only."""
import json, sys, os, itertools, time, argparse
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'gate1'))
from gate1 import instance, structure, ok
VERT = json.load(open(os.path.join(HERE, 'out', 'verts3.json')))
from fractions import Fraction
MARG = np.array([[float(Fraction(x)) if '/' in x and len(x) < 12 else float(x) if '.' in x or 'e' in x else float(Fraction(x)) for x in o['marg']] for o in VERT])
FRAC = np.array([not o['integral'] for o in VERT])
TOL = 1e-9

def sub_vec(A, B, C, X, Y, Z):
    return np.concatenate([A[np.ix_(X, Y)].ravel(), B[np.ix_(X, Z)].ravel(), C[np.ix_(Y, Z)].ravel()])

def sub_lp(A, B, C, X, Y, Z):
    k = len(X)
    v, Mm = structure(A[np.ix_(X, Y)], B[np.ix_(X, Z)], C[np.ix_(Y, Z)])
    r = ok(linprog(-v, A_eq=Mm.T.tocsr(), b_eq=np.ones(3 * k), bounds=(0, None), method='highs-ds'), 'sub')
    return -r.fun

def run(args):
    q, idx, do4 = args
    t0 = time.time()
    A, B, C = instance('triangle', 'unif', q, idx)
    v, Mm = structure(A, B, C)
    n = 3 * q
    lp = ok(linprog(-v, A_ub=Mm.T.tocsr(), b_ub=np.ones(n), bounds=(0, None), method='highs-ds'), 'LP')
    lam, V = lp.x, -lp.fun
    integral = np.max(np.minimum(lam, 1 - lam)) < 1e-6
    if integral:
        sel = lam > 0.5; T = V
    else:
        mi = milp(-v, constraints=LinearConstraint(Mm.T.tocsr(), -np.inf, np.ones(n)), integrality=np.ones(q**3),
                  bounds=Bounds(0, 1), options=dict(time_limit=900, mip_rel_gap=1e-12))
        if mi.status != 0:
            return dict(q=q, idx=idx, fail=True)
        sel = mi.x > 0.5; T = -mi.fun
    tm = np.array(np.unravel_index(np.nonzero(sel)[0], (q, q, q))).T  # rows: (x,y,z)
    empty = (V - T) > 1e-7 * T
    # N3
    W = []
    trip = list(itertools.combinations(range(q), 3))
    for S in trip:
        X, Y, Z = tm[list(S), 0], tm[list(S), 1], tm[list(S), 2]
        W.append(sub_vec(A, B, C, X, Y, Z))
    W = np.array(W)
    vals = W @ MARG.T
    TS = np.array([v.reshape(q, q, q)[tm[list(S), 0], tm[list(S), 1], tm[list(S), 2]].sum() for S in trip])
    excess = vals[:, FRAC].max(axis=1) - TS
    n3 = int((excess > TOL).sum())
    intmax_violation = float((vals[:, ~FRAC].max(axis=1) - TS).max())  # must be <= 0 (M optimal)
    # S*: M-teams with some agent not sitting in its own M-team at weight 1 under the LP vertex
    L = lam.reshape(q, q, q)
    own = np.array([L[x, y, z] for x, y, z in tm])
    Sstar = [s for s in range(q) if own[s] < 1 - 1e-6]
    fr = lam[(lam > 1e-6) & (lam < 1 - 1e-6)]
    dens = sorted(set(int(round(1 / w)) if abs(1 / w - round(1 / w)) < 1e-6 else -1 for w in fr))
    nfrac = int(len(fr))
    rec = dict(q=q, idx=idx, kstar=len(Sstar), nfrac=nfrac, frac_inv=dens, V=V, T=T, empty=bool(empty), n3=n3, max_excess3=float(excess.max()),
               intmax_violation=intmax_violation, gap=V - T)
    if empty and n3 == 0 and do4 and q <= 12:
        n4 = 0; best = -1
        for S in itertools.combinations(range(q), 4):
            S = list(S)
            X, Y, Z = tm[S, 0], tm[S, 1], tm[S, 2]
            e = sub_lp(A, B, C, X, Y, Z) - v.reshape(q, q, q)[X, Y, Z].sum()
            best = max(best, e)
            if e > TOL: n4 += 1
        rec.update(n4=n4, max_excess4=best)
    rec['secs'] = time.time() - t0
    return rec

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--n', type=int)
    ap.add_argument('--do4', type=int, default=1); ap.add_argument('--workers', type=int, default=6)
    a = ap.parse_args()
    for q in a.q:
        with ProcessPoolExecutor(a.workers) as ex, open(os.path.join(HERE, 'out', f'local_q{q:02d}.jsonl'), 'w') as f:
            for r in ex.map(run, [(q, i, a.do4) for i in range(a.n)]):
                f.write(json.dumps(r) + '\n'); f.flush()
        print('done', q, flush=True)
