"""Efron-Stein probe for Z in {LP*-T, mu*}: resample single pair entries and recompute.
Entries are stratified: S1 = pairs used by M or by supp(lam*) (|S1| <= ~6q); S2 = all other pairs.
ES bound estimate: 0.5 * (|S1| E_S1[(dZ)^2] + |S2| E_S2[(dZ)^2]). Floating point; evidence only."""
import json, sys, os, argparse
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
import scipy.sparse as sp
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'gate1')); sys.path.insert(0, HERE)
from gate1 import instance, structure, ok
from mustruct import integer_opt

def solve(A, B, C):
    q = A.shape[0]; v, Mm = structure(A, B, C); n = 3 * q; Ae = Mm.T.tocsr()
    lp = ok(linprog(-v, A_eq=Ae, b_eq=np.ones(n), bounds=(0, None), method='highs-ds'), 'LP')
    lam, LPs = lp.x, -lp.fun; rlp = v - Mm @ (-lp.eqlin.marginals)
    if np.max(np.minimum(lam, 1 - lam)) < 1e-6: sel = lam > 0.5; T = LPs
    else: sel, T, _, _ = integer_opt(v, Mm, n, LPs, rlp)
    sgn = np.where(sel, -1.0, 1.0)
    Ac = sp.vstack([(sp.diags(sgn) @ Mm).T, sp.csr_matrix(sel.astype(float)[None, :])]).tocsr()
    b = np.zeros(n + 1); b[-1] = 1
    cr = ok(linprog(-np.where(sel, -v, v), A_eq=Ac, b_eq=b, bounds=(0, None), method='highs-ds'), 'cone')
    return LPs - T, -cr.fun, lam, sel, LPs, T

def pairs_of(mask, q):
    x, y, z = np.unravel_index(np.nonzero(mask)[0], (q, q, q))
    return {(0, a, b) for a, b in zip(x, y)} | {(1, a, c) for a, c in zip(x, z)} | {(2, b, c) for b, c in zip(y, z)}

def run(args):
    q, idx, K = args
    rng = np.random.default_rng([7, q, idx])
    M3 = list(instance('triangle', 'unif', q, idx))
    g0, m0, lam, sel, L0, T0 = solve(*M3)
    S1 = sorted(pairs_of(sel | (lam > 1e-7), q))
    allp = [(t, a, b) for t in range(3) for a in range(q) for b in range(q)]
    S1set = set(S1); S2 = [p for p in allp if p not in S1set]
    out = {}
    for nm, S in (('S1', S1), ('S2', S2)):
        dg, dm, dL, dT = [], [], [], []
        for j in rng.choice(len(S), K, replace=False):
            t, a, b = S[j]; Mx = [X.copy() for X in M3]; Mx[t][a, b] = rng.random()
            g, m, _, _, L, T = solve(*Mx); dg.append(g - g0); dm.append(m - m0); dL.append(L - L0); dT.append(T - T0)
        out[nm] = dict(size=len(S), msq_gap=float(np.mean(np.square(dg))), msq_mu=float(np.mean(np.square(dm))),
                       frac_nonzero_gap=float(np.mean(np.abs(dg) > 1e-9)), msq_LP=float(np.mean(np.square(dL))), msq_T=float(np.mean(np.square(dT))))
    out['ES_gap'] = 0.5 * sum(out[s]['size'] * out[s]['msq_gap'] for s in ('S1', 'S2'))
    out['ES_T'] = 0.5 * sum(out[s]['size'] * out[s]['msq_T'] for s in ('S1', 'S2'))
    out['ES_LP'] = 0.5 * sum(out[s]['size'] * out[s]['msq_LP'] for s in ('S1', 'S2'))
    out['ES_mu'] = 0.5 * sum(out[s]['size'] * out[s]['msq_mu'] for s in ('S1', 'S2'))
    return dict(q=q, idx=idx, gap=g0, mu=m0, LP=L0, T=T0, **out)

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--n', type=int); ap.add_argument('--K', type=int, default=15); ap.add_argument('--tag', default='influence')
    a = ap.parse_args()
    for q in a.q:
        with ProcessPoolExecutor(4) as ex, open(os.path.join(HERE, 'out', f'{a.tag}_q{q:02d}.jsonl'), 'w') as f:
            for r in ex.map(run, [(q, i, a.K) for i in range(a.n)]):
                f.write(json.dumps(r) + '\n'); f.flush()
        print('done', q, flush=True)
