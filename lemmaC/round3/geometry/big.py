"""Larger q: gap of the LP restricted to M + teams sharing no pair with M (D111), versus the full gap D.
Column generation adapted from ../../scaling/cg.py (solve_lp with a column mask; T by cg.integer_opt). Floating point; evidence only.
Seeds: SeedSequence([20261111, 3, 2, q, idx])."""
import os, sys
os.environ.setdefault('OMP_NUM_THREADS', '1')
import json, time
import numpy as np
from scipy.optimize import linprog
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling'))
from cg import amat, integer_opt, solve_lp, PTOL

def solve_masked(V, q, allowed, first, k0=3):
    n = 3 * q; vf = V.reshape(-1); maxadd = 20 * q
    Vm = np.where(allowed.reshape(q, q, q), V, -np.inf)
    init = [first]
    for ax in range(3):
        top = np.argsort(-Vm, axis=ax).take(np.arange(k0), axis=ax)
        g = np.meshgrid(*[np.arange(q)] * 2, indexing='ij')
        for t in range(k0):
            tt = top.take(t, axis=ax)
            idx = [None, None, None]; rest = [d for d in range(3) if d != ax]
            idx[ax] = tt; idx[rest[0]] = g[0]; idx[rest[1]] = g[1]
            init.append(np.ravel_multi_index(tuple(idx), (q, q, q)).ravel())
    cols = np.unique(np.concatenate(init)); cols = cols[allowed[cols]]
    while True:
        lp = linprog(-vf[cols], A_eq=amat(cols, q), b_eq=np.ones(n), bounds=(0, None), method='highs-ds')
        assert lp.status == 0, lp.message
        p = -lp.eqlin.marginals
        r = (V - p[:q, None, None] - p[None, q:2 * q, None] - p[None, None, 2 * q:]).reshape(-1)
        r[~allowed] = -np.inf
        pos = np.nonzero(r > PTOL)[0]
        if len(pos) == 0: return -lp.fun
        if len(pos) > maxadd: pos = pos[np.argpartition(-r[pos], maxadd)[:maxadd]]
        cols = np.union1d(cols, pos)

def run(q, idx):
    t0 = time.time()
    rng = np.random.default_rng(np.random.SeedSequence([20261111, 3, 2, q, idx]))
    A, B, C = rng.random((q, q)), rng.random((q, q)), rng.random((q, q))
    V = A[:, :, None] + B[:, None, :] + C[None, :, :]; vf = V.reshape(-1)
    LPs, p, r, cols, x, it = solve_lp(V, q)
    if np.max(np.minimum(x, 1 - x)) < 1e-6: T = LPs; Mi = cols[x > 0.5]
    else:
        T, _ = integer_opt(vf, q, LPs, r, 0.15 if q <= 48 else 0.3)
        # recover M: MILP again on the candidate set is hidden in cg.integer_opt, so redo it here cheaply
        from scipy.optimize import milp, LinearConstraint, Bounds
        cand = np.nonzero(r >= -(LPs - T) - 1e-7)[0]
        mi = milp(-vf[cand], constraints=LinearConstraint(amat(cand, q), np.ones(3 * q), np.ones(3 * q)), integrality=np.ones(len(cand)), bounds=Bounds(0, 1), options=dict(mip_rel_gap=1e-12))
        assert abs(-mi.fun - T) < 1e-7; Mi = cand[mi.x > 0.5]
    ma, mb, mc = np.unravel_index(Mi, (q, q, q)); ta = np.empty(q, int); tb = np.empty(q, int); tc = np.empty(q, int)
    ta[ma] = np.arange(q); tb[mb] = np.arange(q); tc[mc] = np.arange(q)
    a, b, c = np.unravel_index(np.arange(q ** 3), (q, q, q)); i, j, k = ta[a], tb[b], tc[c]
    nd = (i != j).astype(np.int8) + (j != k) + (i != k)
    out = dict(q=q, idx=idx, D=LPs - T, D111=solve_masked(V, q, nd != 2, Mi) - T, D21=solve_masked(V, q, nd <= 2, Mi) - T, secs=time.time() - t0)
    return out

if __name__ == '__main__':
    q = int(sys.argv[1]); nn = int(sys.argv[2]); start = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    with open(os.path.join(HERE, f'big_q{q:02d}.jsonl'), 'a') as f:
        for idx in range(start, start + nn):
            f.write(json.dumps(run(q, idx)) + '\n'); f.flush()
