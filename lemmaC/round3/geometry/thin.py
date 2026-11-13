"""Thinning test: keep each pair-keeping team (off-team sharing a pair with M) independently with probability phi; all other teams kept.
D(phi) = LP(restricted) - T.  phi = 1 is the full gap, phi = 0 is D111.  Floating point; evidence only.  Seeds as big.py (nested thinning)."""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
from big import solve_masked
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling'))
from cg import amat, integer_opt, solve_lp
from scipy.optimize import milp, LinearConstraint, Bounds
PHI = (0.5, 0.25, 0.125)
def run(q, idx):
    rng = np.random.default_rng(np.random.SeedSequence([20261111, 3, 2, q, idx]))
    A, B, C = rng.random((q, q)), rng.random((q, q)), rng.random((q, q))
    V = A[:, :, None] + B[:, None, :] + C[None, :, :]; vf = V.reshape(-1)
    LPs, p, r, cols, x, it = solve_lp(V, q)
    if np.max(np.minimum(x, 1 - x)) < 1e-6: T = LPs; Mi = cols[x > 0.5]
    else:
        T, _ = integer_opt(vf, q, LPs, r, 0.15 if q <= 48 else 0.3)
        cand = np.nonzero(r >= -(LPs - T) - 1e-7)[0]
        mi = milp(-vf[cand], constraints=LinearConstraint(amat(cand, q), np.ones(3 * q), np.ones(3 * q)), integrality=np.ones(len(cand)), bounds=Bounds(0, 1), options=dict(mip_rel_gap=1e-12))
        assert abs(-mi.fun - T) < 1e-7; Mi = cand[mi.x > 0.5]
    ma, mb, mc = np.unravel_index(Mi, (q, q, q)); ta = np.empty(q, int); tb = np.empty(q, int); tc = np.empty(q, int)
    ta[ma] = np.arange(q); tb[mb] = np.arange(q); tc[mc] = np.arange(q)
    a, b, c = np.unravel_index(np.arange(q ** 3), (q, q, q)); i, j, k = ta[a], tb[b], tc[c]
    nd = (i != j).astype(np.int8) + (j != k) + (i != k)
    u = np.random.default_rng(np.random.SeedSequence([20261111, 3, 9, q, idx])).random(q ** 3)
    out = dict(q=q, idx=idx, D=LPs - T)
    for phi in PHI: out[f'phi{phi}'] = solve_masked(V, q, (nd != 2) | (u < phi), Mi) - T
    return out
if __name__ == '__main__':
    q = int(sys.argv[1]); nn = int(sys.argv[2])
    with open(os.path.join(HERE, f'thin_q{q:02d}.jsonl'), 'a') as f:
        for idx in range(nn): f.write(json.dumps(run(q, idx)) + '\n'); f.flush()
