"""Interpolation between the pair-additive market and independent team values.
legs = number of pair-additive legs (3: the market; 2: A[a,b] + B[b,c] + fresh U per team; 1: A[a,b] + two fresh U; 0: three fresh U).
Team deficits have the same marginal law in all four models.  D = LP* - T by ../../scaling/cg.py.  Float; evidence only.
Seeds: SeedSequence([20261111, 3, 6, q, idx])."""
import os, sys, json
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling'))
from cg import integer_opt, solve_lp
def gap(V, q):
    LPs, p, r, cols, x, it = solve_lp(V, q)
    if np.max(np.minimum(x, 1 - x)) < 1e-6: return 0.0, 3 * q - LPs
    T, _ = integer_opt(V.reshape(-1), q, LPs, r, 0.15 if q <= 48 else 0.3)
    return LPs - T, 3 * q - T
if __name__ == '__main__':
    q = int(sys.argv[1]); nn = int(sys.argv[2]); out = {k: [] for k in range(4)}; dT = {k: [] for k in range(4)}
    for idx in range(nn):
        rng = np.random.default_rng(np.random.SeedSequence([20261111, 3, 6, q, idx]))
        A, B, C = rng.random((q, q)), rng.random((q, q)), rng.random((q, q)); E = rng.random((3, q, q, q))
        legs = [A[:, :, None] + 0 * E[0], B[None, :, :] + 0 * E[0], C[:, None, :] + 0 * E[0]]
        for k in range(4):
            V = sum(legs[:k]) + E[k:].sum(0) if k < 3 else sum(legs)
            d, t = gap(V, q); out[k].append(d); dT[k].append(t)
    json.dump(dict(q=q, n=nn, D=out, defT=dT), open(os.path.join(HERE, f'hybrid_q{q:02d}.json'), 'w'))
    for k in range(4):
        D = np.array(out[k]); print(f"q={q} n={nn} pair-additive legs={k}: P(D>0)={(D > 1e-7).mean():.3f} ED={D.mean():.4f} E(3q-T)/q^(1/3)={np.mean(dT[k]) / q ** (1 / 3):.3f}", flush=True)
