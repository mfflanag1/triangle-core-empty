"""Check of the theta-ray condition (necessary for LP-optimality of M):
   delta_t - delta_s <= dist_pi(s,t) + dist_kappa(s,t) + dist_rho(t,s)   for all s != t,
with pi, kappa, rho the pair-keeping arc costs.  Must hold whenever D = 0 (in fact whenever D21 = 0).  Float; q = 8, seeds of clone.py."""
import json, numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse.csgraph import floyd_warshall
from clone import amat
q = 8; Ae = amat(q); n = 3 * q; fam = {r['idx']: r for r in map(json.loads, open('fam_q08.jsonl'))}
viol_when_D0 = 0; viol_when_D21pos = 0; nD21pos = 0; nviol = 0; worst = 0.0
for idx in range(300):
    rng = np.random.default_rng(np.random.SeedSequence([20261111, 3, 1, q, idx]))
    A, B, C = rng.random((q, q)), rng.random((q, q)), rng.random((q, q))
    vf = (A[:, :, None] + B[None, :, :] + C[:, None, :]).reshape(-1)
    mi = milp(-vf, constraints=LinearConstraint(Ae, np.ones(n), np.ones(n)), integrality=np.ones(q ** 3), bounds=Bounds(0, 1), options=dict(mip_rel_gap=0.0))
    ma, mb, mc = np.unravel_index(np.nonzero(mi.x > 0.5)[0], (q, q, q))
    dA = 1 - A[np.ix_(ma, mb)]; dB = 1 - B[np.ix_(mb, mc)]; dC = 1 - C[np.ix_(ma, mc)]
    g = np.arange(q); delta = dA[g, g] + dB[g, g] + dC[g, g]
    pi = dB + dC - (dB[g, g] + dC[g, g])[:, None]
    ka = dA + dB.T - (dA[g, g] + dB[g, g])[:, None]
    rho = dA + dC - (dA[g, g] + dC[g, g])[None, :]
    for X in (pi, ka, rho): np.fill_diagonal(X, 0)
    dp, dk, dr = floyd_warshall(pi), floyd_warshall(ka), floyd_warshall(rho)
    slack = dp + dk + dr.T - (delta[None, :] - delta[:, None])      # [s, t]
    np.fill_diagonal(slack, 1); v = slack.min() < -1e-9
    r = fam[idx]; nviol += v; worst = min(worst, slack.min())
    if r['D21'] <= 1e-7 and v: viol_when_D0 += 1
    if r['D21'] > 1e-7: nD21pos += 1; viol_when_D21pos += v
print(f'q={q} n=300: theta condition violated in {nviol} instances; violated although D21 = 0: {viol_when_D0} (must be 0); '
      f'instances with D21 > 0: {nD21pos}, of which theta-violated: {viol_when_D21pos}')
