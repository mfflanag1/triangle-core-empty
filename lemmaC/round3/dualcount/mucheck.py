"""Sanity check of P(no team blocks | flat prices) against exp(-mu) and exp(-mu+Delta). Evidence only.
Flat prices: a team blocks iff its deficit x_ij + y_ik + z_jk < t. mu = q^3 t^3/6, Delta = (3 q^4/2) t^5/20 (pairs sharing one pair)."""
import numpy as np
rng = np.random.default_rng(20261111)
for q, mu in ((20, 2.0), (40, 2.0), (40, 4.0)):
    t = (6 * mu / q ** 3) ** (1 / 3); n = 20000 if q == 20 else 6000; hit = 0
    for _ in range(n):
        X, Y, Z = rng.random((q, q)), rng.random((q, q)), rng.random((q, q))
        hit += (X[:, :, None] + Y[:, None, :] + Z[None, :, :]).min() >= t
    p = hit / n; Delta = 1.5 * q ** 4 * t ** 5 / 20
    print(dict(q=q, mu=mu, t=round(t, 4), n=n, P_emp=round(p, 4), se=round(float(np.sqrt(p * (1 - p) / n)), 4), exp_minus_mu=round(float(np.exp(-mu)), 4),
               harris_lower=round(float((1 - t ** 3 / 6) ** (q ** 3)), 4), janson_upper=round(float(np.exp(-mu + Delta)), 4)), flush=True)
