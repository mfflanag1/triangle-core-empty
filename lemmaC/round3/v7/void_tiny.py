"""v7: Lemma 2 at tiny q (stress: Delta/mu largest while the bound is still < 1). Evidence only."""
import math, numpy as np
rng = np.random.default_rng(3)
for q in (2, 3, 4, 5):
    n = 2000000
    A = rng.random((n, q, q)); B = rng.random((n, q, q)); C = rng.random((n, q, q))
    m = ((A[:, :, :, None] + B[:, None, :, :]).min(axis=2) + C).min(axis=(1, 2))
    for t in (0.3, 0.5, 0.7, 0.85, 1.0):
        mu = q ** 3 * t ** 3 / 6; D = 3 * q ** 3 * (q - 1) * t ** 5 / 40; bd = math.exp(-mu + D)
        p = (m > t).mean(); se = math.sqrt(p * (1 - p) / n)
        print("q=%d t=%.2f P(void)=%.5f +- %.5f bound=%.5f %s" % (q, t, p, se, bd, "ok" if p <= bd + 3 * se else "VIOLATION"))
