"""sc (I-b) numerics, follow-up: per-n table for the s=0 classes at q=32, K=1 with more samples, to
separate Monte Carlo noise from a violation of Lemma 3 (N empty):
P_PA(|H|=n, s=0) <= exp(Delta + K^3 n/(2q) + K^6/(36q)) * P_IT(|H|=n, s=0).  Evidence only."""
import math, numpy as np
from collections import Counter
from cmp import excess
q, K, trials = 32, 1.0, 150000
rng = np.random.default_rng(777); t = K * q ** (-2 / 3); pa = Counter(); it = Counter()
for _ in range(trials):
    A = rng.random((q, q)); B = rng.random((q, q)); C = rng.random((q, q))
    a, b, c = np.nonzero(A[:, :, None] + B[None, :, :] + C[:, None, :] <= t)
    if excess(a, b, c, q)[0] == 0: pa[len(a)] += 1
for _ in range(4 * trials):
    n = rng.binomial(q ** 3, t ** 3 / 6); idx = rng.choice(q ** 3, size=n, replace=False)
    if excess(idx // (q * q), (idx // q) % q, idx % q, q)[0] == 0: it[n] += 1
Delta = 3 * q ** 3 * (q - 1) / 2 * t ** 5 / 20
print("q=%d K=%.1f PA trials=%d IT trials=%d Delta=%.4f ; IT exact void=%.6f" % (q, K, trials, 4 * trials, Delta, (1 - t ** 3 / 6) ** (q ** 3)))
for n in range(0, 13):
    x, y = pa[n], it[n]
    if x < 20 or y < 20: continue
    r = (x / trials) / (y / (4 * trials)); se = r * math.sqrt(1 / x + 1 / y); bd = math.exp(Delta + K ** 3 * n / (2 * q) + K ** 6 / (36 * q))
    print("  n=%2d PA=%6d IT=%6d ratio=%.3f +- %.3f  bound=%.3f  %s" % (n, x, y, r, se, bd, "ok" if r <= bd else "above by %.1f se" % ((r - bd) / se)))
