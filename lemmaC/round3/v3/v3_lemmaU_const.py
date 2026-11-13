"""V3 referee check of Lemma U's Chernoff constants by exact binomial tails, and of the union-bound arithmetic."""
import numpy as np
from scipy.stats import binom
from math import ceil, floor, log, exp
bad_a = 0; n_a = 0; worst_a = np.inf
for n in (50, 200, 1000, 10**4, 10**5, 10**6):
    for p in (0.5, 0.2, 0.05, 0.01, 1e-3, 1e-4, 1e-5):
        if n * p < 1: continue
        ms = sorted(set([1, 2, 3, 5, 10, 30, 100, 1000, floor(1 / p)]))
        for m in ms:
            if m > 1 / p or m > n: continue
            pi = 1 - (1 - p) ** m
            thr = 0.3 * n * p * m                       # event |N(S)| < thr
            lp = binom.logcdf(ceil(thr) - 1, n, pi)
            marg = -0.079 * n * p * m - lp              # want >= 0
            n_a += 1; worst_a = min(worst_a, marg / (n * p * m))
            if marg < 0: bad_a += 1; print('(a) fails', n, p, m, lp, -0.079 * n * p * m)
print('(a) cases', n_a, 'violations', bad_a, 'worst margin per npm', worst_a)
print('KL(0.9 || 1-e^-3) =', 0.9 * log(0.9 / (1 - exp(-3))) + 0.1 * log(0.1 / exp(-3)))
for n in (10**3, 10**4, 10**5, 10**6, 10**7, 10**8):
    for p in (0.1, 1e-2, 1e-3, 1e-5):
        m = ceil(3 / p)
        if m > n: continue
        pi = 1 - (1 - p) ** m
        lp = binom.logcdf(ceil(0.9 * n) - 1, n, pi)
        print('(b)', n, p, 'log P =', lp, ' claimed bound', -0.0209 * n, 'OK' if lp <= -0.0209 * n else 'CLAIM EXCEEDED', ' rate', -lp / n)
# hypotheses at s=q^(-1/2+eps): smallest q with np>=13(K+2)log n and n>=(600/p)log n, K=1
for eps in (0.25, 0.1, 0.05):
    for e10 in np.arange(1, 400, 0.5):
        q = 10.0 ** e10; p = q ** (-1 + 2 * eps) / 2
        if q * p >= 13 * 3 * log(q) and q >= 600 / p * log(q): print('eps', eps, 'hypotheses first hold near q = 1e%.1f' % e10); break
