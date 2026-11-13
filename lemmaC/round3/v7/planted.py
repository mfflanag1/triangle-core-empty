"""v7 referee check of Lemma 3 in planted form (N = identity, d_N fixed). Evidence only.
Classes (n, j): n=|H^-|, j=s(h^-;N)=3n-#(pairs of h^- outside pairs(N)). One bin [0,K], so
W=(36/K) q^(2/3).  Claim: P(class | d_N) <= exp(Delta+K^3(n+q)/(2q)+K^6/(36q)) W^j P'(class).
Also reports the 'true-factor' version with W replaced by 3/t (touch) to see how loose W is."""
import math, sys, numpy as np
from collections import Counter
def stat(a, b, c, q):
    n = len(a)
    if n == 0: return 0, 0, 0
    ab = np.unique((a * q + b)[a != b]); bc = np.unique((b * q + c)[b != c]); ac = np.unique((a * q + c)[a != c])
    touch = int(((a == b) | (b == c) | (a == c)).sum())
    return n, 3 * n - len(ab) - len(bc) - len(ac), touch
def run(q, K, x_frac, trials, seed):
    rng = np.random.default_rng(seed); t = K * q ** (-2 / 3); x = x_frac * t; p = t ** 3 / 6
    pa = Counter(); it = Counter(); Xpa = 0; Xit = 0; idx = np.arange(q)
    for _ in range(trials):
        A = rng.random((q, q)); B = rng.random((q, q)); C = rng.random((q, q))
        A[idx, idx] = x; B[idx, idx] = x; C[idx, idx] = x
        m = (A[:, :, None] + B[None, :, :] + C[:, None, :]) <= t
        m[idx, idx, idx] = False
        a, b, c = np.nonzero(m); n, j, tc = stat(a, b, c, q); pa[(n, j)] += 1; Xpa += tc
    for _ in range(4 * trials):
        n = rng.binomial(q ** 3 - q, p)
        while True:
            ii = rng.choice(q ** 3, size=n, replace=False); a, b, c = ii // (q * q), (ii // q) % q, ii % q
            if not ((a == b) & (b == c)).any(): break
        n, j, tc = stat(a, b, c, q); it[(n, j)] += 1; Xit += tc
    Delta = 3 * q ** 3 * (q - 1) * t ** 5 / 40; W = 36 / K * q ** (2 / 3)
    print("q=%d K=%.2f d_N=%.2f*t each  trials PA=%d IT=%d  E[X]: PA=%.2f (theory %.2f) IT=%.3f  Delta=%.3f W=%.1f 3/t=%.1f" % (q, K, x_frac, trials, 4 * trials, Xpa / trials, 3 * q * (q - 1) * max(t - x, 0) ** 2 / 2, Xit / (4 * trials), Delta, W, 3 / t))
    worst = -1e9; rows = []
    for (n, j) in sorted(pa):
        u, v = pa[(n, j)], it[(n, j)]
        if u < 25 or v < 25: continue
        r = (u / trials) / (v / (4 * trials)); bd = math.exp(Delta + K ** 3 * (n + q) / (2 * q) + K ** 6 / (36 * q)) * W ** j
        worst = max(worst, math.log(r / bd)); rows.append((n, j, r, bd))
    js = sorted(set(j for (n, j, r, bd) in rows))
    for j in js:
        rr = [(r, n, bd) for (n, jj, r, bd) in rows if jj == j]; r, n, bd = max(rr)
        print("   j=%d: max_n ratio PA/IT = %.3g (at n=%d), bound there = %.3g, (3/t)^j=%.3g" % (j, r, n, bd, (3 / t) ** j))
    # class-free totals by j
    for j in range(0, 6):
        u = sum(v for (n, jj), v in pa.items() if jj == j); v = sum(v for (n, jj), v in it.items() if jj == j)
        if u >= 25 and v >= 25: print("   all n, j=%d: P_PA=%.5f P_IT=%.5f ratio=%.3g" % (j, u / trials, v / (4 * trials), (u / trials) / (v / (4 * trials))))
    print("   max log(ratio/bound) over classes with >=25 counts: %.2f (must be <= 0)" % worst); sys.stdout.flush()
run(16, 0.8, 0.0, 60000, 1)
run(16, 0.8, 0.5, 60000, 2)
run(24, 0.8, 0.0, 30000, 3)
run(16, 1.0, 0.3, 60000, 4)
