"""V6 informational: exact minimum total deficit of a perfect assignment for small q (iid deficits U+U'+U''),
DP over (used sector-2 set, used sector-3 set). Reports mean and quantiles of min / q^(1/3); Theorem I is asymptotic
(liminf >= 2.302 whp), so small-q values below 2.302 are not a contradiction; the trend is what is informative."""
import sys, numpy as np
rng = np.random.default_rng(11)
def solve(c, q):
    cur = {(0, 0): 0.0}
    for i in range(q):
        nxt = {}
        ci = c[i]
        for (s2, s3), v in cur.items():
            for j in range(q):
                if s2 >> j & 1: continue
                cj = ci[j]; t2 = s2 | 1 << j
                for k in range(q):
                    if s3 >> k & 1: continue
                    key = (t2, s3 | 1 << k); val = v + cj[k]
                    if val < nxt.get(key, 1e9): nxt[key] = val
        cur = nxt
    return cur[((1 << q) - 1, (1 << q) - 1)]
for q, reps in [(3, 400), (4, 400), (5, 300), (6, 150), (7, 60), (8, 24)]:
    vals = np.array([solve(rng.random((3, q, q, q)).sum(0).tolist(), q) for _ in range(reps)]) / q ** (1 / 3)
    print("q=%d reps=%d: min/q^(1/3) mean %.3f (se %.3f), 5%% quantile %.3f, min %.3f, frac below 2.302: %.2f" % (q, reps, vals.mean(), vals.std() / np.sqrt(reps), np.quantile(vals, 0.05), vals.min(), (vals < 2.302).mean()))
