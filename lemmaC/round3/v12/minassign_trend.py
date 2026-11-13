"""V12 sanity beyond q=3,4: exact minimum assignment deficit by subset DP (q <= 12), and an upper bound on the
minimum by alternating 2-D assignments (q = 20..80).  Reports min / q^(1/3) against 2.302.
Asymptotic statement; no finite-q value can refute it, but a clear downward crossing would be a red flag."""
import sys, time, itertools
import numpy as np
from scipy.optimize import linear_sum_assignment
rng = np.random.default_rng(77)

def exact_min(c):
    """c: (n, q, q, q) costs; returns (n,) exact minimum over perfect assignments (DP over used subsets of sectors 2, 3)."""
    n, q = c.shape[0], c.shape[1]
    pc = np.array([bin(x).count("1") for x in range(1 << q)])
    INF = 1e9
    dp = np.full((n, 1 << q, 1 << q), INF, dtype=np.float32); dp[:, 0, 0] = 0
    for i in range(q):
        lay = np.nonzero(pc == i)[0]
        for j in range(q):
            a = lay[(lay >> j) & 1 == 0]
            for k in range(q):
                b = lay[(lay >> k) & 1 == 0]
                src = dp[:, a[:, None], b[None, :]] + c[:, i, j, k][:, None, None]
                tgt = dp[:, (a | (1 << j))[:, None], (b | (1 << k))[None, :]]
                dp[:, (a | (1 << j))[:, None], (b | (1 << k))[None, :]] = np.minimum(tgt, src)
    return dp[:, -1, -1]

def heuristic(c, restarts=6, sweeps=30):
    q = c.shape[0]; best = np.inf
    I = np.arange(q)
    for r in range(restarts):
        s2 = rng.permutation(q); s3 = rng.permutation(q); cur = c[I, s2, s3].sum()
        for _ in range(sweeps):
            old = cur
            # re-optimise sector 3 given pairs (i, s2[i])
            _, s3 = linear_sum_assignment(c[I, s2, :])
            # re-optimise sector 2 given pairs (i, s3[i])
            _, s2 = linear_sum_assignment(c[I, :, s3])
            # re-optimise sector 1 given pairs (s2, s3): rows = pair index p, columns = agent i
            _, perm = linear_sum_assignment(c[:, s2, s3].T)       # pair p -> agent perm[p]
            n2 = np.empty(q, int); n3 = np.empty(q, int); n2[perm] = s2; n3[perm] = s3; s2, s3 = n2, n3
            cur = c[I, s2, s3].sum()
            if cur >= old - 1e-12: break
        best = min(best, cur)
    return best

if __name__ == "__main__":
    for q, n, chunk in ((4, 4000, 2000), (5, 4000, 2000), (6, 4000, 1000), (7, 2000, 250), (8, 1500, 100), (9, 600, 30), (10, 300, 10), (11, 120, 4), (12, 40, 1)):
        t = time.time(); out = []
        for _ in range(0, n, chunk):
            out.append(exact_min(rng.random((chunk, q, q, q, 3)).sum(axis=4).astype(np.float32)))
        m = np.concatenate(out) / q ** (1 / 3)
        print("exact q=%2d n=%5d: E[min]/q^(1/3) = %.4f (se %.4f)  P(min < 2.302 q^(1/3)) = %.3f   [%.0fs]" % (q, len(m), m.mean(), m.std() / np.sqrt(len(m)), (m < 2.302).mean(), time.time() - t), flush=True)
    for q, n in ((8, 200), (12, 100), (20, 40), (40, 20), (80, 8)):
        t = time.time(); v = np.array([heuristic(rng.random((q, q, q, 3)).sum(axis=3)) for _ in range(n)]) / q ** (1 / 3)
        print("heuristic (upper bound on min) q=%2d n=%3d: mean best/q^(1/3) = %.4f (se %.4f), share below 2.302: %.2f   [%.0fs]" % (q, n, v.mean(), v.std() / np.sqrt(n), (v < 2.302).mean(), time.time() - t), flush=True)
