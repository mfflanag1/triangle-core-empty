"""Recheck of item (c) of split_patch_attack.py with more samples: optimal w-cost on a fixed
subset R (|R_i| = 2) of a q = 5 market versus a fresh 2-agent market. Evidence only."""
import numpy as np
rng = np.random.default_rng(99); n = 200000
def phi2(A, B, C):  # arrays (n,2,2): min over the 4 assignments of a 2-market
    best = np.full(len(A), np.inf)
    for pb in ((0, 1), (1, 0)):
        for pc in ((0, 1), (1, 0)):
            c = sum(A[:, i, pb[i]] + B[:, pb[i], pc[i]] + C[:, i, pc[i]] for i in range(2))
            best = np.minimum(best, c)
    return best
R = ([1, 4], [0, 2], [3, 4])
w = [rng.random((n, 5, 5)) for _ in range(3)]
a = phi2(w[0][:, R[0]][:, :, R[1]], w[1][:, R[1]][:, :, R[2]], w[2][:, R[0]][:, :, R[2]])
b = phi2(*[rng.random((n, 2, 2)) for _ in range(3)])
print(f"on R: {a.mean():.4f} +- {a.std()/n**.5:.4f}   fresh: {b.mean():.4f} +- {b.std()/n**.5:.4f}")
