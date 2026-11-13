"""Sanity check (floating point, evidence only) of the patching step:
split coupling, equal sector slack sums, systematic sampling, loads of x+z, cost bound."""
import itertools, numpy as np
from scipy.optimize import linprog

def team(D):  # D = [AB, BC, AC]
    return D[0][:, :, None] + D[1][None, :, :] + D[2][:, None, :]

def best_assignment(cost, R):
    n = len(R[0]); best, arg = np.inf, None
    for pb in itertools.permutations(R[1]):
        for pc in itertools.permutations(R[2]):
            c = sum(cost[R[0][i], pb[i], pc[i]] for i in range(n))
            if c < best: best, arg = c, [(R[0][i], pb[i], pc[i]) for i in range(n)]
    return best, arg

def run(seed, q=5, delta=0.25):
    rng = np.random.default_rng(seed)
    u = [rng.random((q, q)) for _ in range(3)]; w = [rng.random((q, q)) for _ in range(3)]
    d1 = [1 - (1 - a) ** (1 / (1 - delta)) for a in u]; d2 = [1 - (1 - a) ** (1 / delta) for a in w]
    d = [np.minimum(a, b) for a, b in zip(d1, d2)]
    assert all((d[i] <= u[i] / (1 - delta) + 1e-14).all() and (d[i] <= w[i] / delta + 1e-14).all() for i in range(3))
    Du, Dw, Dd = team(u), team(w), team(d)
    x = rng.random((q, q, q)) * (rng.random((q, q, q)) < 0.08)      # a packing depending on u-side randomness only
    x /= 1.2 * max(1.0, x.sum((1, 2)).max(), x.sum((0, 2)).max(), x.sum((0, 1)).max())
    s = [1 - x.sum((1, 2)), 1 - x.sum((0, 2)), 1 - x.sum((0, 1))]
    Sig = [t.sum() for t in s]
    assert max(Sig) - min(Sig) < 1e-12                               # equal sector sums
    cum = [np.concatenate([[0.0], np.cumsum(t)]) for t in s]
    bps = sorted(set([0.0, 1.0] + [float(c % 1) for cc in cum for c in cc]))
    z = np.zeros((q, q, q)); wcost = 0.0; EN = 0.0
    for lo, hi in zip(bps[:-1], bps[1:]):
        if hi - lo < 1e-15: continue
        U = (lo + hi) / 2
        cnt = [np.maximum(0, np.ceil(cc - U)) for cc in cum]
        R = [[k for k in range(q) if cnt[i][k + 1] > cnt[i][k]] for i in range(3)]
        assert len(R[0]) == len(R[1]) == len(R[2]) and len(R[0]) in (int(np.floor(Sig[0] + 1e-12)), int(np.ceil(Sig[0] - 1e-12)))
        EN += (hi - lo) * len(R[0])
        if R[0]:
            c, M = best_assignment(Dw, R); wcost += (hi - lo) * c
            for e in M: z[e] += hi - lo
    f = x + z
    loads = np.concatenate([f.sum((1, 2)), f.sum((0, 2)), f.sum((0, 1))])
    A = np.zeros((3 * q, q ** 3)); idx = np.arange(q ** 3).reshape(q, q, q)
    for a in range(q):
        A[a, idx[a].ravel()] = 1; A[q + a, idx[:, a].ravel()] = 1; A[2 * q + a, idx[:, :, a].ravel()] = 1
    L = linprog(Dd.ravel(), A_eq=A, b_eq=np.ones(3 * q), bounds=(0, None), method="highs").fun
    cost = (Dd * f).sum(); bound = (Du * x).sum() / (1 - delta) + wcost / delta
    return abs(loads - 1).max(), abs(EN - Sig[0]), L, cost, bound

for seed in range(6):
    le, en, L, cost, bound = run(seed)
    print(f"seed {seed}: max|load-1|={le:.1e} |EN-Sigma|={en:.1e}  L_q={L:.4f} <= cost(x+z)={cost:.4f} <= bound={bound:.4f}  ok={L <= cost + 1e-9 <= bound + 2e-9}")
