"""V5 attack on Lemmas 5-6 (floating point, evidence only). Independent of patch_check.py.
(a) law of d = min(d'(u), d''(w)) and the pointwise inequalities;
(b) systematic sampling with slacks that are exactly 0 or 1 at some agents, Monte Carlo over U;
(c) end to end with a packing that is a genuine range-1 function of u: is E_w[cost_w(M_R)] the
    cost of a fresh market of the same size (independence of R and w)?"""
import itertools, numpy as np
from scipy import stats
from scipy.optimize import linprog

rng = np.random.default_rng(20261111)

# (a)
for delta in (0.05, 0.5, 0.95):
    u = rng.random(500000); w = rng.random(500000)
    d1 = 1 - (1 - u) ** (1 / (1 - delta)); d2 = 1 - (1 - w) ** (1 / delta); d = np.minimum(d1, d2)
    ks = stats.kstest(d, "uniform")
    print(f"(a) delta={delta}: KS p={ks.pvalue:.3f}  d<=u/(1-delta): {(d <= u/(1-delta)+1e-15).all()}  d<=w/delta: {(d <= w/delta+1e-15).all()}")

def team(D):
    return D[0][:, :, None] + D[1][None, :, :] + D[2][:, None, :]

def local_rule(Du, q, K):
    """range-1 rule: each agent spreads weight exp(-xi) over its teams with xi <= K, x = min of three."""
    xi = q ** (2 / 3) * Du; H = xi <= K; wgt = np.exp(-xi) * H
    ya = wgt / np.maximum(wgt.sum((1, 2), keepdims=True), 1e-300)
    yb = wgt / np.maximum(wgt.sum((0, 2), keepdims=True), 1e-300)
    yc = wgt / np.maximum(wgt.sum((0, 1), keepdims=True), 1e-300)
    return np.minimum(np.minimum(ya, yb), yc) * H

def sample_R(s, U):
    R = []
    for t in s:
        c = np.concatenate([[0.0], np.cumsum(t)])
        cnt = np.maximum(0, np.ceil(c - U))            # number of points U+n, n>=0, below c
        R.append([k for k in range(len(t)) if cnt[k + 1] > cnt[k]])
    return R

def best(cost, R):
    n = len(R[0]); b = np.inf; arg = []
    for pb in itertools.permutations(R[1]):
        for pc in itertools.permutations(R[2]):
            c = sum(cost[R[0][i], pb[i], pc[i]] for i in range(n))
            if c < b: b, arg = c, [(R[0][i], pb[i], pc[i]) for i in range(n)]
    return (0.0, []) if n == 0 else (b, arg)

# (b) + (c)
q, delta, K = 5, 0.3, 3.0
worst_marg = 0.0; ok_all = True; sizes_equal = True
for seed in range(4):
    r2 = np.random.default_rng(seed)
    u = [r2.random((q, q)) for _ in range(3)]
    Du = team(u); x = local_rule(Du, q, K)
    s = [np.clip(1 - x.sum(ax), 0, 1) for ax in ((1, 2), (0, 2), (0, 1))]
    Sig = [t.sum() for t in s]
    assert max(Sig) - min(Sig) < 1e-10
    n_exact1 = sum(int((t > 1 - 1e-12).sum()) for t in s); n_exact0 = sum(int((t < 1e-12).sum()) for t in s)
    # Monte Carlo over U: marginals and equal sizes
    hits = [np.zeros(q) for _ in range(3)]; nU = 20000
    for U in r2.random(nU):
        R = sample_R(s, U)
        sizes_equal &= len(R[0]) == len(R[1]) == len(R[2]) and len(R[0]) in (int(np.floor(Sig[0] + 1e-9)), int(np.ceil(Sig[0] - 1e-9)))
        for i in range(3): hits[i][R[i]] += 1
    err = max(np.abs(h / nU - t).max() for h, t in zip(hits, s))
    worst_marg = max(worst_marg, err)
    # one (U, w): feasibility and pointwise cost bound against the true LP
    w = [r2.random((q, q)) for _ in range(3)]
    d = [np.minimum(1 - (1 - a) ** (1 / (1 - delta)), 1 - (1 - b) ** (1 / delta)) for a, b in zip(u, w)]
    Dw, Dd = team(w), team(d)
    z = np.zeros((q, q, q)); wc = 0.0; nU2 = 60
    for U in (np.arange(nU2) + 0.5) / nU2:           # midpoint rule over U (approximate integral)
        c, M = best(Dw, sample_R(s, U)); wc += c / nU2
        for e in M: z[e] += 1 / nU2
    f = x + z
    loads = np.concatenate([f.sum((1, 2)), f.sum((0, 2)), f.sum((0, 1))])
    A = np.zeros((3 * q, q ** 3)); idx = np.arange(q ** 3).reshape(q, q, q)
    for a in range(q):
        A[a, idx[a].ravel()] = 1; A[q + a, idx[:, a].ravel()] = 1; A[2 * q + a, idx[:, :, a].ravel()] = 1
    L = linprog(Dd.ravel(), A_eq=A, b_eq=np.ones(3 * q), bounds=(0, None), method="highs").fun
    feas = (f / np.maximum(loads.max(), 1)).ravel()   # midpoint rule leaves O(1/nU2) load error; compare loosely
    cost = (Dd * f).sum(); bound = (Du * x).sum() / (1 - delta) + wc / delta
    ok = L <= cost + 0.05 and cost <= bound + 1e-9
    ok_all &= ok
    print(f"(b) seed {seed}: Sigma={Sig[0]:.3f} agents with slack exactly 1: {n_exact1}, exactly 0: {n_exact0}; "
          f"max marginal error {err:.4f} (MC sd ~{0.5/np.sqrt(nU):.4f}); max|load-1|={abs(loads-1).max():.3f}; "
          f"L={L:.3f} cost_d(x+z)={cost:.3f} bound={bound:.3f} ok={ok}")
print("(b) sizes equal in all three sectors and in {floor,ceil} for every U:", sizes_equal, " all ok:", ok_all)

# (c) independence of R and w: mean patch cost given (u,U) versus a fresh N-market
r3 = np.random.default_rng(7)
u = [r3.random((q, q)) for _ in range(3)]
x = local_rule(team(u), q, K); s = [np.clip(1 - x.sum(ax), 0, 1) for ax in ((1, 2), (0, 2), (0, 1))]
R = sample_R(s, 0.37); N = len(R[0]); a = []; b = []
for _ in range(600):
    w = [r3.random((q, q)) for _ in range(3)]
    a.append(best(team(w), R)[0])
    w2 = [r3.random((N, N)) for _ in range(3)]
    b.append(best(team(w2), [list(range(N))] * 3)[0])
print(f"(c) N={N}: mean patch cost on R {np.mean(a):.4f} +- {np.std(a)/np.sqrt(len(a)):.4f}; fresh N-market {np.mean(b):.4f} +- {np.std(b)/np.sqrt(len(b)):.4f}")
