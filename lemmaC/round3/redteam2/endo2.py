"""Bivariate population dynamics for X = min_k(xi_k - X_j - X_l), K-truncated Poisson hypertree.
Two copies share tree and costs (same Poisson draws, same child indices); they start independent.
Each generation both marginals are re-centred to the fixed-point mean (removes the -2 shift mode only).
Endogenous <=> the pair converges to the diagonal: within-variance -> 0, P(same team) -> 1."""
import numpy as np, sys
sys.argv = [sys.argv[0], sys.argv[1], '0']
exec(open('endo.py').read().split("rng = np.random.default_rng(5)")[0])
N = int(float(sys.argv[1] if False else 2e5)); rng = np.random.default_rng(9)
h = x[1] - x[0]; mu = (x * f).sum() * h
X = np.interp(rng.random((N, 2)), cdf, x)
lam = K ** 3 / 6
for gen in range(1, 41):
    deg = rng.poisson(lam, N); deg[deg == 0] = 1; M = int(deg.sum())
    xi = K * rng.random(M) ** (1 / 3); j = rng.integers(0, N, M); l = rng.integers(0, N, M)
    val = xi[:, None] - X[j] - X[l]
    starts = np.cumsum(deg) - deg
    Xn = np.minimum.reduceat(val, starts, axis=0)
    par = np.repeat(np.arange(N), deg)
    ismin = val == Xn[par]
    same = (ismin[:, 0] & ismin[:, 1]).sum() / N
    cost = (xi * ismin[:, 0]).sum() / N
    shift = Xn.mean(0) - mu; Xn -= shift; X = Xn
    Z = X[:, 0] - X[:, 1]
    if gen <= 12 or gen % 4 == 0:
        print('gen=%2d  sd(X)=%.3f  shift=(%+.4f,%+.4f)  within_var=%.4f  corr=%.4f  med|Z|=%.3f  P(|Z|<0.01)=%.3f  P(same team)=%.4f  cost=%.4f'
              % (gen, X[:, 0].std(), shift[0], shift[1], Z.var() / 2, 1 - Z.var() / (2 * X[:, 0].var()), np.median(np.abs(Z)), (np.abs(Z) < 0.01).mean(), same, cost), flush=True)
