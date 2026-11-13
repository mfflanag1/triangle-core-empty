"""Sanity checks (floating point, evidence only) of the two finite-dimensional facts used in section 4:
 (a) M := E[diag(pi) - pi pi^T] <= Rt := diag(lam_b rho_b (1-rho_b), eps) in the PSD order (B bins + high option);
 (b) the tangent inequality (*) for H - S0 at random base points, against H(mu) computed by direct minimisation."""
import itertools, math, numpy as np
from scipy.optimize import minimize
from scipy.stats import poisson
rng = np.random.default_rng(1); NMAX = 40
def setup(lam):
    B = len(lam); grids = np.meshgrid(*[np.arange(NMAX)] * B, indexing="ij")
    D = np.stack([g.ravel() for g in grids], 1).astype(float)
    P = np.prod([poisson.pmf(D[:, b], lam[b]) for b in range(B)], axis=0)
    return D, P
def stats(u, D, P):
    w = np.exp(u); W = D @ w[:-1] + w[-1]
    pi = np.concatenate([D * w[:-1] / W[:, None], (w[-1] / W)[:, None]], 1)
    mu = P @ pi; M = np.diag(mu) - (pi * P[:, None]).T @ pi; Phi = P @ np.log(W)
    return Phi, mu, M
def S0(mu, lam):
    r = mu[:-1] / lam; e = mu[-1]
    return (lam * -(r * np.log(r) + (1 - r) * np.log(1 - r))).sum() - e * math.log(e) + e
worst_a = -1; worst_b = -1e9; nb = 0
for trial in range(60):
    B = int(rng.choice([1, 2, 3])); lam = rng.uniform(0.05, 2.5, size=B)
    if B == 3: lam = lam * 0.6
    D, P = setup(lam)
    u0 = rng.normal(0, 2.0, size=B + 1); Phi0, mu0, M0 = stats(u0, D, P)
    rho = mu0[:-1] / lam; Rt = np.diag(np.concatenate([lam * rho * (1 - rho), [mu0[-1]]]))
    worst_a = max(worst_a, np.linalg.eigvalsh(np.diag(np.diag(Rt) ** -0.5) @ M0 @ np.diag(np.diag(Rt) ** -0.5))[-1])
    gS = np.concatenate([np.log((1 - rho) / rho), [-math.log(mu0[-1])]])
    rhs0 = Phi0 + (lam * np.log(1 - rho)).sum() - mu0[-1]
    for j in range(6):
        u1 = rng.normal(0, 3.0, size=B + 1); Phi1, mu1, _ = stats(u1, D, P)   # mu1 = grad Phi(u1), so H(mu1) = Phi1 - <mu1,u1>
        lhs = Phi1 - mu1 @ u1 - S0(mu1, lam); rhs = rhs0 - (u0 + gS) @ mu1
        worst_b = max(worst_b, lhs - rhs); nb += 1
print("(a) max over 60 random models of lambda_max(Rt^-1/2 M Rt^-1/2) = %.6f (claim: <= 1)" % worst_a)
print("(b) max over %d random pairs of [H - S0](mu) - tangent bound = %.3e (claim: <= 0)" % (nb, worst_b))
