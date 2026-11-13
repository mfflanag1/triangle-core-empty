"""V6 check of Proposition 4.5 (direction and tightness) in small binned models with exact (truncated) Poisson sums:
   sup_{mu in grad Phi(R^(B+1))} [3H(mu) - 2S0(mu) + <ell,mu>]  <=  Psi0(u0)  for every u0, with equality at the optimum.
H(mu(u)) = Phi(u) - <mu(u),u> (Legendre), so the left side is evaluated without any use of section 4's inequalities."""
import math, itertools, numpy as np
from scipy.optimize import minimize
from scipy.stats import poisson
rng = np.random.default_rng(7)
def model(K, B, th, NMAX=45):
    eta = K / B; b = np.arange(1, B + 1); lam = eta ** 3 * (3 * b * b - 3 * b + 1) / 6
    G = lambda x: -np.exp(-th * x) * (x * x / (2 * th) + x / th ** 2 + 1 / th ** 3)
    kap = G(eta * b) - G(eta * (b - 1)); J = -G(K)
    grids = np.meshgrid(*[np.arange(NMAX)] * B, indexing="ij"); D = np.stack([g.ravel() for g in grids], 1).astype(float)
    P = np.prod([poisson.pmf(D[:, i], lam[i]) for i in range(B)], axis=0)
    return lam, kap, J, D, P
def stats(u, lam, D, P):
    w = np.exp(u); W = D @ w[:-1] + w[-1]
    Phi = P @ np.log(W); rho = P @ (w[:-1][None, :] / (W[:, None] + w[:-1][None, :])); eps = P @ (w[-1] / W)
    mu_direct = np.concatenate([P @ (D * w[:-1] / W[:, None]), [eps]])        # E pi, without Mecke
    return Phi, rho, eps, mu_direct
def Psi0(u, lam, kap, J, D, P):
    Phi, rho, eps, _ = stats(u, lam, D, P)
    return 3 * Phi + (lam * np.log((1 - rho) ** 3 + np.exp(-3 * u[:-1]) * kap / lam * rho ** 3)).sum() - 3 * eps + math.exp(-3 * u[-1]) * J * eps ** 3
def Lfun(u, lam, kap, J, D, P):
    Phi, rho, eps, mu = stats(u, lam, D, P)
    r = mu[:-1] / lam; H = Phi - mu @ u
    S0 = (lam * -(r * np.log(r) + (1 - r) * np.log(1 - r))).sum() - eps * math.log(eps) + eps
    return 3 * H - 2 * S0 + mu[:-1] @ np.log(kap / lam) + eps * math.log(J), np.abs(r - rho).max()
for (K, B, th) in [(2.0, 1, 1.4), (2.0, 2, 1.7), (3.0, 2, 1.2), (3.0, 3, 1.67), (2.4, 3, 2.3)]:
    lam, kap, J, D, P = model(K, B, th, 45 if B < 3 else 30)
    best = minimize(lambda u: Psi0(u, lam, kap, J, D, P), np.zeros(B + 1), method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-13, maxiter=40000, maxfev=40000))
    psimin = best.fun
    sup = max(-minimize(lambda u: -Lfun(u, lam, kap, J, D, P)[0], rng.normal(0, 1.5, B + 1), method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-13, maxiter=40000, maxfev=40000)).fun for _ in range(6))
    worst = -1e9; mecke = 0.0
    for _ in range(300):
        u0 = rng.normal(0, 2.0, B + 1); u1 = rng.normal(0, 2.5, B + 1)
        L1, dm = Lfun(u1, lam, kap, J, D, P); mecke = max(mecke, dm)
        worst = max(worst, L1 - Psi0(u0, lam, kap, J, D, P))
    print("K=%g B=%d theta=%g: sup L-theta*beta = %.9f ; min_u0 Psi0 = %.9f ; diff %.2e ; max over 300 random (mu,u0) of L - Psi0(u0) = %.3e ; Mecke rho check %.1e" % (K, B, th, sup, psimin, sup - psimin, worst, mecke))
