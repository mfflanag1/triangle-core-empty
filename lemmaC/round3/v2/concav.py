"""V2 referee: is Sigma[rho] = 3H - 2 int h(rho) concave?  Hessian in rho is -3 M^+ + 2 R^-1 with
M = E[diag(p) - p p^T] (choice covariance), R = diag(rho(1-rho) lam).  Concave at rho iff
lambda_max(R^-1/2 M R^-1/2) <= 3/2.  (i) two-bin models, exact Poisson sums; (ii) random profiles in the
continuum model, autograd Jacobian."""
import numpy as np, math, torch, sys
from scipy.stats import poisson
from var import Model
torch.set_default_dtype(torch.float64)

print("== two-bin scan: lambda = m (1/R_A + 1/R_B)")
n = np.arange(0, 90); worst = (0,)
for a in (0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20):
    for b in (0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20):
        pa = poisson.pmf(n, a)[:, None]; pb = poisson.pmf(n, b)[None, :]; P = pa * pb
        na = n[:, None].astype(float); nb = n[None, :].astype(float)
        for lw in np.linspace(-10, 10, 41):
            w = math.exp(lw); tot = na + w * nb
            with np.errstate(divide="ignore", invalid="ignore"):
                m = np.where(tot > 0, na * w * nb / tot ** 2, 0.0)
            m = (P * m).sum(); rA = (P / (1 + tot)).sum(); rB = (P * w / (w + tot)).sum()
            lam_ = m * (1 / (a * rA * (1 - rA)) + 1 / (b * rB * (1 - rB)))
            if lam_ > worst[0]: worst = (lam_, a, b, lw, rA, rB)
print("  max lambda over scan = %.4f at a=%g b=%g ln w=%g rhoA=%.3f rhoB=%.3f  (concave iff <= 1.5)" % worst)

print("== continuum model, random (non-stationary) profiles")
M = Model(nx=140, nt=700); xi = M.xi.numpy(); rng = np.random.default_rng(5)
knots = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 6, 9, 14.0]); mx = 0
for trial in range(int(sys.argv[1]) if len(sys.argv) > 1 else 24):
    kind = trial % 4
    if kind == 0: u0 = np.interp(xi, knots, np.sort(rng.uniform(-30, 5, size=len(knots)))[::-1].copy())
    elif kind == 1: u0 = -rng.uniform(0.2, 6) * xi ** rng.choice([1, 2]) + rng.normal(0, 2, size=len(xi))
    elif kind == 2: c = rng.uniform(0.3, 3); u0 = np.where(xi < c, rng.uniform(0, 25), 0) - rng.uniform(0.3, 4) * xi
    else: u0 = np.where(np.floor(xi / rng.uniform(0.1, 1)) % 2 == 0, 0.0, -rng.uniform(2, 30)) - rng.uniform(0.3, 4) * xi
    ut = torch.tensor(u0); S, e, rho = M.SE(ut); rho = rho.detach().numpy()
    J = torch.autograd.functional.jacobian(lambda v: M.stats(v)[1] * M.lam, ut).numpy()
    R = rho * (1 - rho) * M.lam.numpy(); keep = R > 1e-12
    A = J[np.ix_(keep, keep)] / np.sqrt(np.outer(R[keep], R[keep])); l = np.linalg.eigvalsh((A + A.T) / 2)[-1]; mx = max(mx, l)
    print("  trial %2d kind %d: e=%.3f Sigma=%+.3f rho_max=%.4f lambda_max=%.4f" % (trial, kind, e.item(), S.item(), rho.max(), l), flush=True)
print("  max lambda_max over random profiles = %.4f" % mx)
