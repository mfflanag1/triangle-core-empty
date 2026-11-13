"""Floating-point search (evidence only): discrete B-bin model with high option.
Minimise Psi0(theta,u)/theta, where the certified bound is beta_cert = -Psi0/theta.
Psi0 = 3 Phi + sum_b lam_b ln((1-rho_b)^3 + C_b rho_b^3) - 3 eps + C_inf eps^3,
C_b = exp(-3u_b) kappa_b/lam_b, C_inf = exp(-3u_inf) J  (tau = 0 by shift invariance)."""
import sys, json, math, numpy as np, torch
from scipy.optimize import minimize
torch.set_default_dtype(torch.float64)
K = float(sys.argv[1]); eta = float(sys.argv[2]); B = int(round(K / eta))
edges = torch.arange(B + 1) * eta; lo, hi = edges[:-1], edges[1:]
lam = (hi ** 3 - lo ** 3) / 6
sg = torch.linspace(-32.0, 40.0, 3000); t = torch.exp(sg); ds = (sg[1] - sg[0]).item()
def G(x, th): return -torch.exp(-th * x) * (x * x / (2 * th) + x / th ** 2 + 1 / th ** 3)
def psi0(x):
    u, uinf, th = x[:B], x[B], x[B + 1]
    kap = G(hi, th) - G(lo, th); J = -G(torch.tensor(K), th)
    w = torch.exp(u); winf = torch.exp(uinf); tw = t[:, None] * w[None, :]
    psi = torch.exp(-t * winf - ((-torch.expm1(-tw)) @ lam))
    Phi = (torch.exp(-t) - psi).sum() * ds
    rho = (tw * torch.exp(-tw) * psi[:, None]).sum(0) * ds
    eps = (t * winf * psi).sum() * ds
    Cb = torch.exp(-3 * u) * kap / lam; Cinf = torch.exp(-3 * uinf) * J
    P = 3 * Phi + (lam * torch.log((1 - rho) ** 3 + Cb * rho ** 3)).sum() - 3 * eps + Cinf * eps ** 3
    return P, rho, eps, Phi
def f(x):
    xt = torch.tensor(x, requires_grad=True); v = psi0(xt)[0] / xt[B + 1]; v.backward(); return v.item(), xt.grad.numpy().copy()
mid = ((lo + hi) / 2).numpy()
x0 = np.concatenate([-0.557 * mid + (2 / 3) * np.log(0.77 * np.exp(-1.3 * mid) / (1 - 0.77 * np.exp(-1.3 * mid))), [-4.6, 1.672]])
r = minimize(f, x0, jac=True, method="L-BFGS-B", options=dict(maxiter=20000, maxfun=40000, ftol=1e-16, gtol=1e-11))
P, rho, eps, Phi = psi0(torch.tensor(r.x))
print("K=%g eta=%g B=%d: beta=%.6f theta=%.6f Phi=%.6f eps=%.5f norm=%.8f rho0=%.4f uinf=%.4f it=%d" % (K, eta, B, -r.fun, r.x[B + 1], Phi.item(), eps.item(), (lam * rho).sum().item() + eps.item(), rho[0].item(), r.x[B], r.nit))
if len(sys.argv) > 3:
    json.dump(dict(K=K, eta=eta, B=B, theta=round(float(r.x[B + 1]), 4), u=[round(float(v), 5) for v in r.x[:B]], uinf=round(float(r.x[B]), 5)), open(sys.argv[3], "w"))
