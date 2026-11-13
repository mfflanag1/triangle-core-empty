"""V2 referee, front 4: rigorous K-truncated version.  Condition only on stars of teams with cost <= K;
teams above K are handled by the plain union bound (mixed first moment):
  Sigma_K = 3 H' - 2 int_{<=K} h(rho) dLambda + 3 eps ln eps - 2 eps - int_{>K} rho ln rho dLambda,
  H' = E ln(W + w_inf) - int_{<=K} rho ln w dLambda - eps ln w_inf,  eps = E[w_inf/(W + w_inf)] = int_{>K} rho dLambda.
K = 0 is the plain first moment (2.1496); K -> infinity should recover 2.3025."""
import numpy as np, torch, math
from scipy.optimize import minimize
torch.set_default_dtype(torch.float64)
Xmax, nx, nt = 14.0, 280, 1000
xi = (torch.arange(nx) + 0.5) * Xmax / nx; lam = xi ** 2 / 2 * (Xmax / nx)
sg = torch.linspace(-30.0, 55.0, nt); t = torch.exp(sg); ds = (sg[1] - sg[0]).item()

def parts(x, beta, K):
    lo = xi <= K; u = x[:-1][lo]; uinf = x[-1]
    mx = torch.maximum(u.max(), uinf).detach() if lo.any() else uinf.detach()
    u = u - mx; uinf = uinf - mx
    w = torch.exp(u); winf = torch.exp(uinf); l = lam[lo]; x_ = xi[lo]
    tw = t[:, None] * w[None, :]
    phi = torch.exp(-((-torch.expm1(-tw)) @ l)); phip = phi * torch.exp(-t * winf)
    rho = ((tw * torch.exp(-tw)) * phip[:, None]).sum(0) * ds
    eps = (t * winf * torch.exp(-t * winf) * phi).sum() * ds
    ElnW = (torch.exp(-t) - phip).sum() * ds
    rho = rho.clamp(1e-300, 1 - 1e-13); eps = eps.clamp(1e-300, 1.0)
    H = ElnW - (u * rho * l).sum() - eps * uinf
    h = -(rho * torch.log(rho) + (1 - rho) * torch.log1p(-rho))
    hi = ~lo; Z = (lam[hi] * torch.exp(-3 * beta * xi[hi])).sum(); ehi = (lam[hi] * xi[hi] * torch.exp(-3 * beta * xi[hi])).sum() / Z
    L = 3 * H - 2 * (h * l).sum() - 3 * beta * (x_ * rho * l).sum() + eps * (2 * torch.log(eps) - 2 + torch.log(Z))
    e = (x_ * rho * l).sum() + eps * ehi
    return L, e, eps, (rho * l).sum() + eps

def solve(beta, K, x0):
    def f(x):
        xt = torch.tensor(x, requires_grad=True); v = -parts(xt, beta, K)[0]; v.backward(); return v.item(), xt.grad.numpy().copy()
    r = minimize(f, x0, jac=True, method="L-BFGS-B", options=dict(maxiter=3000, ftol=1e-15, gtol=1e-10))
    L, e, eps, norm = parts(torch.tensor(r.x), beta, K)
    return r.x, L.item() + 3 * beta * e.item(), e.item(), eps.item(), norm.item()

for K in (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0):
    x = np.concatenate([-1.5 * xi.numpy(), [-1.5 * max(K, 0.5) - 1.0]])
    b0, b1 = 0.50, 0.60
    x, S0, e0, _, _ = solve(b0, K, x); x, S1, e1, eps, norm = solve(b1, K, x)
    for _ in range(7):
        if abs(S1 - S0) < 1e-13: break
        b2 = b1 - S1 * (b1 - b0) / (S1 - S0); b0, S0 = b1, S1; b1 = b2; x, S1, e1, eps, norm = solve(b1, K, x)
    print("K=%.1f  e_1*(K)=%.5f  beta=%.5f  Sigma_K=%+.1e  eps(used above K)=%.2e  norm=%.6f" % (K, e1, b1, S1, eps, norm), flush=True)
