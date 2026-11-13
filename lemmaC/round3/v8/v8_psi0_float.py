"""V8 referee: independent floating-point evaluation of Psi_0 (K=5,B=50,theta=1.6712, weights w_K5_e10.json).
Different route from certify.py: Gauss-Legendre on a log-t grid via mpmath-free numpy; plus a direct Monte Carlo-free
check of Phi through truncated Poisson is impossible at B=50, so only the Laplace route is used here. Evidence only."""
import json, math, numpy as np
d = json.load(open("../firstmoment/w_K5_e10.json"))
B = d["B"]; eta = d["eta"]; th = d["theta"]; K = B * eta
u = np.array(d["u"]); uinf = d["uinf"]; w = np.exp(u); winf = math.exp(uinf)
b = np.arange(1, B + 1)
lam = eta ** 3 * (3 * b * b - 3 * b + 1) / 6
G = lambda x: -np.exp(-th * x) * (x * x / (2 * th) + x / th ** 2 + 1 / th ** 3)
kap = G(eta * b) - G(eta * (b - 1)); J = -G(K)
# substitution t = exp(s), composite Gauss-Legendre (8 nodes) on s in [ln 1e-14, ln 5e3]
xs, ws = np.polynomial.legendre.leggauss(8)
edges = np.linspace(math.log(1e-14), math.log(5e3), 6001)
mid = (edges[1:] + edges[:-1]) / 2; half = (edges[1:] - edges[:-1]) / 2
s = (mid[:, None] + half[:, None] * xs[None, :]).ravel(); wt = (half[:, None] * ws[None, :]).ravel()
t = np.exp(s)
E = np.exp(-np.outer(t, w))                       # exp(-t w_b)
psi = np.exp(-t * winf - (1 - E) @ lam)
# Phi = int_0^inf (exp(-t) - psi) dt/t = int (exp(-t)-psi) ds
Phi = np.sum(wt * (np.exp(-t) - psi))
rho = (wt * t * psi) @ (E * w[None, :])
eps = np.sum(wt * t * psi) * winf
C = np.exp(-3 * u) * kap / lam; Cinf = math.exp(-3 * uinf) * J
Psi0 = 3 * Phi + lam @ np.log((1 - rho) ** 3 + C * rho ** 3) - 3 * eps + Cinf * eps ** 3
print("Phi = %.10f  eps0 = %.10f  rho_1 = %.10f  sum lam rho + eps = %.10f" % (Phi, eps, rho[0], lam @ rho + eps))
print("Psi0 = %.9f   -Psi0/theta = %.8f   (certificate: Psi0 <= -3.84751328)" % (Psi0, -Psi0 / th))
for bt in (2.30, 2.302, 2.30224): print("  beta=%.5f theta*beta+Psi0 = %.3e" % (bt, th * bt + Psi0))
