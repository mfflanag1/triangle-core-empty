"""Refine the zero of Sigma and test that the stationary point is a maximum of the Lagrangian
L = Sigma - 3*beta*e over choice weights w (random and structured perturbations, other families)."""
import numpy as np
from star_first_moment import xi, lam, stats, solve

def L_of(lnw, beta):
    rho, ElnW = stats(lnw); rho = np.clip(rho, 1e-300, 1 - 1e-12)
    H = ElnW - (lnw * rho * lam).sum()
    h = -(rho * np.log(rho) + (1 - rho) * np.log1p(-rho))
    Sig = 3 * H - 2 * (h * lam).sum(); e = (xi * rho * lam).sum()
    return Sig - 3 * beta * e, Sig, e

# refine zero by secant on beta
b0, b1 = 0.54, 0.58
r0, r1 = solve(b0), solve(b1)
for _ in range(4):
    b2 = b1 - r1["Sigma"] * (b1 - b0) / (r1["Sigma"] - r0["Sigma"])
    b0, r0, b1, r1 = b1, r1, b2, solve(b2)
print("zero: beta=%.5f e=%.5f Sigma=%.2e ann_at_e=%.4f rho0=%.4f" % (b1, r1["e"], r1["Sigma"], r1["ann_at_e"], r1["rho0"]), flush=True)
beta = b1
# reconstruct lnw at the solution
lnw = -3 * beta * xi
for it in range(300):
    rho, _ = stats(lnw); rho = np.clip(rho, 1e-300, 1 - 1e-12)
    new = -beta * xi + (2 / 3) * (np.log(rho) - np.log1p(-rho)); new -= new.max()
    lnw = 0.5 * new + 0.5 * lnw
Ls, Sg, e = L_of(lnw, beta)
print("L*=%.6f Sigma=%.2e e=%.5f" % (Ls, Sg, e), flush=True)
rng = np.random.default_rng(1)
shapes = {"xi": xi, "xi2": xi**2 / 4, "bump1": np.exp(-(xi - 1) ** 2), "bump2": np.exp(-(xi - 2.5) ** 2),
          "bump0": np.exp(-xi**2 * 4), "step": (xi < 1.5) * 1.0, "sin": np.sin(2 * xi), "rand": np.interp(xi, np.linspace(0, 14, 15), rng.normal(size=15))}
worst = -1e9
for name, g in shapes.items():
    for eps in (-1.0, -0.3, -0.05, 0.05, 0.3, 1.0):
        Lp, _, _ = L_of(lnw + eps * g, beta); worst = max(worst, Lp - Ls)
        print("pert %-6s eps=%+.2f  L-L*=%+.3e" % (name, eps, Lp - Ls), flush=True)
for b in (0.6, 1.0, 1.5, 2.0, 3.0, 5.0):
    for fam, lw in (("exp", -b * xi), ("gauss", -b * xi**2 / 2)):
        Lp, _, _ = L_of(lw, beta); worst = max(worst, Lp - Ls)
        print("family %-5s b=%.1f  L-L*=%+.3e" % (fam, b, Lp - Ls), flush=True)
print("max over tests of L-L* = %+.3e" % worst)
