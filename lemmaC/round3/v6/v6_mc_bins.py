"""V6: Monte Carlo check of the tilting inputs of section 2.1 at tiny q (true cost law xi = q^(2/3)(U+U'+U'')):
P(bin b) = lam_b/q^2, E[exp(-theta xi) | bin b] = kappa_b/lam_b, E[exp(-theta xi); high] <= J/q^2."""
import numpy as np, math
rng = np.random.default_rng(3); n = 20_000_000
for q, K, B, th in [(3, 2.0, 2, 1.67), (4, 2.5, 5, 1.67), (30, 5.0, 5, 1.67)]:
    xi = q ** (2 / 3) * rng.random((3, n)).sum(0); eta = K / B
    G = lambda x: -math.exp(-th * x) * (x * x / (2 * th) + x / th ** 2 + 1 / th ** 3)
    out = []
    for b in range(1, B + 1):
        sel = (xi > (b - 1) * eta) & (xi <= b * eta); lam = eta ** 3 * (3 * b * b - 3 * b + 1) / 6; kap = G(b * eta) - G((b - 1) * eta)
        out.append("b=%d P*q^2/lam=%.4f LT/(kap/lam)=%.4f" % (b, sel.mean() * q * q / lam, np.exp(-th * xi[sel]).mean() / (kap / lam)))
    hi = xi > K
    print("q=%d K=%g:" % (q, K), "; ".join(out), "; high: E[e^-th xi;high]*q^2/J = %.4f (must be <= 1)" % ((np.exp(-th * xi) * hi).mean() * q * q / (-G(K))))
