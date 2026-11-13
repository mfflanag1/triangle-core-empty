"""V6 independent evaluation of Psi_0(K,B,theta,u) from Definition 1.1 of THEOREM_I.md.
(a) composite Gauss-Legendre in s = ln t (float64) and in mpmath (30 digits, trapezoid in s at two steps);
(b) direct Monte Carlo over Poisson stars (no Laplace representation at all)."""
import json, sys, math, numpy as np
from fractions import Fraction
import mpmath as mp
F = "/Users/maxf/projects/economics-research/game-theory/coopgpt/average_case_protection_2026-11-09/lemmaC/round3/firstmoment/" + (sys.argv[1] if len(sys.argv) > 1 else "w_K5_e10.json")
d = json.load(open(F)); B = d["B"]; eta = d["eta"]; K = B * eta; th = d["theta"]
u = np.array(d["u"]); uinf = d["uinf"]; w = np.exp(u); winf = math.exp(uinf)
bb = np.arange(1, B + 1)
lam = eta ** 3 * (3 * bb ** 2 - 3 * bb + 1) / 6
# kappa_b by Gauss-Legendre on each bin (independent of the closed-form antiderivative)
xg, wg = np.polynomial.legendre.leggauss(40)
def gl(f, a, b): x = (a + b) / 2 + (b - a) / 2 * xg; return (b - a) / 2 * (wg * f(x)).sum()
kap = np.array([gl(lambda x: np.exp(-th * x) * x * x / 2, (b - 1) * eta, b * eta) for b in bb])
J = sum(gl(lambda x: np.exp(-th * x) * x * x / 2, K + 2 * j, K + 2 * j + 2) for j in range(40))
Jc = math.exp(-th * K) * (K * K / (2 * th) + K / th ** 2 + 1 / th ** 3)
print("lam sum", lam.sum(), "K^3/6", K ** 3 / 6, "| J quad", J, "closed", Jc)
def psi(t): return np.exp(-t * winf - (lam * -np.expm1(-np.outer(t, w))).sum(1))
# (a) Gauss-Legendre composite in s = ln t
S0, S1, NP = -46.0, 10.0, 448
Phi = 0.0; rho = np.zeros(B); eps = 0.0
for j in range(NP):
    a = S0 + (S1 - S0) * j / NP; b = S0 + (S1 - S0) * (j + 1) / NP
    s = (a + b) / 2 + (b - a) / 2 * xg; ww = (b - a) / 2 * wg; t = np.exp(s); p = psi(t)
    # exp(-t) - psi with cancellation-safe form for small t
    g = -np.expm1(-t * winf - (lam * -np.expm1(-np.outer(t, w))).sum(1)) + np.expm1(-t)
    Phi += (ww * g).sum()
    rho += ((ww * t * p)[:, None] * np.exp(-np.outer(t, w))).sum(0) * w
    eps += (ww * t * p).sum() * winf
def Psi0(Phi, rho, eps):
    C = np.exp(-3 * u) * kap / lam; Cinf = math.exp(-3 * uinf) * J
    return 3 * Phi + (lam * np.log((1 - rho) ** 3 + C * rho ** 3)).sum() - 3 * eps + Cinf * eps ** 3
P = Psi0(Phi, rho, eps)
print("GL float64: Phi=%.12f eps0=%.12f rho1=%.12f norm=%.12f" % (Phi, eps, rho[0], (lam * rho).sum() + eps))
print("GL float64: Psi0=%.10f  -Psi0/theta=%.10f" % (P, -P / th))
# mp trapezoid in s, 30 digits, two step sizes
mp.mp.dps = 30
R = lambda x: mp.mpf(Fraction(str(x)).numerator) / Fraction(str(x)).denominator
mth = R(d["theta"]); meta = R(d["eta"]); mK = meta * B
mu_ = [R(v) for v in d["u"]]; muinf = R(d["uinf"]); mw = [mp.exp(v) for v in mu_]; mwinf = mp.exp(muinf)
mlam = [meta ** 3 * (3 * b * b - 3 * b + 1) / 6 for b in range(1, B + 1)]
def Gm(x): return -mp.exp(-mth * x) * (x * x / (2 * mth) + x / mth ** 2 + 1 / mth ** 3)
mkap = [Gm(meta * b) - Gm(meta * (b - 1)) for b in range(1, B + 1)]; mJ = -Gm(mK)
print("kappa closed vs GL max rel diff:", max(abs(float(mkap[i]) / kap[i] - 1) for i in range(B)))
def trap(h):
    Phi = mp.mpf(0); rho = [mp.mpf(0)] * B; eps = mp.mpf(0)
    n0, n1 = int(-60 / h), int(10 / h)
    for k in range(n0, n1 + 1):
        t = mp.exp(k * h); e = [mp.exp(-t * x) for x in mw]
        E = -t * mwinf - sum(l * (1 - x) for l, x in zip(mlam, e)); p = mp.exp(E)
        Phi += (mp.expm1(-t) - mp.expm1(E)) * h
        for i in range(B): rho[i] += mw[i] * t * e[i] * p * h
        eps += mwinf * t * p * h
    tot = 3 * Phi + sum(mlam[i] * mp.log((1 - rho[i]) ** 3 + mp.exp(-3 * mu_[i]) * mkap[i] / mlam[i] * rho[i] ** 3) for i in range(B))
    tot += -3 * eps + mp.exp(-3 * muinf) * mJ * eps ** 3
    return Phi, eps, rho[0], tot
for h in (mp.mpf(1) / 4, mp.mpf(1) / 8):
    Phi_, eps_, r1, tot = trap(h)
    print("mp trapezoid h=%s: Phi=%s eps0=%s rho1=%s Psi0=%s beta=%s" % (h, mp.nstr(Phi_, 18), mp.nstr(eps_, 18), mp.nstr(r1, 18), mp.nstr(tot, 18), mp.nstr(-tot / mth, 14)))
# (b) Monte Carlo straight from Definition 1.1
rng = np.random.default_rng(20261111); NS = int(float(sys.argv[2])) if len(sys.argv) > 2 else 4_000_000; CH = 500_000
sP = sP2 = 0.0; sr = np.zeros(B); se = 0.0; sPsi_terms = []
for c in range(NS // CH):
    D = rng.poisson(lam, size=(CH, B)); W = winf + D @ w
    lw = np.log(W); sP += lw.sum(); sP2 += (lw ** 2).sum()
    sr += (w[None, :] / (W[:, None] + w[None, :])).sum(0); se += (winf / W).sum()
n = (NS // CH) * CH; PhiM = sP / n; sd = math.sqrt(sP2 / n - PhiM ** 2) / math.sqrt(n); rhoM = sr / n; epsM = se / n
PM = Psi0(PhiM, rhoM, epsM)
print("MC n=%d: Phi=%.6f +- %.6f  eps0=%.6f rho1=%.6f norm=%.6f  Psi0=%.5f +- ~%.5f  beta=%.5f" % (n, PhiM, sd, epsM, rhoM[0], (lam * rhoM).sum() + epsM, PM, 3 * sd, -PM / th))
