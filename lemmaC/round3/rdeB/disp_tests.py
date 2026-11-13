"""(1) Theorem M on random IFR pairs Q0 <=_disp Q1 = Q0 + lam*Qb (both concave in theta).
(2) Candidate super-solution rho^0: G0(x) = A*E(x+sig*N)_+^3  (IFR; upper tail exp(-A x^3), Gaussian lower tail):
    test T(rho^0) <=_disp rho^0, i.e. s_out(theta) >= s_in(theta), then iterate downwards.
Floating point: evidence only."""
import numpy as np
from scipy.stats import norm
from disp_iter import T, slope_at_level
L, h = 24.0, 0.01
n = int(round(2 * L / h)) + 1
x = -L + h * np.arange(n)
th = np.exp(np.linspace(np.log(1e-8), np.log(30.0), 240))
rng = np.random.default_rng(7)

def dens_from_G(G, G1): return G1 * np.exp(-G)
def randG():
    """random convex nondecreasing G on the grid: sum of a_i (x-b_i)_+^{p_i}."""
    G = np.zeros_like(x); G1 = np.zeros_like(x)
    for _ in range(rng.integers(1, 4)):
        a, b, p = rng.uniform(0.05, 1.0), rng.uniform(-2, 1), rng.uniform(1.0, 4.0)
        u = np.clip(x - b, 0, None); G += a * u ** p; G1 += a * p * u ** (p - 1) * (u > 0)
    return G, G1
def quant(G, tg):
    ok = G > 0; i = np.argmax(ok)
    return np.interp(np.log(tg), np.log(G[i:]), x[i:])
tg = np.exp(np.linspace(np.log(1e-13), np.log(200.0), 20001))
worst = -1
for trial in range(40):
    Ga, G1a = randG(); Gb, G1b = randG(); lam = rng.uniform(0.1, 1.5)
    Qa, Qb = quant(Ga, tg), quant(Gb, tg)
    Q1 = Qa + lam * Qb
    # density of the law with quantile Q1 (in theta): G_1(x) = Q1^{-1}(x)
    G1x = np.exp(np.interp(x, Q1, np.log(tg), left=-700, right=np.log(tg[-1]) + 50))
    F1 = 1 - np.exp(-G1x); f1 = np.gradient(F1, h); f1 = np.clip(f1, 0, None); f1 /= f1.sum() * h
    f0 = dens_from_G(Ga, G1a); f0 /= f0.sum() * h
    out = []
    for f in (f0, f1):
        m = (x * f).sum() * h; fc = np.interp(x + m, x, f, left=0, right=0)
        G, Gp = T(fc, x, h); out.append(slope_at_level(G, Gp, th))
    r = (out[1] - out[0]) / out[0]          # Theorem M predicts <= 0
    worst = max(worst, r.max())
print('Theorem M, 40 random IFR pairs: max over pairs and theta of (s_1 - s_0)/s_0 = %+.3e (prediction <= 0)' % worst)

def Gsup(A, sig):
    z = x / sig
    G = A * sig ** 3 * ((z ** 3 + 3 * z) * norm.cdf(z) + (z ** 2 + 2) * norm.pdf(z))
    G1 = A * 3 * sig ** 2 * ((z ** 2 + 1) * norm.cdf(z) + z * norm.pdf(z))
    return np.maximum(G, 1e-300), np.maximum(G1, 1e-300)
th2 = np.exp(np.linspace(np.log(1e-7), np.log(30.0), 240))
for A, sig in [(1 / 12, 1.5), (1 / 12, 2.0), (1 / 10, 2.0), (1 / 8, 1.5), (1 / 20, 3.0)]:
    G0, G10 = Gsup(A, sig); s_in = slope_at_level(G0, G10, th2)
    f = dens_from_G(G0, G10); f /= f.sum() * h
    m = (x * f).sum() * h; fc = np.interp(x + m, x, f, left=0, right=0)
    G, Gp = T(fc, x, h); s_out = slope_at_level(G, Gp, th2)
    r = s_out / s_in; i = np.argmin(r)
    line = 'A=%.4f sig=%.1f: min s_out/s_in = %.4f at theta=%.2e (need >= 1)' % (A, sig, r.min(), th2[i])
    # iterate downwards, report costs (lower bounds on c if rho^0 is a super-solution)
    cs = []; f = G1p = None; fcur = dens_from_G(G0, G10); fcur /= fcur.sum() * h
    for it in range(8):
        m = (x * fcur).sum() * h; fc = np.interp(x + m, x, fcur, left=0, right=0)
        G, Gp = T(fc, x, h); fcur = Gp * np.exp(-G); cs.append(3 * (x * fcur).sum() * h)
    print(line, '| 3*Gamma(T^k rho^0), k=0..7:', ' '.join('%.4f' % c for c in cs))
