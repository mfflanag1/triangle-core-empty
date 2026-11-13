"""Theorem M on random IFR pairs Q0 <=_disp Q1 = Q0 + lam*Qb (all concave in theta), laws kept well inside the grid.
Also a control: a NON-IFR pair (bimodal input), to see whether the hypothesis matters.  Evidence only."""
import numpy as np
from disp_iter import T, slope_at_level
L, h = 24.0, 0.01
n = int(round(2 * L / h)) + 1
x = -L + h * np.arange(n)
th = np.exp(np.linspace(np.log(1e-6), np.log(20.0), 200))
tg = np.exp(np.linspace(np.log(1e-14), np.log(150.0), 40001))
rng = np.random.default_rng(11)
def randG():
    G = np.zeros_like(x); G1 = np.zeros_like(x)
    for _ in range(rng.integers(1, 4)):
        a, b, p = rng.uniform(0.3, 1.0), rng.uniform(-1, 1), rng.uniform(2.5, 4.0)
        u = np.clip(x - b, 0, None); G += a * u ** p; G1 += a * p * u ** (p - 1)
    return G, G1
def quant(G):
    i = np.argmax(G > 0); return np.interp(np.log(tg), np.log(G[i:]), x[i:])
def dens_from_Q(Q):
    assert Q[-1] < L - 1 and np.all(np.diff(Q) >= 0)
    lG = np.interp(x, Q, np.log(tg), left=-800, right=800)
    F = -np.expm1(-np.exp(np.minimum(lG, 700))); f = np.clip(np.gradient(F, h), 0, None); return f / (f.sum() * h)
def slopes(f):
    m = (x * f).sum() * h; fc = np.interp(x + m, x, f, left=0, right=0)
    G, Gp = T(fc, x, h); return slope_at_level(G, Gp, th)
worst = -9; worst_same = 0
for trial in range(60):
    Ga, _ = randG(); Gb, _ = randG(); lam = rng.uniform(0.05, 0.5)
    Qa, Qb = quant(Ga), quant(Gb)
    s0 = slopes(dens_from_Q(Qa)); s1 = slopes(dens_from_Q(Qa + lam * Qb))
    r = (s1 - s0) / s0; worst = max(worst, r.max())
print('IFR pairs (60): max over pairs, theta of (s_1 - s_0)/s_0 = %+.3e  (Theorem M: <= 0)' % worst)
# control: non-IFR input.  Q0 = quantile of a two-bump law (not concave in theta); Q1 = Q0 + lam*theta-linear expansion
cnt = 0; w2 = -9
for trial in range(60):
    d = rng.uniform(2, 5); w = rng.uniform(0.2, 0.8); sc = rng.uniform(0.15, 0.5)
    f = w * np.exp(-(x + d / 2) ** 2 / (2 * sc ** 2)) + (1 - w) * np.exp(-(x - d / 2) ** 2 / (2 * sc ** 2)); f /= f.sum() * h
    F = np.cumsum(f) * h; F = np.clip(F, 1e-15, 1 - 1e-15)
    # expansion phi(x) = x + lam*psi(x), psi increasing: stretch only one bump
    lam = rng.uniform(0.2, 1.0); x0 = rng.uniform(-d / 2, d / 2)
    phi = x + lam * np.clip(x - x0, 0, None) * (rng.random() < 0.5) + lam * np.clip(x - x0, None, 0) * (rng.random() < 0.5) + 0.3 * lam * x
    f1 = np.interp(x, phi, f / np.gradient(phi, h), left=0, right=0); f1 /= f1.sum() * h
    r = (slopes(f1) - slopes(f)) / slopes(f); w2 = max(w2, r.max()); cnt += r.max() > 1e-4
print('non-IFR control (60 two-bump inputs with an expansion): max (s_1-s_0)/s_0 = %+.3e, pairs with a violation: %d' % (w2, cnt))
