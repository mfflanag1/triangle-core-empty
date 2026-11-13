"""v7 referee check of Lemma 2 (Janson-type void bound). Evidence only.
(1) full family: P(min team deficit > t) vs exp(-mu+Delta), mu=q^3 t^3/6, Delta=3q^3(q-1)t^5/40.
(2) sub-family = all q^2 teams at one agent: Delta_F = q^2(q-1) t^5/20.
(3) sub-family = q teams through one pair: exact integral vs exp(-q t^3/6 + C(q,2) t^5/20)."""
import math, numpy as np
from scipy.integrate import quad
rng = np.random.default_rng(7)
Ks = (0.6, 0.8, 1.0, 1.3, 1.6, 2.0)
print("(1) full family")
for q, trials in ((8, 400000), (12, 300000), (16, 300000), (24, 120000), (32, 80000)):
    bs = max(1, 20000000 // q ** 3); mins = []
    done = 0
    while done < trials:
        b = min(bs, trials - done); done += b
        A = rng.random((b, q, q)); B = rng.random((b, q, q)); C = rng.random((b, q, q))
        M = (A[:, :, :, None] + B[:, None, :, :]).min(axis=2) + C   # [a,c]
        mins.append(M.min(axis=(1, 2)))
    m = np.concatenate(mins)
    for K in Ks:
        t = K * q ** (-2 / 3)
        if t > 1: continue
        mu = q ** 3 * t ** 3 / 6; D = 3 * q ** 3 * (q - 1) * t ** 5 / 40
        p = (m > t).mean(); se = math.sqrt(p * (1 - p) / trials); bd = math.exp(-mu + D)
        it = (1 - t ** 3 / 6) ** (q ** 3)
        if p * trials < 20: continue
        print("  q=%2d K=%.1f  P_PA(void)=%.5f +- %.5f  bound=%.5f  IT=%.5f  %s" % (q, K, p, se, min(bd, 9.99999), it, "ok" if p <= bd else "VIOLATION by %.1f se" % ((p - bd) / se)))
print("(2) star of one agent (q^2 teams)")
for q, trials in ((8, 2000000), (16, 2000000), (32, 1000000), (64, 300000)):
    bs = max(1, 20000000 // q ** 2); mins = []; done = 0
    while done < trials:
        b = min(bs, trials - done); done += b
        A = rng.random((b, q)); B = rng.random((b, q, q)); C = rng.random((b, q))
        mins.append((A[:, :, None] + B + C[:, None, :]).min(axis=(1, 2)))
    m = np.concatenate(mins)
    for K in (1.0, 1.5, 2.0, 3.0, 4.0):
        t = K * q ** (-2 / 3)
        if t > 1: continue
        mu = q ** 2 * t ** 3 / 6; D = q ** 2 * (q - 1) * t ** 5 / 20
        p = (m > t).mean(); se = math.sqrt(p * (1 - p) / trials); bd = math.exp(-mu + D)
        print("  q=%2d K=%.1f  P(void)=%.5f +- %.5f  bound=%.5f  %s" % (q, K, p, se, bd, "ok" if p <= bd else "VIOLATION by %.1f se" % ((p - bd) / se)))
print("(3) q teams through one pair, exact")
for q in (8, 32, 128, 1024):
    for K in (1.0, 2.0, 4.0):
        t = K * q ** (-2 / 3)
        if t > 1: continue
        ex = (1 - t) + quad(lambda u: (1 - (t - u) ** 2 / 2) ** q, 0, t)[0]
        bd = math.exp(-q * t ** 3 / 6 + q * (q - 1) / 2 * t ** 5 / 20)
        print("  q=%4d K=%.1f exact=%.8f bound=%.8f %s" % (q, K, ex, bd, "ok" if ex <= bd else "VIOLATION"))
