"""Referee check of Theorem A at r = 1 (floating point, evidence only).

Root o in sector A. Teams through o: (o,b,c), deficit a_b + c_c + B_bc. Only b with a_b <= s and
c with c_c <= s can be near-tight, so sample those directly (exact law).
1-ball is a tree iff no two near-tight teams through o share b or share c.
Compare mu(tree with n teams) with Poisson(lambda) and with the per-shape lower bound (2.2).
"""
import sys, math
import numpy as np

def run(q, K, N, seed):
    rng = np.random.default_rng(seed)
    s = K * q ** (-2 / 3); lam = K ** 3 / 6
    cnt_tree = np.zeros(12); cnt_all = np.zeros(60); nontree = 0
    nb_bins = 4; mark_hist = np.zeros(nb_bins)
    nbs = rng.binomial(q, s, size=N); ncs = rng.binomial(q, s, size=N)
    for i in range(N):
        nb, nc = nbs[i], ncs[i]
        if nb == 0 or nc == 0:
            cnt_tree[0] += 1; cnt_all[0] += 1; continue
        a = rng.random(nb) * s; c = rng.random(nc) * s
        d = a[:, None] + c[None, :] + rng.random((nb, nc))
        t = d <= s
        n = int(t.sum()); cnt_all[min(n, 59)] += 1
        if n == 0:
            cnt_tree[0] += 1; continue
        if t.sum(axis=1).max() > 1 or t.sum(axis=0).max() > 1:
            nontree += 1; continue
        if n < 12: cnt_tree[n] += 1
        if n == 1:
            xi = d[t][0] / s * K
            mark_hist[min(int(xi / K * nb_bins), nb_bins - 1)] += 1
    # epsilon_1 at r=1
    Q1 = 1 + 2 * (2 * lam) + (2 * lam) ** 2 + 4 * lam   # E(1+Z1)^2, Z1 = 2 Poisson(lam)
    eps1 = 0.9 * K * K * q ** (-1 / 3) * lam + (2 * lam + lam * lam) * K * q ** (-2 / 3) + Q1 / q + lam * lam / q ** 2
    EX = q * q * (q - 1) * s ** 5 / 20
    print(f"q={q} K={K} N={N} lam={lam:.4f}  P(nontree)={nontree/N:.4f}  E X_o={EX:.4f}  eps1={eps1:.4f}")
    tv = 0.0
    for n in range(6):
        nu = math.exp(-lam) * lam ** n / math.factorial(n)
        mu = cnt_tree[n] / N; se = math.sqrt(max(mu * (1 - mu), 1e-12) / N)
        bound = 2 * n * (2 * n + 1) / q + lam * lam / q ** 2 + n * (n - 1) * K * q ** (-2 / 3) + 0.9 * n * K * K * q ** (-1 / 3) * (1 - 1 / q)
        f = mu / nu
        flag = "VIOLATION" if f + 3 * se / nu < 1 - bound else ""
        tv += max(nu - mu, 0)
        print(f"   n={n}: mu(tree,n)={mu:.5f}+-{se:.5f}  nu(n)={nu:.5f}  f={f:.4f}  claimed lower bound on f={1-bound:.4f} {flag}")
    print(f"   TV on the partition {{nontree, tree with n teams}} (lower estimate of d_TV) = {tv:.4f}  vs eps1 = {eps1:.4f} {'VIOLATION' if tv > eps1 + 0.01 else ''}")
    # degree law (all balls) vs Poisson: mean, var
    ks = np.arange(60); tot = cnt_all.sum(); mean = (ks * cnt_all).sum() / tot; var = (ks * ks * cnt_all).sum() / tot - mean ** 2
    print(f"   degree of o: mean={mean:.4f} (lam={lam:.4f}), var={var:.4f} (lam + 2 E X_o = {lam + 2*EX:.4f}), P(deg=0)={cnt_all[0]/tot:.5f} vs e^-lam={math.exp(-lam):.5f}")
    # marks of the single team in 1-team trees, by bin: f(bin) vs lower bound with E[(s-x_p)^2/2 | d] = (s^2 - 2 s d/3 + d^2/6)/2
    edges = np.linspace(0, K, nb_bins + 1)
    for j in range(nb_bins):
        lo, hi = edges[j], edges[j + 1]
        nu = math.exp(-lam) * (hi ** 3 - lo ** 3) / 6
        # average of the conditional F1 term over the bin under density prop. to xi^2
        xs = np.linspace(lo, hi, 2001); w = xs ** 2; dd = xs * q ** (-2 / 3)
        term = (s * s - 2 * s * dd / 3 + dd * dd / 6) / 2
        avg = (w * term).sum() / w.sum()
        bound = 6 / q + lam * lam / q ** 2 + (q - 1) * 3 * avg
        mu = mark_hist[j] / N; se = math.sqrt(max(mu * (1 - mu), 1e-12) / N)
        print(f"   1-team trees, mark in [{lo:.2f},{hi:.2f}]: f={mu/nu:.4f}+-{se/nu:.4f}  claimed lower bound={1-bound:.4f} {'VIOLATION' if mu/nu + 3*se/nu < 1-bound else ''}")
    sys.stdout.flush()

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    for K in (1.0, 2.0):
        for q in (64, 128, 256, 4096, 262144):
            run(q, K, N, 1234 + q)
