"""sc (I-b) numerics: near-tight hypergraph H (xi<=K) in the pair-additive market (PA) versus the
independent-teams model (IT, 'Poissonised').  Compares the joint law of (|H|, s(H)), s = 3|H| - #pairs
used (pair-sharing excess), the void probability against the Janson bound exp(-mu+Delta), and the
star-size law against Poisson(K^3/6).  Evidence only."""
import sys, math, numpy as np
from collections import Counter
from scipy.stats import poisson

def excess(a, b, c, q):
    n = len(a)
    if n == 0: return 0, 0
    y = 0; used = 0
    for key in (a * q + b, b * q + c, a * q + c):
        u, m = np.unique(key, return_counts=True); used += len(u); y += int((m * (m - 1) // 2).sum())
    return 3 * n - used, y

def run(q, Ks, trials, seed):
    rng = np.random.default_rng(seed); out = {}
    for K in Ks: out[K] = dict(pa=Counter(), it=Counter(), Ypa=0.0, degpa=Counter(), degit=Counter(), tvpa=0.0, tvit=0.0, cov=[0.0, 0.0, 0.0])
    for _ in range(trials):
        A = rng.random((q, q)); B = rng.random((q, q)); C = rng.random((q, q))
        cost = A[:, :, None] + B[None, :, :] + C[:, None, :]
        for K in Ks:
            t = K * q ** (-2 / 3); lam = K ** 3 / 6; o = out[K]
            a, b, c = np.nonzero(cost <= t); s, y = excess(a, b, c, q)
            o["pa"][(len(a), min(s, 3))] += 1; o["Ypa"] += y
            deg = np.bincount(a, minlength=q); h = np.bincount(deg, minlength=40)[:40] / q
            o["tvpa"] += 0.5 * np.abs(h - poisson.pmf(np.arange(40), lam)).sum()
            for d in deg: o["degpa"][int(d)] += 1
            z0 = (deg[0] == 0); z1 = (deg[1] == 0); o["cov"][0] += z0; o["cov"][1] += z1; o["cov"][2] += z0 and z1
            n = rng.binomial(q ** 3, t ** 3 / 6); idx = rng.choice(q ** 3, size=n, replace=False)
            a, b, c = idx // (q * q), (idx // q) % q, idx % q; s, y = excess(a, b, c, q)
            o["it"][(n, min(s, 3))] += 1
            deg = np.bincount(a, minlength=q); h = np.bincount(deg, minlength=40)[:40] / q
            o["tvit"] += 0.5 * np.abs(h - poisson.pmf(np.arange(40), lam)).sum()
    for K in Ks:
        o = out[K]; t = K * q ** (-2 / 3); mu = q ** 3 * t ** 3 / 6; Delta = 3 * q ** 3 * (q - 1) / 2 * t ** 5 / 20
        print("q=%d K=%.2f trials=%d  mu=%.3f Delta=%.3f  E[Y]_PA=%.3f (Y = sharing pairs of teams; theory Delta)" % (q, K, trials, mu, Delta, o["Ypa"] / trials))
        p0 = sum(v for (n, s), v in o["pa"].items() if n == 0) / trials
        print("  void: P_PA(|H|=0)=%.5f  IT exact (1-p)^(q^3)=%.5f  Janson bound exp(-mu+Delta)=%.5f" % (p0, (1 - t ** 3 / 6) ** (q ** 3), math.exp(-mu + Delta)))
        for s in (0, 1, 2):
            pa = sum(v for (n, ss), v in o["pa"].items() if ss == s) / trials; it = sum(v for (n, ss), v in o["it"].items() if ss == s) / trials
            print("  P(s=%d): PA=%.4f IT=%.4f ratio=%.2f" % (s, pa, it, pa / it if it else float("nan")))
        ns = sorted(set(n for (n, s) in o["pa"]) | set(n for (n, s) in o["it"]))
        worst = {0: (0, None), 1: (0, None)}
        for s in (0, 1):
            for n in ns:
                x, y = o["pa"][(n, s)], o["it"][(n, s)]
                if x >= 30 and y >= 30:
                    r = x / y
                    if r > worst[s][0]: worst[s] = (r, n)
        b0 = math.exp(Delta + K ** 3 * (ns[len(ns) // 2]) / (2 * q) + K ** 6 / (36 * q))
        print("  max over n (cells with >=30 counts) of P_PA(|H|=n,s)/P_IT(|H|=n,s): s=0: %.3f at n=%s (bound exp(Delta+K^3 n/2q) ~ %.2f);  s=1: %.2f at n=%s (one-sharing factor 1.8 q^(2/3)/K = %.2f)" % (worst[0][0], worst[0][1], b0, worst[1][0], worst[1][1], 1.8 * q ** (2 / 3) / K))
        tot = sum(o["degpa"].values()); lam = K ** 3 / 6
        tvmean = 0.5 * sum(abs(o["degpa"][d] / tot - poisson.pmf(d, lam)) for d in range(40))
        c0, c1, c01 = [x / trials for x in o["cov"]]
        print("  star size: TV(mean PA law, Poisson(%.3f))=%.4f ; mean per-instance TV(empirical, Poisson): PA=%.4f IT=%.4f ; Cov(1{deg_a=0},1{deg_a'=0})=%.5f" % (lam, tvmean, o["tvpa"] / trials, o["tvit"] / trials, c01 - c0 * c1))
    sys.stdout.flush()

if __name__ == "__main__":
    for q, trials in ((16, 40000), (24, 40000), (32, 30000), (48, 15000)):
        run(q, (0.8, 1.0, 1.3, 2.0), trials, 11 + q)
