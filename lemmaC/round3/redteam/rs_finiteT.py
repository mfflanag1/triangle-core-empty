"""Finite-temperature replica-symmetric (Bethe) entropy for 3-index matching on the Poisson hypertree, intensity t^r/r!.
Cavity: exp(-b X) = sum_k exp(-b (xi_k - X'_k - X''_k)).  Per sector-1 agent:
 f = 3 E[X] - 2 E int t^r/r! * (-(1/b) log(1+exp(-b (t - S3)))) dt,   e = E int t^(r+1)/r! / (1+exp(b (t - S3))) dt,   s = b (e - f).
Question: does s(b) cross zero at finite b (entropy crisis, as in Martin-Mezard-Rivoire for uniform costs)? Monte Carlo; evidence only."""
import numpy as np, sys
from math import factorial
rng = np.random.default_rng(20261111)
def sweep(pop, b, r, Lam, chunk=4000):
    N = len(pop); new = np.empty(N); mean = Lam ** (r + 1) / factorial(r + 1)
    for a in range(0, N, chunk):
        n = min(chunk, N - a)
        K = rng.poisson(mean, n); K = np.maximum(K, 1)
        tot = int(K.sum()); st = np.concatenate([[0], np.cumsum(K)[:-1]])
        xi = Lam * rng.random(tot) ** (1 / (r + 1))
        arg = -b * (xi - pop[rng.integers(0, N, tot)] - pop[rng.integers(0, N, tot)])
        mx = np.maximum.reduceat(arg, st)
        ssum = np.add.reduceat(np.exp(arg - np.repeat(mx, K)), st)
        new[a:a + n] = -(mx + np.log(ssum)) / b
    return new
def observables(pop, b, r, nS=20000):
    N = len(pop)
    S = pop[rng.integers(0, N, nS)] + pop[rng.integers(0, N, nS)] + pop[rng.integers(0, N, nS)]
    t = np.linspace(0, S.max() + 40 / b + 2, 1500); dt = t[1] - t[0]
    fe = en = nm = 0.0
    for a in range(0, nS, 2000):
        u = b * (t[None, :] - S[a:a + 2000, None])
        h = -np.logaddexp(0, -u) / b
        occ = 0.5 * (1 - np.tanh(u / 2))
        w = t[None, :] ** r / factorial(r)
        fe += (w * h).sum() * dt; en += (w * t[None, :] * occ).sum() * dt; nm += (w * occ).sum() * dt
    return fe / nS, en / nS, nm / nS
def run(b, r=2, N=12000, sweeps=70, avg=20):
    Lam = max(7.0, 5.0 + 9.0 / b)
    pop = rng.normal(0.25, 0.9, N); acc = []
    for it in range(sweeps):
        new = sweep(pop, b, r, Lam); m = rng.random(N) < 0.3; pop = np.where(m, new, pop)   # partial replacement: the uniform-shift mode has multiplier -2 under full sweeps
        if it >= sweeps - avg:
            fe, en, nm = observables(pop, b, r)
            acc.append((3 * pop.mean() - 2 * fe, en, nm, pop.mean(), pop.std()))
    a = np.array(acc); m = a.mean(0); se = a.std(0, ddof=1) / len(a) ** .5
    f, e = m[0], m[1]
    return dict(beta=b, f=round(f, 4), e=round(e, 4), s=round(b * (e - f), 4), s_se=round(b * (se[0] ** 2 + se[1] ** 2) ** .5, 4),
                norm=round(m[2], 4), meanX=round(m[3], 4), sdX=round(m[4], 4))
if __name__ == '__main__':
    r = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    for b in [float(x) for x in sys.argv[2:]] or [1, 2, 3, 4, 6, 8, 12, 20, 40]:
        print(run(b, r=r), flush=True)
