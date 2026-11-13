"""V2 referee, front 2: brute-force checks of the finite configuration-model formula
E#{PM with bin profile m} = prod_s N_s(m) * prod_b C(n_b, m_b)^(-2)
(N_s(m) = number of choice functions in sector s with profile m)."""
import itertools, math, random, sys
from fractions import Fraction
import numpy as np
from scipy.optimize import minimize

def choice_counts(deg):
    B = len(deg[0]); poly = {tuple([0] * B): 1}
    for d in deg:
        new = {}
        for m, c in poly.items():
            for b in range(B):
                if d[b]:
                    mm = m[:b] + (m[b] + 1,) + m[b + 1:]
                    new[mm] = new.get(mm, 0) + c * d[b]
        poly = new
    return poly

def formula(degs):
    B = len(degs[0][0]); nb = [sum(d[b] for d in degs[0]) for b in range(B)]
    N = [choice_counts(degs[s]) for s in range(3)]; out = {}
    for m in N[0]:
        if m in N[1] and m in N[2]:
            v = Fraction(N[0][m] * N[1][m] * N[2][m])
            for b in range(B): v /= math.comb(nb[b], m[b]) ** 2
            out[m] = v
    return out

def owners(degs, s, b):
    return [v for v, d in enumerate(degs[s]) for _ in range(d[b])]

def count_pms(edges_at, n, B):
    """edges_at[v1] = list of (v2, v3, bin); returns dict profile -> number of perfect matchings."""
    res = {}
    def rec(v, used2, used3, prof):
        if v == n:
            res[prof] = res.get(prof, 0) + 1; return
        for (a, c, b) in edges_at[v]:
            if not (used2 >> a) & 1 and not (used3 >> c) & 1:
                rec(v + 1, used2 | (1 << a), used3 | (1 << c), prof[:b] + (prof[b] + 1,) + prof[b + 1:])
    rec(0, 0, 0, tuple([0] * B)); return res

def build(degs, perms):
    n = len(degs[0]); B = len(degs[0][0]); edges_at = [[] for _ in range(n)]
    for b in range(B):
        o1, o2, o3 = owners(degs, 0, b), owners(degs, 1, b), owners(degs, 2, b)
        p2, p3 = perms[b]
        for i in range(len(o1)): edges_at[o1[i]].append((o2[p2[i]], o3[p3[i]], b))
    return edges_at

def exact(degs):
    n = len(degs[0]); B = len(degs[0][0]); nb = [sum(d[b] for d in degs[0]) for b in range(B)]
    tot = {}; cnt = 0
    spaces = [list(itertools.product(itertools.permutations(range(nb[b])), repeat=2)) for b in range(B)]
    for perms in itertools.product(*spaces):
        cnt += 1
        for m, c in count_pms(build(degs, perms), n, B).items(): tot[m] = tot.get(m, 0) + c
    return {m: Fraction(c, cnt) for m, c in tot.items()}, cnt

def mc(degs, S, rng):
    n = len(degs[0]); B = len(degs[0][0]); nb = [sum(d[b] for d in degs[0]) for b in range(B)]
    acc = {}
    for _ in range(S):
        perms = [(rng.permutation(nb[b]), rng.permutation(nb[b])) for b in range(B)]
        for m, c in count_pms(build(degs, perms), n, B).items(): acc.setdefault(m, []).append(c)
    return {m: (sum(v) / S, math.sqrt(max(sum(x * x for x in v) / S - (sum(v) / S) ** 2, 0) / S)) for m, v in acc.items()}

def rand_degs(n, B, rng, mean=0.9):
    while True:
        d0 = rng.poisson(mean, size=(n, B))
        if (d0.sum(1) == 0).any() or (d0.sum(0) == 0).any(): continue
        degs = [[tuple(int(x) for x in r) for r in d0]]
        ok = True
        for s in (1, 2):
            d = np.zeros((n, B), dtype=int)
            for b in range(B):
                for v in rng.integers(0, n, size=d0[:, b].sum()): d[v, b] += 1
            if (d.sum(1) == 0).any(): ok = False; break
            degs.append([tuple(int(x) for x in r) for r in d])
        if ok: return degs

if __name__ == "__main__":
    rng = np.random.default_rng(7)
    print("== A. exact enumeration over ALL pairings (n=3 per part, 2 bins)")
    degs = [[(2, 1), (1, 1), (1, 1)], [(1, 1), (2, 1), (1, 1)], [(1, 1), (1, 1), (2, 1)]]
    ex, cnt = exact(degs); fo = formula(degs)
    keys = sorted(set(ex) | set(k for k, v in fo.items() if v))
    print(" pairings enumerated:", cnt)
    for m in keys: print("  profile", m, "exact", ex.get(m, 0), "formula", fo.get(m, 0), "EQUAL" if ex.get(m, 0) == fo.get(m, 0) else "DIFFER")
    print("== B. Monte Carlo, random degree sequences")
    for (n, B, S) in ((6, 2, 20000), (7, 3, 20000), (8, 3, 12000), (10, 2, 6000)):
        degs = rand_degs(n, B, rng); fo = formula(degs); est = mc(degs, S, rng)
        worst = 0; tf = sum(float(v) for v in fo.values()); te = sum(v[0] for v in est.values())
        for m, v in fo.items():
            if v == 0: continue
            mu, se = est.get(m, (0.0, 0.0))
            if se > 0: worst = max(worst, abs(mu - float(v)) / se)
        print("  n=%d bins=%d samples=%d profiles=%d  total formula=%.5f  MC=%.5f  worst |z| over profiles=%.2f" % (n, B, S, len(fo), tf, te, worst))
    print("== C. d-regular single bin")
    for (n, d, S) in ((6, 3, 20000), (8, 3, 8000), (6, 4, 8000)):
        degs = [[(d,)] * n] * 3; fo = float(formula(degs)[(n,)]); est = mc(degs, S, rng)[(n,)]
        rate = math.log(d) - 2 * (d - 1) * math.log(d / (d - 1))
        print("  n=%d d=%d formula=%.5f MC=%.5f +- %.5f | (1/n)ln formula=%.4f, limit rate ln d-2(d-1)ln(d/(d-1))=%.4f" % (n, d, fo, est[0], est[1], math.log(fo) / n, rate))
    for n in (100, 1000, 10000):
        for d in (6, 7):
            lf = 3 * n * math.log(d) - 2 * (math.lgamma(d * n + 1) - math.lgamma(n + 1) - math.lgamma((d - 1) * n + 1))
            print("  n=%d d=%d (1/n)ln E = %.5f  vs rate %.5f" % (n, d, lf / n, math.log(d) - 2 * (d - 1) * math.log(d / (d - 1))))
    print("== D. entropy rate: (1/q) ln N(m) vs inf_z [mean ln W(z) - sum p ln z], Poisson stars, 3 bins")
    lam = np.array([0.3, 0.8, 1.5])
    for q in (60, 120, 240):
        deg = []
        while len(deg) < q:
            d = rng.poisson(lam)
            if d.sum() > 0: deg.append(tuple(int(x) for x in d))
        N = choice_counts(deg); D = np.array(deg, dtype=float)
        for frac in ((0.25, 0.4, 0.35), (0.15, 0.35, 0.5), (0.05, 0.3, 0.65)):
            m = [round(frac[0] * q), round(frac[1] * q)]; m.append(q - sum(m)); m = tuple(m)
            if m not in N: print("  q=%d m=%s infeasible" % (q, m)); continue
            p = np.array(m) / q
            f = lambda u: np.log(D @ np.exp(u)).mean() - p @ u
            r = minimize(f, np.zeros(3), method="BFGS")
            print("  q=%d m=%s exact (1/q)lnN=%.5f  dual inf=%.5f  gap=%.5f  (ln q)/q=%.5f" % (q, m, math.log(N[m]) / q, r.fun, r.fun - math.log(N[m]) / q, math.log(q) / q))
