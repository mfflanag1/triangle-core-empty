"""V12 finite-q checks.
 P1  Lemma 2.1 by enumeration of ALL configurations (non-simple ones included), q=2, B=2, random degree arrays.
 P2  Lemma 2.2 direction by enumeration of ALL labelled low hypergraphs: q=2,B=2 (3^8) and q=3,B=1 (<=5 low teams):
     sum_{g in G(D)} Y_m(g) <= (prod n_b!)^2/prod D! * prod_s N_s * prod_b C(n_b,m_b)^-2 for every D and m;
     and sum_D Q(D) <= exp(q^3 p^2/(1-p)).
 P3  chain at q=3 (K=2,B=2): Monte Carlo E[Z], P(Z>=1)  vs  exact tilt  vs  the split bound (2.1) averaged over G.
 P4  P(min assignment deficit < 2.302 q^(1/3)) for q=3,4,5 by simulation."""
import itertools, math, sys, time
from fractions import Fraction
from collections import defaultdict
import numpy as np
rng = np.random.default_rng(2026)

def Ns(m, Ds):
    """number of choice functions: coefficient of z^m z_inf^(q-|m|) in prod_v (z_inf + sum_b D_vb z_b)."""
    poly = {tuple([0] * len(m)): 1}
    for row in Ds:
        new = defaultdict(int)
        for k, c in poly.items():
            new[k] += c
            for b, dvb in enumerate(row):
                if dvb and k[b] < m[b]:
                    kk = list(k); kk[b] += 1; new[tuple(kk)] += c * dvb
        poly = new
    return poly.get(tuple(m), 0)

def profiles(q, B):
    return [m for m in itertools.product(range(q + 1), repeat=B) if sum(m) <= q]

# ---------- P1 ----------
def p1(q=2, B=2, ntry=40):
    bad = 0; cnt = 0
    for _ in range(ntry):
        n = [int(rng.integers(1, 4)), int(rng.integers(0, 3))]
        # random arrays with these column sums
        D = [[[0] * B for _ in range(q)] for _ in range(3)]
        for s in range(3):
            for b in range(B):
                for _ in range(n[b]): D[s][int(rng.integers(q))][b] += 1
        he = [[[v for v in range(q) for _ in range(D[s][v][b])] for b in range(B)] for s in range(3)]   # half-edges -> agent
        tot = defaultdict(int)
        perms = [list(itertools.permutations(range(n[b]))) for b in range(B)]
        for cfg in itertools.product(*[itertools.product(perms[b], perms[b]) for b in range(B)]):
            triples = [(b, he[0][b][i], he[1][b][cfg[b][0][i]], he[2][b][cfg[b][1][i]]) for b in range(B) for i in range(n[b])]
            for r in range(0, q + 1):
                for S in itertools.combinations(range(len(triples)), r):
                    ok = all(len({triples[i][c] for i in S}) == r for c in (1, 2, 3))
                    if ok:
                        m = [0] * B
                        for i in S: m[triples[i][0]] += 1
                        tot[tuple(m)] += 1
        for m in profiles(q, B):
            rhs = 1
            for s in range(3): rhs *= Ns(m, D[s])
            for b in range(B):
                rhs *= (math.factorial(m[b]) * math.factorial(n[b] - m[b])) ** 2 if m[b] <= n[b] else 0
            cnt += 1; bad += (rhs != tot.get(m, 0))
    print("P1 Lemma 2.1 (q=%d,B=%d): %d mismatches in %d (array, profile) pairs" % (q, B, bad, cnt), flush=True)

# ---------- P2 ----------
def p2(q, B, maxlow=None, pb=None):
    teams = list(itertools.product(range(q), repeat=3)); T = len(teams)
    lhs = defaultdict(lambda: defaultdict(int)); ngraph = defaultdict(int)
    def graphs():
        if maxlow is None:
            for lab in itertools.product(range(B + 1), repeat=T):
                yield [(teams[i], l - 1) for i, l in enumerate(lab) if l]
        else:
            for k in range(maxlow + 1):
                for S in itertools.combinations(range(T), k):
                    for lab in itertools.product(range(B), repeat=k):
                        yield [(teams[i], l) for i, l in zip(S, lab)]
    for g in graphs():
        D = [[[0] * B for _ in range(q)] for _ in range(3)]
        for e, b in g:
            for s in range(3): D[s][e[s]][b] += 1
        key = tuple(tuple(tuple(r) for r in Ds) for Ds in D)
        ngraph[key] += 1
        for r in range(0, q + 1):
            for S in itertools.combinations(g, r):
                if all(len({e[c] for e, _ in S}) == r for c in range(3)):
                    m = [0] * B
                    for _, b in S: m[b] += 1
                    lhs[key][tuple(m)] += 1
    viol = 0; cnt = 0; minratio = None; tight = 0
    for key in ngraph:
        n = [sum(key[0][v][b] for v in range(q)) for b in range(B)]
        pref = Fraction(math.prod(math.factorial(x) ** 2 for x in n), math.prod(math.factorial(x) for Ds in key for r in Ds for x in r))
        assert ngraph[key] <= pref        # number of simple graphs <= number of configurations / prod D!
        for m in profiles(q, B):
            L = lhs[key].get(m, 0)
            if any(m[b] > n[b] for b in range(B)):
                assert L == 0; continue
            rhs = pref * math.prod(Ns(m, key[s]) for s in range(3)) / math.prod(math.comb(n[b], m[b]) ** 2 for b in range(B))
            cnt += 1; viol += (L > rhs); tight += (L == rhs)
            if L > 0: minratio = min(minratio, rhs / L) if minratio is not None else rhs / L
    print("P2 Lemma 2.2 (q=%d,B=%d,maxlow=%s): degree arrays %d, (D,m) pairs %d, violations %d, equalities %d, min RHS/LHS = %.4f"
          % (q, B, maxlow, len(ngraph), cnt, viol, tight, float(minratio)), flush=True)
    if pb is not None and maxlow is None:
        p = sum(pb); totQ = Fraction(0); totP = Fraction(0)
        for key in ngraph:
            n = [sum(key[0][v][b] for v in range(q)) for b in range(B)]
            base = math.prod(Fraction(pb[b]) ** n[b] for b in range(B)) * (1 - p) ** (T - sum(n))
            totQ += base * Fraction(math.prod(math.factorial(x) ** 2 for x in n), math.prod(math.factorial(x) for Ds in key for r in Ds for x in r))
            totP += base * ngraph[key]
        print("   with p_b=%s: sum_D P(D)=%s (must be 1), sum over REALISED arrays of Q(D) = %.6f <= exp(q^3 p^2/(1-p)) = %.6f"
              % ([str(x) for x in pb], totP, float(totQ), math.exp(T * float(p) ** 2 / (1 - float(p)))), flush=True)

# ---------- P3 ----------
def p3(q=3, K=2.0, B=2, nsim=400000):
    eta = K / B; s = q ** (2 / 3); assert s >= K
    lam = [eta ** 3 * (3 * b * b - 3 * b + 1) / 6 for b in range(1, B + 1)]; p_b = [l / q ** 2 for l in lam]; p = sum(p_b)
    perms = list(itertools.permutations(range(q)))
    idx = np.array([[i * q * q + a[i] * q + c[i] for i in range(q)] for a in perms for c in perms])      # (q!^2, q)
    xi = s * rng.random((nsim, q ** 3, 3)).sum(axis=2)
    X = xi[:, idx].sum(axis=2)                                                                            # (nsim, q!^2)
    for beta in (0.8, 1.2, 1.6, 2.0, 2.302):
        Z = (X < beta * q).sum(axis=1); EZ = Z.mean(); PZ = (Z > 0).mean(); se = Z.std() / math.sqrt(nsim)
        best = None
        for th in np.linspace(0.2, 6, 117):
            A = lambda x: -math.exp(-th * x) * (x * x / (2 * th) + x / th ** 2 + 1 / th ** 3)
            kap = [A(eta * b) - A(eta * (b - 1)) for b in range(1, B + 1)]; J = -A(K)
            a = th * s; exact = math.exp(th * beta * q) * len(idx) * ((1 - math.exp(-a)) / a) ** (3 * q)
            split = 0.0
            for m in profiles(q, B):
                r = q - sum(m)
                EY = math.comb(q, r) ** 3 * math.factorial(q - r) ** 3 / math.prod(math.factorial(x) for x in m) * math.prod(pp ** x for pp, x in zip(p_b, m))
                split += EY * math.prod((k / l) ** x for k, l, x in zip(kap, lam, m)) * math.factorial(r) ** 2 * (J / (q * q * (1 - p))) ** r
            split *= math.exp(th * beta * q)
            assert split >= exact * (1 - 1e-12), (th, split, exact)
            if best is None or split < best[1]: best = (th, split, exact)
        print("P3 q=%d beta=%.3f: MC P(Z>=1)=%.5f  MC E[Z]=%.5f (se %.5f)  <=  exact tilt %.5f  <=  split bound (2.1) %.5f  (theta=%.2f; split>=exact asserted on the theta grid)  direction ok: %s"
              % (q, beta, PZ, EZ, se, best[2], best[1], best[0], EZ <= best[2] <= best[1]), flush=True)

# ---------- P4 ----------
def p4():
    for q, nsim in ((3, 200000), (4, 100000), (5, 4000)):
        perms = list(itertools.permutations(range(q)))
        idx = np.array([[i * q * q + a[i] * q + c[i] for i in range(q)] for a in perms for c in perms])
        mins = []
        for _ in range(0, nsim, 2000):
            d = rng.random((2000, q ** 3, 3)).sum(axis=2)
            mins.append(d[:, idx].sum(axis=2).min(axis=1))
        mins = np.concatenate(mins); thr = 2.302 * q ** (1 / 3)
        print("P4 q=%d: threshold 2.302 q^(1/3) = %.4f; P(min deficit < threshold) = %.4f (n=%d); E[min]/q^(1/3) = %.4f"
              % (q, thr, (mins < thr).mean(), len(mins), mins.mean() / q ** (1 / 3)), flush=True)

if __name__ == "__main__":
    t = time.time()
    p1(); p2(2, 2, None, [Fraction(1, 24), Fraction(7, 24)]); p2(2, 1, None, [Fraction(1, 3)]); p2(3, 1, 5)
    p3(); p4()
    print("total %.0fs" % (time.time() - t))
