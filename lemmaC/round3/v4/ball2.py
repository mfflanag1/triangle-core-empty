"""Referee check of Theorems A (r = 2) and B (two roots, r = 1). Floating point, evidence only.

Part 1: 2-ball of a sector-A root in a fresh market; unmarked shape = multiset of (child-team counts of
the B-child, of the C-child) over the root's teams. Hypertree iff #agents = 1 + 2 #teams.
Compare with the exact Poisson-tree shape law and with eps_1(q, K, 2).
Part 2: two roots o, o' (same sector, and different sectors), r = 1: joint law of (tree?, degree) against
the product of Poissons, compared with 4 eps_1(q, K, 1); and Var of N = #agents with a given 1-ball type.
"""
import sys, math
from collections import Counter
import numpy as np

def pois(lam, n): return math.exp(-lam) * lam ** n / math.factorial(n)

def eps1(q, K, r):
    lam = K ** 3 / 6; S = sum((2 * lam) ** j for j in range(r))
    # exact Q_r = E V^2 for V = sum_{j<=r} Z_j, GW with offspring 2*Poisson(lam): by simulation-free recursion
    # E[V_r] and E[V_r^2] via V_r = 1 + sum_{i<=X} V_{r-1}^{(i)}, X = 2 Poisson(lam)
    m1, m2 = 1.0, 1.0
    mu, var = 2 * lam, 4 * lam
    for _ in range(r):
        s1 = mu * m1; s2 = mu * (m2 - m1 * m1) + (var + mu * mu) * m1 * m1
        m1, m2 = 1 + s1, 1 + 2 * s1 + s2
    return 0.9 * K * K * q ** (-1 / 3) * lam * S + (2 * lam + lam * lam) * S * K * q ** (-2 / 3) + m2 / q + lam * lam * S / q ** 2

def part1(q, K, N, seed):
    rng = np.random.default_rng(seed); s = K * q ** (-2 / 3); lam = K ** 3 / 6
    shapes = Counter(); nontree = 0
    for _ in range(N):
        A = rng.random((q, q)); B = rng.random((q, q)); C = rng.random((q, q))  # A[a,b], B[b,c], C[a,c]
        o = 0
        t0 = np.argwhere(A[o][:, None] + B + C[o][None, :] <= s)      # (b,c)
        teams = set((o, int(b), int(c)) for b, c in t0)
        bs = set(int(b) for b, c in t0); cs = set(int(c) for b, c in t0)
        kids = {}
        for b in bs:
            tb = np.argwhere(A[:, b][:, None] + B[b][None, :] + C <= s)   # (a,c)
            kids[('b', b)] = len(tb) - 1
            teams.update((int(a), b, int(c)) for a, c in tb)
        for c in cs:
            tc = np.argwhere(A + B[:, c][None, :] + C[:, c][:, None] <= s)  # (a,b)
            kids[('c', c)] = len(tc) - 1
            teams.update((int(a), int(b), c) for a, b in tc)
        ag = set()
        for a, b, c in teams: ag.add(('a', a)); ag.add(('b', b)); ag.add(('c', c))
        ag.add(('a', o))
        if len(ag) != 1 + 2 * len(teams):
            nontree += 1; continue
        shapes[tuple(sorted((kids[('b', int(b))], kids[('c', int(c))]) for b, c in t0))] += 1
    e1 = eps1(q, K, 2)
    print(f"[r=2] q={q} K={K} N={N}: P(nontree)={nontree/N:.4f}  eps1={e1:.4f}")
    tv = 0.0; rows = []
    for sh, cnt in shapes.items():
        n0 = len(sh); mult = math.factorial(n0)
        for v in Counter(sh).values(): mult //= math.factorial(v)
        nu = pois(lam, n0) * mult * math.prod(pois(lam, x) * pois(lam, y) for x, y in sh)
        rows.append((nu, sh, cnt / N))
    rows.sort(reverse=True)
    seen = 0.0
    for nu, sh, mu in rows:
        seen += nu
        if nu * N >= 50:   # shapes with enough expected count to compare
            tv += max(nu - mu, 0)
    for nu, sh, mu in rows[:6]:
        se = math.sqrt(mu * (1 - mu) / N)
        E = len(sh) + sum(x + y for x, y in sh); V = 1 + 2 * E; m = 1 + 2 * len(sh)
        P2 = len(sh) * (len(sh) - 1) // 2 + sum(x * (x + 1) // 2 + y * (y + 1) // 2 for x, y in sh)
        bound = 2 * E * V / q + m * lam * lam / q ** 2 + 2 * P2 * K * q ** (-2 / 3) + 0.9 * E * K * K * q ** (-1 / 3) * (1 - 1 / q)
        print(f"   shape {sh}: mu={mu:.5f}+-{se:.5f} nu={nu:.5f} f={mu/nu:.4f} claimed lower bound={1-bound:.4f} {'VIOLATION' if mu/nu + 3*se/nu < 1-bound else ''}")
    print(f"   sum over well-sampled shapes of (nu-mu)_+ = {tv:.4f} (lower estimate of d_TV; nu-mass of unseen shapes {1-seen:.4f}) vs eps1={e1:.4f} {'VIOLATION' if tv > e1 + 0.01 else ''}")
    sys.stdout.flush()

def ball1(aval, cval, Bsub, s):
    """1-ball of a root given its two deficit rows and the matrix of the opposite pairs. returns (deg, tree?)"""
    t = aval[:, None] + cval[None, :] + Bsub <= s
    n = int(t.sum())
    tree = n == 0 or (t.sum(axis=1).max() <= 1 and t.sum(axis=0).max() <= 1)
    return n, tree

def part2(q, K, N, seed, same_sector):
    rng = np.random.default_rng(seed); s = K * q ** (-2 / 3); lam = K ** 3 / 6
    joint = Counter(); 
    for _ in range(N):
        A = rng.random((q, q)); B = rng.random((q, q)); C = rng.random((q, q))
        n1, t1 = ball1(A[0], C[0], B, s)                       # root a=0: teams (0,b,c)
        if same_sector:
            n2, t2 = ball1(A[1], C[1], B, s)                   # root a=1
        else:
            n2, t2 = ball1(A[:, 0], B[0], C, s)                # root b=0: teams (a,0,c): A[a,0]+B[0,c]+C[a,c]
        joint[(n1 if t1 else -1, n2 if t2 else -1)] += 1
    e1 = eps1(q, K, 1)
    tv = 0.0
    for i in range(0, 6):
        for j in range(0, 6):
            tv += max(pois(lam, i) * pois(lam, j) - joint[(i, j)] / N, 0)
    # covariance of indicators h = 1{1-ball is a tree with exactly one team}, and 1{deg = 0}
    out = []
    for typ in (0, 1):
        p1 = sum(v for (i, j), v in joint.items() if i == typ) / N
        p2 = sum(v for (i, j), v in joint.items() if j == typ) / N
        p12 = joint[(typ, typ)] / N
        out.append((typ, p12 - p1 * p2, math.sqrt(p12 * (1 - p12) / N)))
    print(f"[two roots, r=1, {'same sector' if same_sector else 'different sectors'}] q={q} K={K} N={N}: TV lower estimate on (tree?,deg) partition = {tv:.4f} vs 4*eps1 = {4*e1:.4f} {'VIOLATION' if tv > 4*e1 + 0.01 else ''}")
    for typ, cov, se in out:
        print(f"   type 'tree with {typ} teams': Cov = {cov:+.5f} (se ~{se:.5f}); claimed |Cov| <= 6 eps1 = {6*e1:.4f}")
    sys.stdout.flush()

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    for K in (1.0, 1.26):
        for q in (64, 128):
            part1(q, K, N, 77 + q)
    for K in (1.0, 2.0):
        for q in (64, 128):
            part2(q, K, N, 99 + q, True); part2(q, K, N, 199 + q, False)
