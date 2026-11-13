"""V6 brute force for the finite-q inequalities of THEOREM_I.md section 2 (exact integers).
For ALL labelled low hypergraphs g on q+q+q agents with at most NMAX low teams and B bins:
 (L2.2)  sum_{g in G(D)} Y_m(g) * prod D!  <=  prod_s N_s(m;D^s) * prod_b (m_b! (n_b-m_b)!)^2     for every D, m
 (2.1c)  A_m(g) := #{perfect assignments N with profile(N cap g) = m}  <=  Y_m(g) * (r!)^2        for every g, m
 total-mass identity: sum over arrays D with column sums n of |G(D)| = number of hypergraphs (sanity of the enumeration)."""
import sys, itertools, math
from collections import defaultdict
q = int(sys.argv[1]); B = int(sys.argv[2]); NMAX = int(sys.argv[3]); doA = len(sys.argv) > 4
teams = list(itertools.product(range(q), repeat=3)); T = len(teams)
bit = [(1 << a) | (1 << (q + b)) | (1 << (2 * q + c)) for a, b, c in teams]   # agent masks
perfect = []
if doA:
    for s2 in itertools.permutations(range(q)):
        for s3 in itertools.permutations(range(q)):
            perfect.append(frozenset(teams.index((i, s2[i], s3[i])) for i in range(q)))
def Ncoef(deg):  # deg: list over v of tuple over b ; returns dict profile -> coefficient (z_inf exponent implicit)
    poly = {(0,) * B: 1}
    for dv in deg:
        new = dict(poly)
        for m, c in poly.items():
            for b in range(B):
                if dv[b]:
                    mm = m[:b] + (m[b] + 1,) + m[b + 1:]; new[mm] = new.get(mm, 0) + c * dv[b]
        poly = new
    return poly
SY = defaultdict(lambda: defaultdict(int)); CNT = defaultdict(int); violA = 0; ng = 0; slackA = 0
for n in range(NMAX + 1):
    for es in itertools.combinations(range(T), n):
        for bs in itertools.product(range(B), repeat=n):
            ng += 1
            D = [[[0] * B for _ in range(q)] for _ in range(3)]
            for e, b in zip(es, bs):
                for s in range(3): D[s][teams[e][s]][b] += 1
            key = tuple(tuple(tuple(x) for x in D[s]) for s in range(3))
            CNT[key] += 1
            Y = defaultdict(int)
            for k in range(n + 1):
                for sub in itertools.combinations(range(n), k):
                    msk = 0; ok = True
                    for i in sub:
                        if msk & bit[es[i]]: ok = False; break
                        msk |= bit[es[i]]
                    if ok:
                        m = [0] * B
                        for i in sub: m[bs[i]] += 1
                        Y[tuple(m)] += 1
            for m, c in Y.items(): SY[key][m] += c
            if doA:
                lab = dict(zip(es, bs)); A = defaultdict(int)
                for N in perfect:
                    m = [0] * B
                    for e in N:
                        if e in lab: m[lab[e]] += 1
                    A[tuple(m)] += 1
                for m, c in A.items():
                    r = q - sum(m); bound = Y[m] * math.factorial(r) ** 2
                    if c > bound: violA += 1
                    slackA += bound - c
viol = 0; eq = 0; tot = 0; minratio = 1.0; worst = None; psmax = 0
for key, S in SY.items():
    nb = [sum(key[0][v][b] for v in range(q)) for b in range(B)]
    prodD = 1
    for s in range(3):
        for v in range(q):
            for b in range(B): prodD *= math.factorial(key[s][v][b])
    cfg = 1
    for b in range(B): cfg *= math.factorial(nb[b]) ** 2
    assert CNT[key] * prodD <= cfg          # P_config(simple | D) = |G(D)| prod D! / prod (n_b!)^2 <= 1
    Ns = [Ncoef(key[s]) for s in range(3)]
    for m, lhs0 in S.items():
        lhs = lhs0 * prodD
        rhs = Ns[0].get(m, 0) * Ns[1].get(m, 0) * Ns[2].get(m, 0)
        for b in range(B): rhs *= (math.factorial(m[b]) * math.factorial(nb[b] - m[b])) ** 2
        tot += 1
        if lhs > rhs: viol += 1
        elif lhs == rhs: eq += 1
        if rhs and lhs / rhs < minratio: minratio = lhs / rhs; worst = (key, m, lhs, rhs)
print("q=%d B=%d NMAX=%d: hypergraphs %d, degree arrays %d, (D,m) pairs %d" % (q, B, NMAX, ng, len(SY), tot))
print("  Lemma 2.2 inequality: violations %d, equalities %d, smallest lhs/rhs %.4f" % (viol, eq, minratio))
if doA: print("  (2.1) count step A_m(g) <= Y_m(g) (r!)^2: violations %d (total slack %d)" % (violA, slackA))
