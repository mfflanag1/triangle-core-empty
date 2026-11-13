"""Exact brute-force check of the extreme-ray list at q <= 3: positive circuits of the balance matrix N.
A set S of off-M teams is the support of an extreme ray iff N_S has nullity 1 and its null vector is strictly positive."""
import sys, json, itertools
from flint import fmpz_mat
from rays import off_teams
q = int(sys.argv[1]); F = off_teams(q); n = len(F)
col = [[(1 if f[0] == t else 0) - (1 if f[s] == t else 0) for t in range(q - 1) for s in (1, 2)] for f in F]
r = 2 * (q - 1); found = set()
for s in range(1, r + 2):
    for S in itertools.combinations(range(n), s):
        Mx = fmpz_mat([[col[j][i] for j in S] for i in range(r)])
        X, nul = Mx.nullspace()
        if nul != 1: continue
        v = [int(X[i, 0]) for i in range(s)]
        if all(x > 0 for x in v) or all(x < 0 for x in v):
            g = 0
            from math import gcd
            for x in v: g = gcd(g, abs(x))
            found.add(tuple(sorted((F[j], abs(x) // g) for j, x in zip(S, v))))
print('q', q, 'positive circuits (exact):', len(found))
R = json.load(open(f'rays_q{q}.json'))
cddset = set(tuple(sorted((tuple(e[:3]), e[3]) for e in ray)) for ray in R)
print('cdd list size', len(cddset), 'equal:', cddset == found)
