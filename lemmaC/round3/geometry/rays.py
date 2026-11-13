"""Exact extreme rays of the tangent cone of the axial 3-index assignment polytope at M = identity.
Cone C = { y >= 0 on off-M teams (i,j,k) : load(a_t) = load(b_t) = load(c_t) for every t }  (z_t = common load).
Uses pycddlib (floating point double description); every ray is rationalised and its balance equations are checked exactly.
Completeness of the list is NOT certified by this script (see brute.py for the exact check at q <= 3). Output: rays as primitive integer vectors, saved to rays_q{q}.json."""
import sys, json, itertools, time
from fractions import Fraction
from math import gcd
from functools import reduce
import cdd

def off_teams(q):
    return [t for t in itertools.product(range(q), repeat=3) if not (t[0] == t[1] == t[2])]

def rays(q):
    F = off_teams(q); n = len(F)
    rows = []
    # equalities: alpha_t - beta_t = 0, alpha_t - gamma_t = 0  (drop t = q-1: dependent)
    for t in range(q - 1):
        for s in (1, 2):
            rows.append([0] + [(1 if f[0] == t else 0) - (1 if f[s] == t else 0) for f in F])
    neq = len(rows)
    for j in range(n):
        r = [0] * (n + 1); r[j + 1] = 1; rows.append(r)
    mat = cdd.matrix_from_array([[float(x) for x in r] for r in rows], rep_type=cdd.RepType.INEQUALITY, lin_set=set(range(neq)))
    poly = cdd.polyhedron_from_matrix(mat)
    gen = cdd.copy_generators(poly)
    out = []
    assert len(gen.lin_set) == 0
    for r in gen.array:
        if r[0] != 0: continue  # the apex
        v = [Fraction(x).limit_denominator(1000) for x in r[1:]]
        assert all(abs(float(a) - b) < 1e-9 for a, b in zip(v, r[1:]))
        den = reduce(lambda a, b: a * b // gcd(a, b), [x.denominator for x in v], 1)
        iv = [int(x * den) for x in v]
        # exact check of the balance equations
        for t in range(q):
            al = sum(c for f, c in zip(F, iv) if f[0] == t); be = sum(c for f, c in zip(F, iv) if f[1] == t); ga = sum(c for f, c in zip(F, iv) if f[2] == t)
            assert al == be == ga and min(iv) >= 0
        g = reduce(gcd, iv)
        iv = [x // g for x in iv]
        out.append({F[j]: iv[j] for j in range(n) if iv[j]})
    return F, out

if __name__ == '__main__':
    q = int(sys.argv[1]); t0 = time.time()
    F, R = rays(q)
    print('q', q, 'off teams', len(F), 'rays', len(R), 'secs', round(time.time() - t0, 1))
    json.dump([[list(k) + [v] for k, v in r.items()] for r in R], open(f'rays_q{q}.json', 'w'))
