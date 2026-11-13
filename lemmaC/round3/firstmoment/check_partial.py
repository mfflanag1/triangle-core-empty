"""Exact check of the configuration-model identity with the high option (partial matchings):
E_config #{partial matchings with m_b configuration triples of bin b, all vertices distinct}
   = prod_s N_s(m) * prod_b C(n_b, m_b)^(-2),   N_s(m) = [z^m] prod_v (1 + sum_b D_{v,b} z_b).
Enumerates ALL pairings for a 3+3+3 instance with 2 bins."""
import itertools, math
from fractions import Fraction
degs = [[(2, 1), (1, 1), (1, 1)], [(1, 1), (2, 1), (1, 1)], [(1, 0), (1, 2), (2, 1)]]
B = 2; nb = [sum(d[b] for d in degs[0]) for b in range(B)]
assert all([sum(d[b] for d in degs[s]) for b in range(B)] == nb for s in range(3))
def N(deg):
    poly = {(0,) * B: 1}
    for d in deg:
        new = dict(poly)
        for m, c in poly.items():
            for b in range(B):
                if d[b]:
                    mm = m[:b] + (m[b] + 1,) + m[b + 1:]; new[mm] = new.get(mm, 0) + c * d[b]
        poly = new
    return poly
Ns = [N(degs[s]) for s in range(3)]
owners = [[[v for v, d in enumerate(degs[s]) for _ in range(d[b])] for b in range(B)] for s in range(3)]
tot = {}; cnt = 0
spaces = [list(itertools.product(itertools.permutations(range(nb[b])), repeat=2)) for b in range(B)]
for perms in itertools.product(*spaces):
    cnt += 1
    triples = [(owners[0][b][i], owners[1][b][perms[b][0][i]], owners[2][b][perms[b][1][i]], b) for b in range(B) for i in range(nb[b])]
    def rec(i, u1, u2, u3, prof):
        if i == len(triples):
            tot[prof] = tot.get(prof, 0) + 1; return
        rec(i + 1, u1, u2, u3, prof)
        a, bb, c, b = triples[i]
        if not (u1 >> a) & 1 and not (u2 >> bb) & 1 and not (u3 >> c) & 1:
            rec(i + 1, u1 | 1 << a, u2 | 1 << bb, u3 | 1 << c, prof[:b] + (prof[b] + 1,) + prof[b + 1:])
    rec(0, 0, 0, 0, (0,) * B)
ok = True
for m in sorted(set(tot) | set(Ns[0])):
    ex = Fraction(tot.get(m, 0), cnt)
    fo = Fraction(Ns[0].get(m, 0) * Ns[1].get(m, 0) * Ns[2].get(m, 0))
    for b in range(B): fo /= math.comb(nb[b], m[b]) ** 2
    ok &= ex == fo
    print(" profile", m, "exact", ex, "formula", fo, "EQUAL" if ex == fo else "DIFFER")
print("pairings enumerated:", cnt, "| all equal:", ok)
