"""q = 4: are the fractional extreme rays made only of teams sharing no pair with M essential in pair space? (float LP; evidence)
Also: essential fraction by number of pair-keeping teams in the support, and count of hub rays z = (2,1,1,1)."""
import json, collections, random
import numpy as np
from scipy.optimize import linprog
q = 4
R = [{tuple(e[:3]): e[3] for e in ray} for ray in json.load(open('rays_q4.json'))]
def zvec(y):
    z = [0] * q
    for f, c in y.items(): z[f[0]] += c
    return z
def proj(y):
    U = np.zeros((3, q, q))
    for (i, j, k), c in y.items(): U[0, i, j] += c; U[1, j, k] += c; U[2, i, k] += c
    z = zvec(y)
    for s in range(3):
        for t in range(q): U[s, t, t] -= z[t]
    return U.reshape(-1)
integ = [y for y in R if max(y.values()) == 1 and max(zvec(y)) == 1]
G = np.array([proj(y) for y in integ]).T
ess = lambda y: linprog(np.zeros(G.shape[1]), A_eq=G, b_eq=proj(y), bounds=(0, None), method='highs').status != 0
frac4 = [y for y in R if not (max(y.values()) == 1 and max(zvec(y)) == 1) and min(zvec(y)) > 0]
pure = [y for y in frac4 if all(len(set(f)) == 3 for f in y)]
c = collections.Counter((tuple(sorted(zvec(y), reverse=True)), len(y), ess(y)) for y in pure)
print('fractional rays, k=4, no pair-keeping team:', len(pure)); [print('  z profile, support, essential:', k, v) for k, v in sorted(c.items())]
prof = collections.Counter(tuple(sorted(zvec(y), reverse=True)) for y in frac4)
print('k=4 fractional rays by z profile (top 12):', prof.most_common(12))
random.seed(2); by = collections.defaultdict(list)
for y in frac4: by[sum(1 for f in y if len(set(f)) == 2)].append(y)
for nk in sorted(by):
    s = random.sample(by[nk], min(300, len(by[nk]))); print('  pair-keeping teams in support =', nk, 'count', len(by[nk]), 'essential share (sample %d)' % len(s), round(np.mean([ess(y) for y in s]), 3))
