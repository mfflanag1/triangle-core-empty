"""q = 4: exact tabulation of the 356,000 extreme rays; essentiality (pair projection outside the cone of integral directions)
tested by float LP on a random sample of fractional rays (evidence)."""
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
tab = collections.Counter(); meta = []
for y in R:
    z = zvec(y); k = sum(1 for t in z if t); integ = max(y.values()) == 1 and max(z) == 1
    keep = sum(1 for f in y if len(set(f)) == 2)
    meta.append((k, integ, max(z), len(y), keep))
    tab[(k, 'int' if integ else 'frac')] += 1
print('rays', len(R), sorted(tab.items()))
fr = [m for m in meta if not m[1]]
print('fractional k=4: zmax distribution', sorted(collections.Counter(m[2] for m in fr if m[0] == 4).items()))
print('fractional k=4: support-size distribution', sorted(collections.Counter(m[3] for m in fr if m[0] == 4).items()))
print('fractional k=4: number of pair-keeping teams in support', sorted(collections.Counter(m[4] for m in fr if m[0] == 4).items()))
print('integral k=4: number of pair-keeping teams', sorted(collections.Counter(m[4] for m in meta if m[1] and m[0] == 4).items()))
G = np.array([proj(y) for y, m in zip(R, meta) if m[1]]).T
random.seed(1); samp = random.sample([i for i, m in enumerate(meta) if not m[1] and m[0] == 4], 3000)
ess = collections.Counter()
for i in samp:
    res = linprog(np.zeros(G.shape[1]), A_eq=G, b_eq=proj(R[i]), bounds=(0, None), method='highs')
    ess[(meta[i][2], res.status != 0)] += 1
print('sample of 3000 fractional k=4 rays: (zmax, essential) counts', sorted(ess.items()))
print('essential fraction', sum(v for (z, e), v in ess.items() if e) / 3000)
