"""Classify the extreme rays (y, z) of the tangent cone at M = identity.
Per ray: k = number of M-teams touched, support size, zmax (denominator of the adjacent vertex), integral or not,
pair projection (U,V,W), and whether the projection lies in the cone of the integral directions (float LP, evidence)."""
import sys, json, itertools, collections
import numpy as np
from scipy.optimize import linprog
q = int(sys.argv[1])
R = [{tuple(e[:3]): e[3] for e in ray} for ray in json.load(open(f'rays_q{q}.json'))]

def zvec(y):
    z = [0] * q
    for f, c in y.items(): z[f[0]] += c
    return z

def proj(y):
    """pair projection: off-diagonal pair loads of U (AB), V (BC), W (AC), flattened; plus the released M-pair weights"""
    U = np.zeros((3, q, q))
    for (i, j, k), c in y.items():
        U[0, i, j] += c; U[1, j, k] += c; U[2, i, k] += c
    z = zvec(y)
    for s in range(3):
        for t in range(q): U[s, t, t] -= z[t]
    return U

def is_integral(y):
    return max(y.values()) == 1 and max(zvec(y)) == 1

stats = collections.Counter(); P = []; meta = []
for y in R:
    z = zvec(y); k = sum(1 for t in z if t > 0); integ = is_integral(y)
    kept = sum(1 for f in y if len(set(f)) == 2)
    meta.append(dict(k=k, supp=len(y), zmax=max(z), integ=integ, ymax=max(y.values()), kept=kept, z=tuple(sorted(z, reverse=True))))
    P.append(proj(y).reshape(-1))
P = np.array(P)
integ_idx = [i for i, m in enumerate(meta) if m['integ']]
G = P[integ_idx].T  # generators of the integral cone (projected)
ess = []
for i, m in enumerate(meta):
    if m['integ']: m['ess'] = True; continue
    res = linprog(np.zeros(G.shape[1]), A_eq=G, b_eq=P[i], bounds=(0, None), method='highs')
    m['ess'] = (res.status != 0)
# distinct projections
print('q', q, 'rays', len(R), 'integral', len(integ_idx), 'fractional', len(R) - len(integ_idx))
tab = collections.Counter((m['k'], m['integ'], m['ess'], m['supp'], m['zmax'], m['ymax'], m['z']) for m in meta)
print('k integ essential supp zmax ymax zprofile : count')
for key in sorted(tab): print(' ', key, tab[key])
byk = collections.Counter((m['k'], 'int' if m['integ'] else ('frac-ess' if m['ess'] else 'frac-red')) for m in meta)
print('by k:', sorted(byk.items()))
fe = [i for i, m in enumerate(meta) if (not m['integ']) and m['ess']]
dp = {tuple(P[i]) for i in fe}
print('essential fractional rays', len(fe), 'distinct pair projections', len(dp))
json.dump(meta, open(f'meta_q{q}.json', 'w'))
if len(sys.argv) > 2:
    seen = set()
    for i in fe:
        key = (meta[i]['k'], meta[i]['supp'], meta[i]['z'])
        if key in seen: continue
        seen.add(key); print('example', key, sorted(R[i].items()))
