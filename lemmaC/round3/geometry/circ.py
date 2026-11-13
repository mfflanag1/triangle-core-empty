"""Extreme rays of the tangent cone at M = identity for q = 3, 4 by positive-circuit enumeration with symmetry.
Rays = positive circuits of the balance matrix N (columns = off-M teams).  Symmetry group G = S_q (relabel M-teams) x S_3 (sectors).
Part A: circuits containing the pair-keeping team f0 = (0,0,1) (every circuit with a pair-keeping team is G-equivalent to one).
Part B: circuits made of teams with three distinct indices only.
Search uses float64 batch solves on small integer matrices (sign tests with tolerance); EVERY circuit found is then verified in exact
integer arithmetic (flint: nullity 1, strictly positive null vector).  Completeness therefore rests on the float sign tests."""
import sys, json, itertools, time
import numpy as np
from math import gcd
from flint import fmpz_mat
from rays import off_teams
q = int(sys.argv[1]); F = off_teams(q); n = len(F); r = 2 * (q - 1); pos = {f: i for i, f in enumerate(F)}
col = np.array([[(1 if f[0] == t else 0) - (1 if f[s] == t else 0) for t in range(q - 1) for s in (1, 2)] for f in F], float).T  # r x n
TOL = 1e-7

def search(cands, must):
    """all positive circuits within column set cands containing column `must` (or any, if must is None)"""
    out = set(); cands = list(cands)
    others = [c for c in cands if c != must]
    for s in range(1, r + 1):                      # |I| = s, circuit size s + 1
        ks = s - (0 if must is None else 1)
        it = itertools.combinations(others, ks)
        while True:
            chunk = list(itertools.islice(it, 40000))
            if not chunk: break
            idx = np.array(chunk, int).reshape(len(chunk), ks)
            if must is not None: idx = np.concatenate([np.full((len(chunk), 1), must), idx], 1)
            B = col[:, idx].transpose(1, 0, 2)                      # m x r x s
            G = B.transpose(0, 2, 1) @ B                           # m x s x s (integer entries)
            ok = np.abs(np.linalg.det(G)) > 0.5
            if not ok.any(): continue
            B = B[ok]; G = G[ok]; idx = idx[ok]
            X = np.linalg.solve(G, B.transpose(0, 2, 1) @ col[None, :, cands])   # m x s x |cands|
            res = np.abs(B @ X - col[None, :, cands]).max(1)       # m x |cands|
            hit = (res < TOL) & (X.max(1) < -TOL)
            for a, b in zip(*np.nonzero(hit)):
                out.add(frozenset(idx[a].tolist() + [cands[b]]))
    return out

def exact(S):
    S = sorted(S)
    Mx = fmpz_mat([[int(col[i, j]) for j in S] for i in range(r)]); X, nul = Mx.nullspace()
    assert nul == 1
    v = [int(X[i, 0]) for i in range(len(S))]
    assert all(x > 0 for x in v) or all(x < 0 for x in v)
    g = 0
    for x in v: g = gcd(g, abs(x))
    return tuple(sorted((F[j], abs(x) // g) for j, x in zip(S, v)))

t0 = time.time()
A = search(range(n), pos[(0, 0, 1)])
B3 = search([i for i, f in enumerate(F) if len(set(f)) == 3], None)
print('part A (contain f0):', len(A), ' part B (no pair-keeping team):', len(B3), 'secs', round(time.time() - t0, 1), flush=True)
allr = set()
perms = list(itertools.permutations(range(q))); secs = list(itertools.permutations(range(3)))
for S in A:
    ray = exact(S)
    for p in perms:
        for s in secs:
            allr.add(tuple(sorted((tuple(p[f[s[0]]], p[f[s[1]]], p[f[s[2]]]) if False else tuple(p[f[x]] for x in s), c) for f, c in ray)))
for S in B3: allr.add(exact(S))
print('q', q, 'total extreme rays', len(allr), 'secs', round(time.time() - t0, 1))
json.dump([[list(f) + [c] for f, c in ray] for ray in sorted(allr)], open(f'rays_q{q}.json' if q > 3 else f'rays_q{q}_circ.json', 'w'))
