"""Exact-enumeration check (q=3, floating point values) of the deletion identity
E sum_{i in M} (T - T^{i->0}) = E(3q - T).  Evidence only."""
import numpy as np, itertools, sys
q = 3; perms = list(itertools.permutations(range(q)))
def T_of(A, B, C):
    best = -1; arg = None
    for s in perms:
        for t in perms:
            v = sum(A[a, s[a]] + B[a, t[a]] + C[s[a], t[a]] for a in range(q))
            if v > best: best, arg = v, (s, t)
    return best, arg
rng = np.random.default_rng(5); n = int(sys.argv[1]); reg = []; defi = []
for _ in range(n):
    A, B, C = rng.random((q, q)), rng.random((q, q)), rng.random((q, q))
    T, (s, t) = T_of(A, B, C); tot = 0.0
    for a in range(q):
        for mat, u, w in ((A, a, s[a]), (B, a, t[a]), (C, s[a], t[a])):
            old = mat[u, w]; mat[u, w] = 0.0; tot += T - T_of(A, B, C)[0]; mat[u, w] = old
    reg.append(tot); defi.append(3 * q - T)
reg, defi = np.array(reg), np.array(defi)
print('n', n, 'E sum regret', reg.mean(), '+-', reg.std() / n ** .5, ' E(3q-T)', defi.mean(), '+-', defi.std() / n ** .5,
      ' paired diff', (reg - defi).mean(), '+-', (reg - defi).std() / n ** .5)
