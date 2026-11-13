"""Post-gate extension (2026-11-11): exact certificates of empty cores, triangle, q = 3.

Instances: gate1 seeds (SeedSequence([20261111, 0, 0, q, idx])), pair values ROUNDED to
two decimals (so the certified instance is a readable rational instance; it is a different
instance from the float one, certified on its own).  For each:
  * fractional primal lam (exact rationals) feasible for  sum_{e ni i} lam_e <= 1, lam >= 0;
  * dual y >= 0 with y(e) >= v(e) for EVERY team e and sum y = sum lam_e v(e)  (so LP* is exact);
  * integer optimum T by exhaustive enumeration of all perfect 3D matchings ((q!)^2 of them;
    values are nonnegative so some perfect matching is optimal among all integer matchings);
  * gap LP* - T > 0 exactly.  The core {p >= 0: p(N) = T, p(e) >= v(e)} is then empty, because
    any core point would be dual feasible with sum T < LP*.
Basis found by float HiGHS, then everything recomputed and checked in Fractions.
"""
import itertools, json, os, sys
from fractions import Fraction as F
import numpy as np
import sympy
from scipy.optimize import linprog
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gate1'))
from gate1 import instance, structure


def exact_instance(q, idx):
    A, B, C = instance('triangle', 'unif', q, idx)
    r = lambda X: [[F(round(float(x) * 100), 100) for x in row] for row in X]
    return r(A), r(B), r(C)


def gj(M, b):
    sol, params = M.gauss_jordan_solve(b)  # raises ValueError if inconsistent
    return sol.subs({p: 0 for p in params})  # any particular solution; checked exactly below


def certify(q, idx):
    A, B, C = exact_instance(q, idx)
    teams = list(itertools.product(range(q), repeat=3))
    v = {e: A[e[0]][e[1]] + B[e[0]][e[2]] + C[e[1]][e[2]] for e in teams}
    players = lambda e: (e[0], q + e[1], 2 * q + e[2])
    n = 3 * q
    vf = np.array([float(v[e]) for e in teams])
    _, M = structure(*(np.array([[float(x) for x in row] for row in X]) for X in (A, B, C)))
    lp = linprog(-vf, A_ub=M.T.tocsr(), b_ub=np.ones(n), bounds=(0, None), method='highs-ds')
    if np.max(np.minimum(lp.x, 1 - lp.x)) < 1e-6:
        return None  # float LP integral: not a candidate
    y_f = -lp.ineqlin.marginals
    S = [teams[k] for k in range(len(teams)) if lp.x[k] > 1e-9]
    P = [i for i in range(n) if y_f[i] > 1e-9]
    # exact dual: y(e) = v(e) on S, y_i = 0 off P
    Md = sympy.Matrix([[1 if i in players(e) else 0 for i in P] for e in S])
    yP = gj(Md, sympy.Matrix([sympy.Rational(v[e].numerator, v[e].denominator) for e in S]))
    y = [F(0)] * n
    for i, val in zip(P, yP):
        y[i] = F(int(sympy.fraction(val)[0]), int(sympy.fraction(val)[1]))
    # exact primal: sum_{e in S, e ni i} lam_e = 1 for i in P
    Mp = sympy.Matrix([[1 if i in players(e) else 0 for e in S] for i in P])
    lamS = gj(Mp, sympy.ones(len(P), 1))
    lam = {e: F(int(sympy.fraction(x)[0]), int(sympy.fraction(x)[1])) for e, x in zip(S, lamS)}
    # exact checks
    assert all(x >= 0 for x in lam.values())
    load = [sum((lam[e] for e in S if i in players(e)), F(0)) for i in range(n)]
    assert all(l <= 1 for l in load)
    assert all(yi >= 0 for yi in y)
    assert all(y[a] + y[q + b] + y[2 * q + c] >= v[(a, b, c)] for (a, b, c) in teams)
    LP = sum(lam[e] * v[e] for e in S)
    assert LP == sum(y), (LP, sum(y))
    best, arg = None, None
    for pb in itertools.permutations(range(q)):
        for pc in itertools.permutations(range(q)):
            val = sum(v[(a, pb[a], pc[a])] for a in range(q))
            if best is None or val > best:
                best, arg = val, [(a, pb[a], pc[a]) for a in range(q)]
    gap = LP - best
    return dict(q=q, idx=idx, A=A, B=B, C=C, lam=lam, y=y, LP=LP, T=best, T_matching=arg, gap=gap)


def fmt(d):
    s = lambda x: str(x)
    out = [f"## q={d['q']}, gate1 seed idx={d['idx']} (values rounded to hundredths)", '']
    for name in 'ABC':
        out.append(f"{name} = " + json.dumps([[s(x) for x in row] for row in d[name]]))
    out.append("teams (a,b,c) indexed by sector-1 agent a, sector-2 agent b, sector-3 agent c; "
               "v(a,b,c) = A[a][b] + B[a][c] + C[b][c]")
    out.append("fractional primal lam: " + ", ".join(f"{e}:{s(x)}" for e, x in d['lam'].items()))
    out.append("dual y (sector 1 | 2 | 3): " + " | ".join(", ".join(s(x) for x in d['y'][k*d['q']:(k+1)*d['q']]) for k in range(3)))
    out.append(f"LP* = {s(d['LP'])}   integer optimum T = {s(d['T'])} at {d['T_matching']}   "
               f"gap LP* - T = {s(d['gap'])} > 0  (relative {float(d['gap']/d['T']):.5f})")
    out.append('')
    return '\n'.join(out)


if __name__ == '__main__':
    q = int(sys.argv[1]); want = int(sys.argv[2])
    got, idx, lines = 0, 0, []
    while got < want:
        try:
            d = certify(q, idx)
            if d is not None and d['gap'] > 0:
                lines.append(fmt(d)); got += 1
        except (AssertionError, ValueError) as ex:
            print(f'idx {idx}: basis recovery failed ({type(ex).__name__}); skipped', file=sys.stderr)
        idx += 1
    print('\n'.join(lines))
    print(f'(scanned idx 0..{idx - 1})')
