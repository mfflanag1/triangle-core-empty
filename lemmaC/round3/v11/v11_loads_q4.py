"""Lemma V at q = 4, 5 on instances with D > 0: optimal loads unique although x is not; gradient of LP*."""
import numpy as np, itertools
from scipy.optimize import linprog
import v11_checks as V

rng = np.random.default_rng(7)
for q, want in [(4, 4), (5, 3)]:
    found = tried = 0
    spread = gerr = xspread = 0.0
    while found < want and tried < 4000:
        tried += 1
        E = rng.exponential(size=(3, q, q)); d = np.exp(-E)
        T, v, Aeq, res, dual = V.solve_lp(q, d[0], d[1], d[2])
        x = res.x
        if not np.any((x > 1e-6) & (x < 1 - 1e-6)):
            continue
        found += 1
        LP = -res.fun
        l = V.loads(q, T, x)
        A2 = np.vstack([Aeq, v]); b2 = np.concatenate([np.ones(3 * q), [LP]])
        for k in range(len(T)):                      # non-uniqueness of x itself
            c = np.zeros(len(T)); c[k] = 1
            hi = -linprog(-c, A_eq=A2, b_eq=b2, bounds=(0, None), method="highs").fun
            lo = linprog(c, A_eq=A2, b_eq=b2, bounds=(0, None), method="highs").fun
            xspread = max(xspread, hi - lo)
        for typ in range(3):
            for i in range(q):
                for j in range(q):
                    c = np.array([1.0 if ((a, b), (b, cc), (a, cc))[typ] == (i, j) else 0.0 for a, b, cc in T])
                    hi = -linprog(-c, A_eq=A2, b_eq=b2, bounds=(0, None), method="highs").fun
                    lo = linprog(c, A_eq=A2, b_eq=b2, bounds=(0, None), method="highs").fun
                    spread = max(spread, hi - lo)
                    h = 1e-6
                    Ep, Em = E.copy(), E.copy(); Ep[typ, i, j] += h; Em[typ, i, j] -= h
                    fp = -V.solve_lp(q, *np.exp(-Ep))[3].fun; fm = -V.solve_lp(q, *np.exp(-Em))[3].fun
                    gerr = max(gerr, abs((fp - fm) / (2 * h) - l[typ][i, j] * d[typ, i, j]))
    print(f"q={q}: {found} fractional instances in {tried} draws; load spread {spread:.2e}; "
          f"x spread {xspread:.2e}; LP* gradient error {gerr:.2e}")
