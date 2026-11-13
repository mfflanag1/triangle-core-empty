"""V11 adversarial checks of THEOREM_1.md (deterministic steps and constants).

A  Theorem A Steps 1-3 by brute force at q = 3, 4 (threshold, (6.2), rotations).
B  Lemma U / Lemma G constants and arithmetic.
C  Theorem A Step 6: the displayed bound on p|G| and the constant kappa.
D  Lemma P: Poincare constant 4 on bounded Lipschitz test functions.
E  Lemma V: gradients of T and LP*, uniqueness of loads (q = 3).
F  Theorem B: dual-price chain (5.1)-(5.3), matching step (M), closing step.
"""
import itertools, math, sys
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from scipy.integrate import quad

rng = np.random.default_rng(20261112)
TOL = 1e-9
fails = []


def note(ok, msg):
    if not ok:
        fails.append(msg)
        print("FAIL:", msg)


# ---------------------------------------------------------------- A
def all_assignments(A, B, C):
    """All perfect assignments on agent lists A, B, C (as lists of teams)."""
    out = []
    for sb in itertools.permutations(B):
        for sc in itertools.permutations(C):
            out.append([(A[i], sb[i], sc[i]) for i in range(len(A))])
    return out


def cost(teams, d12, d23, d13):
    return sum(d12[a, b] + d23[b, c] + d13[a, c] for a, b, c in teams)


def check_A(q, trials):
    nrot = 0
    worst_slack = math.inf
    for _ in range(trials):
        d12, d23, d13 = rng.random((q, q)), rng.random((q, q)), rng.random((q, q))
        full = all_assignments(list(range(q)), list(range(q)), list(range(q)))
        costs = np.array([cost(t, d12, d23, d13) for t in full])
        best = costs.min()
        Mopt = full[int(costs.argmin())]
        for a0 in range(q):
            for b0 in range(q):
                uses = np.array([(a0, b0) in {(a, b) for a, b, _ in t} for t in full])
                W0 = (costs[uses] - d12[a0, b0]).min()
                W1 = costs[~uses].min()
                note(abs(best - min(W0 + d12[a0, b0], W1)) < TOL, "A: min cost != min(W0+d, W1)")
                rho = W1 - W0
                if any((a, b) == (a0, b0) for a, b, _ in Mopt):
                    note(d12[a0, b0] <= rho + TOL, "A: M uses e but d(e) > rho_e")
                else:
                    note(d12[a0, b0] >= rho - TOL, "A: M avoids e but d(e) < rho_e")
                # (6.2) and rotations
                vals, Mc_list = [], []
                for c0 in range(q):
                    A_ = [a for a in range(q) if a != a0]
                    B_ = [b for b in range(q) if b != b0]
                    C_ = [c for c in range(q) if c != c0]
                    sub = all_assignments(A_, B_, C_)
                    sc = [cost(t, d12, d23, d13) for t in sub]
                    Mc = sub[int(np.argmin(sc))]
                    cMc = min(sc)
                    vals.append(d13[a0, c0] + d23[b0, c0] + cMc)
                    Mc_list.append((Mc, cMc))
                note(abs(W0 - min(vals)) < TOL, "A: (6.2) fails")
                cmin = int(np.argmin(vals))
                n = q - 1
                for c0 in range(q):
                    Mc, cMc = Mc_list[c0]
                    for k in range(1, n + 1):
                        for seq in itertools.permutations(range(n), k):
                            teams = [Mc[j] for j in range(n) if j not in seq]
                            teams.append((Mc[seq[0]][0], b0, Mc[seq[0]][2]))
                            new = [(1, Mc[seq[0]][0], b0), (2, b0, Mc[seq[0]][2])]
                            for i in range(k - 1):
                                aj, bj, cj = Mc[seq[i + 1]][0], Mc[seq[i]][1], Mc[seq[i + 1]][2]
                                teams.append((aj, bj, cj))
                                new += [(1, aj, bj), (2, bj, cj)]
                            bk = Mc[seq[-1]][1]
                            teams.append((a0, bk, c0))
                            new += [(1, a0, bk), (2, bk, c0)]
                            # valid assignment avoiding e
                            ok = all(sorted(t[s] for t in teams) == list(range(q)) for s in range(3))
                            note(ok and len(teams) == q, "A: rotation is not an assignment")
                            note((a0, b0) not in {(a, b) for a, b, _ in teams}, "A: rotation uses e")
                            N = sum(d12[x, y] if typ == 1 else d23[x, y] for typ, x, y in new)
                            cr = cost(teams, d12, d23, d13)
                            note(cr <= d13[a0, c0] + cMc + N + TOL, "A: rotation cost bound fails")
                            note(W1 <= d13[a0, c0] + d23[b0, c0] + cMc + N + TOL, "A: W1 bound fails")
                            if c0 == cmin:
                                note(rho <= N + TOL, "A: rho_e > N at the minimiser c")
                                worst_slack = min(worst_slack, N - rho)
                            nrot += 1
    print(f"A  q={q}: {trials} instances, {nrot} rotations checked, min(N - rho_e) = {worst_slack:.3e}")


# ---------------------------------------------------------------- B
def check_B():
    e1 = 1 - math.exp(-1)
    delta = 1 - 0.3 / e1
    rate = delta ** 2 * e1 / 2
    print(f"B  delta = {delta:.6f}, rate = {rate:.6f} (doc 0.52540, 0.08724)")
    note(rate >= 0.079 and 0.079 * 13 >= 1.0, "B: U(a) rate")
    pi = 1 - math.exp(-3)
    kl = 0.9 * math.log(0.9 / pi) + 0.1 * math.log(0.1 / (1 - pi))
    print(f"B  KL(0.9 || 1-e^-3) = {kl:.6f} (doc 0.020879); exponent {0.01 + 1/600 - 0.0208:.6f}")
    note(kl >= 0.0208 and 0.01 + 1 / 600 - 0.0208 <= -0.009, "B: U(b) rate")
    # (1-e^-t)/t decreasing on (0,1]; pi_m >= (1-e^-1) p m for p m <= 1
    for p in np.concatenate([np.logspace(-8, np.log10(0.5), 400), [0.5]]):
        m0 = math.floor(1 / p)
        note(m0 >= 2 or abs(p - 0.5) < 1e-12 and m0 == 2, "B: floor(1/p) < 2")
        for m in {1, 2, max(1, m0 // 2), m0}:
            if m * p <= 1:
                note(1 - (1 - p) ** m >= e1 * p * m - 1e-15, "B: pi_m bound")
        note(0.3 * p * m0 >= 0.15 - 1e-15, f"B: 0.3 p floor(1/p) < 0.15 at p={p}")
        mb = math.ceil(3 / p)
        note(1 - (1 - p) ** mb >= pi - 1e-15, "B: pi_{m_b} bound")
    # exact binomial tail vs Chernoff (C1) at a few points
    from scipy.stats import binom
    for n, pp in [(50, 0.1), (200, 0.3), (1000, 0.05)]:
        for dl in [0.2, 0.5, 0.9]:
            note(binom.cdf(math.floor((1 - dl) * n * pp), n, pp) <= math.exp(-dl ** 2 * n * pp / 2) + 1e-15, "B: C1")
    # count of vertex-disjoint position sets <= n^(2m)
    for n in [3, 7, 20]:
        for m in range(1, n + 1):
            note(math.comb(n, m) ** 2 * math.factorial(m) <= n ** (2 * m), "B: set count")
    # (Q2) step of Theorem B: lam^(h+1) >= 0.15 q whenever q^(2eps) >= (20/3)^h
    for eps in [0.01, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.49]:
        h = math.ceil(1 / (2 * eps))
        lq = h * math.log(20 / 3) / (2 * eps)          # log q at equality in (Q2)
        for lq_ in [lq, lq * 1.5, lq * 4]:
            lhs = (h + 1) * (math.log(0.15) + 2 * eps * lq_)
            note(lhs >= math.log(0.15) + lq_ - 1e-9, f"B: (Q2) step fails eps={eps}")
    print("B  constants and Lemma G arithmetic checked")


# ---------------------------------------------------------------- C
def check_C():
    print("C  Theorem A Step 6: guaranteed p|G| versus the document's display and kappa*q^(2eps)")
    bad_display = bad_kappa_doc_inputs = bad_kappa_exact = 0
    worst = math.inf
    for eps in [0.05, 0.1, 0.125, 1 / 6, 0.2, 0.25, 0.3, 0.4, 0.45]:
        rA = math.ceil(1 / (2 * eps) - 1e-12) - 1
        kappa = 0.3 ** rA / (2 * 4 ** (rA + 1))
        for lq in np.linspace(math.log(1e3), math.log(1e60), 300):
            q = math.exp(lq)
            p = q ** (-1 + 2 * eps) / 2
            for n in (q - 1, q / 2):               # true n, and the weakest n the text assumes
                npv = n * p
                lam = 0.3 * npv
                G = min(lam ** rA * npv / 2, 0.15 * n)     # guaranteed |G| on Good
                pG = p * G
                display = min(0.3 ** rA * npv ** (rA + 1) / 2, 0.15 * npv)
                if pG < display * (1 - 1e-12):
                    bad_display += 1
                if n == q / 2 and pG < kappa * q ** (2 * eps) * (1 - 1e-12):
                    bad_kappa_doc_inputs += 1
                if n == q - 1:
                    worst = min(worst, pG / (kappa * q ** (2 * eps)))
                    if pG < kappa * q ** (2 * eps):
                        bad_kappa_exact += 1
    print(f"   grid points where guaranteed p|G| < document's displayed minimum: {bad_display}")
    print(f"   with only n*p >= q^(2eps)/4 (n = q/2): p|G| < kappa q^(2eps) at {bad_kappa_doc_inputs} points")
    print(f"   with n = q-1 exactly: p|G| < kappa q^(2eps) at {bad_kappa_exact} points; min ratio {worst:.4f}")
    # the fixed constant kappa/2 under the document's own inputs
    ok = True
    for eps in [0.05, 0.1, 0.125, 1 / 6, 0.2, 0.25, 0.3, 0.4, 0.45]:
        rA = math.ceil(1 / (2 * eps) - 1e-12) - 1
        kappa = 0.3 ** rA / (2 * 4 ** (rA + 1))
        for lq in np.linspace(math.log(1e3), math.log(1e60), 300):
            q = math.exp(lq); p = q ** (-1 + 2 * eps) / 2; n = q / 2
            pG = p * min((0.3 * n * p) ** rA * n * p / 2, 0.15 * n)
            ok &= pG >= (kappa / 2) * q ** (2 * eps) * (1 - 1e-9)
    print(f"   kappa/2 suffices under n*p >= q^(2eps)/4: {ok}")


# ---------------------------------------------------------------- D
def check_D():
    worst = 0.0
    tests = []
    for a in [0.1, 0.3, 0.45, 0.49, 0.499]:
        for Bcap in [5.0, 20.0, 60.0, 200.0]:
            tests.append((f"min(exp({a}x),exp({a}*{Bcap}))",
                          lambda x, a=a, Bc=Bcap: math.exp(a * min(x, Bc)),
                          lambda x, a=a, Bc=Bcap: a * math.exp(a * x) if x < Bc else 0.0, Bcap))
    for c in [0.5, 2.0, 10.0]:
        tests.append((f"min(x,{c})", lambda x, c=c: min(x, c), lambda x, c=c: 1.0 if x < c else 0.0, c))
    tests.append(("sin x", math.sin, math.cos, None))
    tests.append(("1-exp(-x)", lambda x: 1 - math.exp(-x), lambda x: math.exp(-x), None))
    for name, g, dg, brk in tests:
        pts = [brk] if brk else None
        I = lambda h: quad(lambda x: h(x) * math.exp(-x), 0, 400, points=pts, limit=400)[0]
        m1, m2 = I(g), I(lambda x: g(x) ** 2)
        var, en = m2 - m1 ** 2, I(lambda x: dg(x) ** 2)
        ratio = var / en
        worst = max(worst, ratio)
        note(ratio <= 4 + 1e-6, f"D: Poincare ratio {ratio} > 4 for {name}")
    print(f"D  largest Var/E g'^2 over {len(tests)} bounded Lipschitz test functions: {worst:.4f} (<= 4)")


# ---------------------------------------------------------------- E, F helpers
def teams_of(q):
    return [(a, b, c) for a in range(q) for b in range(q) for c in range(q)]


def solve_lp(q, d12, d23, d13):
    T = teams_of(q)
    v = np.array([3 - d12[a, b] - d23[b, c] - d13[a, c] for a, b, c in T])
    Aeq = np.zeros((3 * q, len(T)))
    for k, (a, b, c) in enumerate(T):
        Aeq[a, k] = Aeq[q + b, k] = Aeq[2 * q + c, k] = 1
    res = linprog(-v, A_eq=Aeq, b_eq=np.ones(3 * q), bounds=(0, None), method="highs")
    dual = linprog(np.ones(3 * q), A_ub=-Aeq.T, b_ub=-v, bounds=(None, None), method="highs")
    return T, v, Aeq, res, dual


def loads(q, T, x):
    l12, l23, l13 = np.zeros((q, q)), np.zeros((q, q)), np.zeros((q, q))
    for k, (a, b, c) in enumerate(T):
        l12[a, b] += x[k]; l23[b, c] += x[k]; l13[a, c] += x[k]
    return l12, l23, l13


def check_E(q=3, trials=6):
    gapcount = 0
    maxerr = 0.0
    spread = 0.0
    for _ in range(trials):
        E = rng.exponential(size=(3, q, q))
        d = np.exp(-E)

        def vals(dd):
            T, v, Aeq, res, dual = solve_lp(q, dd[0], dd[1], dd[2])
            full = all_assignments(list(range(q)), list(range(q)), list(range(q)))
            cs = [cost(t, dd[0], dd[1], dd[2]) for t in full]
            return 3 * q - min(cs), -res.fun, full[int(np.argmin(cs))], T, v, Aeq, res, dual

        Tval, LPval, M, T, v, Aeq, res, dual = vals(d)
        note(abs(LPval - dual.fun) < 1e-8, "E: strong duality")
        gapcount += LPval - Tval > 1e-7
        l = loads(q, T, res.x)
        # uniqueness of loads: max and min of each load over the optimal face
        A2 = np.vstack([Aeq, v]); b2 = np.concatenate([np.ones(3 * q), [LPval]])
        for typ in range(3):
            for i in range(q):
                for j in range(q):
                    cvec = np.array([1.0 if ((typ == 0 and (a, b) == (i, j)) or (typ == 1 and (b, c) == (i, j))
                                             or (typ == 2 and (a, c) == (i, j))) else 0.0 for a, b, c in T])
                    hi = -linprog(-cvec, A_eq=A2, b_eq=b2, bounds=(0, None), method="highs").fun
                    lo = linprog(cvec, A_eq=A2, b_eq=b2, bounds=(0, None), method="highs").fun
                    spread = max(spread, hi - lo)
        # finite differences
        Mpairs = [set((a, b) for a, b, c in M), set((b, c) for a, b, c in M), set((a, c) for a, b, c in M)]
        h = 1e-6
        for typ in range(3):
            for i in range(q):
                for j in range(q):
                    Ep, Em = E.copy(), E.copy()
                    Ep[typ, i, j] += h; Em[typ, i, j] -= h
                    Tp, Lp = vals(np.exp(-Ep))[:2]; Tm, Lm = vals(np.exp(-Em))[:2]
                    gT = (Tp - Tm) / (2 * h); gL = (Lp - Lm) / (2 * h)
                    maxerr = max(maxerr, abs(gT - ((i, j) in Mpairs[typ]) * d[typ, i, j]),
                                 abs(gL - l[typ][i, j] * d[typ, i, j]))
    note(maxerr < 1e-4, f"E: gradient formula error {maxerr}")
    note(spread < 1e-6, f"E: optimal loads not unique, spread {spread}")
    print(f"E  q={q}: {trials} instances ({gapcount} with D>0), max gradient error {maxerr:.2e}, "
          f"max load spread on the optimal face {spread:.2e}")


def check_F(q, s, r, trials):
    """Deterministic content of Theorem B Steps 1-6 for arbitrary s (no event needed)."""
    closed = 0
    gap = 0
    minmargin = math.inf
    for _ in range(trials):
        d12, d23, d13 = rng.random((q, q)), rng.random((q, q)), rng.random((q, q))
        T, v, Aeq, res, dual = solve_lp(q, d12, d23, d13)
        x, u = res.x, dual.x
        note(abs(-res.fun - dual.fun) < 1e-8, "F: strong duality")
        full_cost = None
        X = x.reshape(q, q, q)
        supp = X > 1e-9
        note(np.all(Aeq.T @ u >= v - 1e-8), "F: dual infeasible")
        note(np.all(np.abs((Aeq.T @ u - v)[x > 1e-9]) < 1e-7), "F: complementary slackness")
        gap += int(np.any((x > 1e-6) & (x < 1 - 1e-6)))
        # deficits arranged per sector role: sector 'mid' with outer positions
        # role 2 (as in the document): positions (a,c), mid b
        D = {  # D[mid][0][x, z] = deficit between first outer agent x and mid agent z; [1] second outer
            1: (d12, d23.T),          # mid b: d(a,b), d(c,b)
            0: (d12.T, d13.T),        # mid a: d(b,a), d(c,a)
            2: (d13, d23),            # mid c: d(a,c), d(b,c)
        }
        Xr = {1: X.transpose(0, 2, 1), 0: X.transpose(1, 2, 0), 2: X}   # [outer1, outer2, mid]
        off = {0: 0, 1: q, 2: 2 * q}
        for f in np.argwhere(supp):
            a, b, c = f
            df = d12[a, b] + d23[b, c] + d13[a, c]
            sets = {}
            for mid in (0, 1, 2):
                o1, o2 = {1: (a, c), 0: (b, c), 2: (a, b)}[mid]
                z0 = (a, b, c)[mid]
                d1, d2 = D[mid]
                base = u[off[mid] + z0] + d1[o1, z0] + d2[o2, z0]
                Bi = set(np.nonzero(d1[o1, :] + d2[o2, :] < s)[0])
                for i in range(r + 1):
                    for z in Bi:                                          # (5.3)
                        m = u[off[mid] + z] - (base - (i + 1) * s)
                        minmargin = min(minmargin, m)
                        note(m > -1e-8, "F: (5.3) fails")
                    if i == r or not Bi:
                        break
                    y = Xr[mid][:, :, sorted(Bi)].sum(axis=2)
                    H = csr_matrix((y > 1e-9).astype(int))
                    match = maximum_bipartite_matching(H, perm_type="column")
                    S = [(row, col) for row, col in enumerate(match) if col >= 0]
                    note(len(S) >= len(Bi), f"F: matching smaller than |B_i| ({len(S)} < {len(Bi)})")
                    Bi = set()
                    for (x1, x2) in S:
                        Bi |= set(np.nonzero(d1[x1, :] + d2[x2, :] < s)[0])
                sets[mid] = Bi
            if any(supp[a2, b2, c2] for a2 in sets[0] for b2 in sets[1] for c2 in sets[2]):
                closed += 1
                note(df < 3 * (r + 1) * s + 1e-8, f"F: closing bound fails d(f)={df}")
    print(f"F  q={q}, s={s}, r={r}: {trials} instances ({gap} fractional), closing step reached {closed} times, "
          f"min margin in (5.3) = {minmargin:.3e}")


if __name__ == "__main__":
    check_A(3, 40)
    check_A(4, 6)
    check_B()
    check_C()
    check_D()
    check_E()
    check_F(6, 0.5, 2, 40)
    check_F(8, 0.35, 3, 25)
    check_F(8, 0.2, 1, 25)
    print("FAILURES:", len(fails))
    for m in sorted(set(fails)):
        print("  ", m, "x", fails.count(m))
