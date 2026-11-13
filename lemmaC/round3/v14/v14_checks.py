"""V14 second read of LP_LOCALITY: numerical checks (floating point, evidence only).

T1 split: d = min(d', d'') is U[0,1]; d <= u/(1-delta), d <= w/delta pointwise.
T2 rounding: exact integration over the shared uniform U; marginals = s, equal sizes.
T3 patch on tiny markets, brute force: loads of x+z are 1, LP(d) <= cost_d(x+z)
   <= cost_u(x)/(1-delta) + cost_w(z)/delta.
T4 the sub-market is genuine: E cost_w(M_R) = E Phi(N) when the rule reads u only;
   negative control: rule reads d (which depends on w).
"""
import itertools
import math
import numpy as np
from scipy.optimize import linprog
from scipy.stats import kstest

rng = np.random.default_rng(20261112)


def split(u, w, delta):
    d1 = 1 - (1 - u) ** (1 / (1 - delta))
    d2 = 1 - (1 - w) ** (1 / delta)
    return np.minimum(d1, d2), d1, d2


def t1():
    print("T1 split")
    for delta in (0.05, 0.25, 0.5, 0.9):
        u = rng.random(2_000_000)
        w = rng.random(2_000_000)
        d, d1, d2 = split(u, w, delta)
        ks = kstest(d, "uniform").statistic
        t = 0.3
        print(f"  delta={delta}: KS(d,U)={ks:.5f} (1.36/sqrt(n)={1.36/math.sqrt(len(d)):.5f})"
              f"  max(d-u/(1-delta))={np.max(d - u/(1-delta)):.2e}"
              f"  max(d-w/delta)={np.max(d - w/delta):.2e}"
              f"  P(d'>{t})={np.mean(d1 > t):.4f} vs {(1-t)**(1-delta):.4f}"
              f"  P(d''>{t})={np.mean(d2 > t):.4f} vs {(1-t)**delta:.4f}")


def team_cost(m12, m13, m23):
    return m12[:, :, None] + m13[:, None, :] + m23[None, :, :]


def rule(D, q, K):
    """A packing that is a function of the given deficits: near-tight teams, y^v_e = 1/deg_v,
    x_e = min over the three agents."""
    H = (q ** (2 / 3)) * D <= K
    deg = [H.sum(axis=(1, 2)), H.sum(axis=(0, 2)), H.sum(axis=(0, 1))]
    x = np.zeros_like(D)
    for a, b, c in zip(*np.nonzero(H)):
        x[a, b, c] = 1.0 / max(deg[0][a], deg[1][b], deg[2][c])
    return x


def slacks(x):
    return [1 - x.sum(axis=(1, 2)), 1 - x.sum(axis=(0, 2)), 1 - x.sum(axis=(0, 1))]


def pieces(s):
    """Exact decomposition of U in [0,1): list of (weight, R_1, R_2, R_3)."""
    cs = [np.concatenate([[0.0], np.cumsum(si)]) for si in s]
    brk = sorted(set([0.0, 1.0] + [float(c % 1.0) for ci in cs for c in ci]))
    out = []
    for lo, hi in zip(brk[:-1], brk[1:]):
        if hi - lo < 1e-13:
            continue
        U = (lo + hi) / 2
        R = []
        for ci in cs:
            pts = U + np.arange(0, int(math.floor(ci[-1] - U)) + 2)
            pts = pts[pts < ci[-1]]
            R.append(tuple(np.searchsorted(ci, pts, side="right") - 1))
        out.append((hi - lo, R))
    return out


PERMS = {n: list(itertools.permutations(range(n))) for n in range(0, 5)}


def best_assignment(w12, w13, w23, R):
    n = len(R[0])
    if n == 0:
        return 0.0, []
    R1, R2, R3 = (np.array(r) for r in R)
    best, arg = np.inf, None
    for sg in PERMS[n]:
        c12 = w12[R1, R2[list(sg)]].sum()
        for tau in PERMS[n]:
            c = c12 + w13[R1, R3[list(tau)]].sum() + w23[R2[list(sg)], R3[list(tau)]].sum()
            if c < best:
                best, arg = c, (sg, tau)
    sg, tau = arg
    return best, [(R1[i], R2[sg[i]], R3[tau[i]]) for i in range(n)]


def t2():
    print("T2 rounding, exact in U")
    worst_marg = worst_mean = 0.0
    sizes_ok = True
    for _ in range(300):
        q = int(rng.integers(2, 9))
        x = rule(team_cost(rng.random((q, q)), rng.random((q, q)), rng.random((q, q))), q, 2.5)
        s = slacks(x)
        Sig = s[0].sum()
        assert abs(s[1].sum() - Sig) < 1e-12 and abs(s[2].sum() - Sig) < 1e-12
        marg = [np.zeros(q) for _ in range(3)]
        EN = 0.0
        for wt, R in pieces(s):
            n = len(R[0])
            sizes_ok &= len(R[1]) == n and len(R[2]) == n and n in (math.floor(Sig), math.ceil(Sig))
            sizes_ok &= all(len(set(r)) == len(r) for r in R)
            EN += wt * n
            for i in range(3):
                for v in R[i]:
                    marg[i][v] += wt
        worst_marg = max(worst_marg, max(np.abs(marg[i] - s[i]).max() for i in range(3)))
        worst_mean = max(worst_mean, abs(EN - Sig))
    print(f"  300 instances: sizes equal and in {{floor,ceil}} and no repeats: {sizes_ok};"
          f" max|P(v in R)-s_v|={worst_marg:.2e}; max|EN-Sigma|={worst_mean:.2e}")


def t3(q=4, delta=0.25, K=2.0, runs=40):
    print(f"T3 patch, brute force, q={q}, delta={delta}")
    A = np.zeros((3 * q, q ** 3))
    for a, b, c in itertools.product(range(q), repeat=3):
        j = (a * q + b) * q + c
        A[a, j] = A[q + b, j] = A[2 * q + c, j] = 1
    worst_load = 0.0
    ok = True
    for _ in range(runs):
        u = [rng.random((q, q)) for _ in range(3)]
        w = [rng.random((q, q)) for _ in range(3)]
        d = [split(ui, wi, delta)[0] for ui, wi in zip(u, w)]
        Du, Dw, Dd = team_cost(*u), team_cost(*w), team_cost(*d)
        x = rule(Du, q, K)
        z = np.zeros_like(x)
        for wt, R in pieces(slacks(x)):
            _, teams = best_assignment(*w, R)
            for t in teams:
                z[t] += wt
        tot = x + z
        worst_load = max(worst_load, max(np.abs(sl).max() for sl in slacks(tot)))
        lp = linprog(Dd.ravel(), A_eq=A, b_eq=np.ones(3 * q), bounds=(0, None)).fun
        cd = (Dd * tot).sum()
        ub = (Du * x).sum() / (1 - delta) + (Dw * z).sum() / delta
        ok &= lp <= cd + 1e-9 <= ub + 2e-9
    print(f"  {runs} runs: max|load(x+z)-1|={worst_load:.2e}; LP(d)<=cost_d(x+z)<=bound in all: {ok}")


def t4(q=4, delta=0.5, K=2.0, runs=6000):
    print(f"T4 genuineness of the patch market, q={q}, delta={delta}")
    phi = {0: (0.0, 0.0)}
    for n in range(1, q + 1):
        v = [best_assignment(rng.random((n, n)), rng.random((n, n)), rng.random((n, n)),
                             [tuple(range(n))] * 3)[0] for _ in range(1500 if n == 4 else 6000)]
        phi[n] = (float(np.mean(v)), float(np.std(v) / math.sqrt(len(v))))
    print("  Phi(n) by simulation:", {n: (round(a, 4), round(b, 4)) for n, (a, b) in phi.items()})
    for label in ("rule reads u (as in the proof)", "rule reads d (negative control)"):
        diff = []
        for _ in range(runs):
            u = [rng.random((q, q)) for _ in range(3)]
            w = [rng.random((q, q)) for _ in range(3)]
            d = [split(ui, wi, delta)[0] for ui, wi in zip(u, w)]
            x = rule(team_cost(*(u if "reads u" in label else d)), q, K)
            ps = pieces(slacks(x))
            wts = np.array([p[0] for p in ps])
            wt, R = ps[rng.choice(len(ps), p=wts / wts.sum())]
            diff.append(best_assignment(*w, R)[0] - phi[len(R[0])][0])
        print(f"  {label}: mean[cost_w(M_R) - Phi(N)] = {np.mean(diff):+.4f}"
              f" +- {np.std(diff)/math.sqrt(runs):.4f}")


if __name__ == "__main__":
    t1()
    t2()
    t3()
    t4()
