"""Pair-sharing frequency among near-tight teams (xi <= K), against exact formulas.

Floating point / Monte Carlo: evidence only.
Formulas (s = K q^(-2/3) <= 1):
  mean degree of an agent                    = K^3/6
  mean # near-tight pair-sharing partners of a near-tight team
                                             = (9/10) K^2 q^(-1/3) (1 - 1/q)
  P(a near-tight team has >= 1 partner)      = 1 - E prod_i (1 - (s-x_i)^2/2)^(q-1),
                                               x uniform on the solid simplex of size s
  E # unordered pair-sharing pairs through a fixed agent (sharing a pair AT that agent)
                                             = q^2 (q-1) s^5 / 20  ~ K^5 q^(-1/3)/20
  P(1-ball of an agent is not a tree)       <= the line above (first moment)
"""
import sys
import numpy as np

rng = np.random.default_rng(20261111)


def simplex_noPartner(s, q, n=400000):
    # x uniform on solid simplex {x_i >= 0, sum <= s}: first three of four Dirichlet(1,1,1,1) coords times s
    g = rng.dirichlet(np.ones(4), size=n)[:, :3] * s
    return np.mean(np.prod((1 - (s - g) ** 2 / 2) ** (q - 1), axis=1))


def one(q, K):
    s = K * q ** (-2 / 3)
    A = rng.random((q, q)); B = rng.random((q, q)); C = rng.random((q, q))
    d = A[:, :, None] + B[None, :, :] + C[:, None, :]
    tight = d <= s
    n = int(tight.sum())
    nab = tight.sum(axis=2); nbc = tight.sum(axis=0); nac = tight.sum(axis=1)
    part = (nab[:, :, None] - 1) + (nbc[None, :, :] - 1) + (nac[:, None, :] - 1)
    p = part[tight]
    # sector-1 agents: pairs of tight teams sharing a pair at that agent, and non-tree 1-balls
    shareA = (nab * (nab - 1) // 2).sum(axis=1) + (nac * (nac - 1) // 2).sum(axis=1)
    return n, p.sum(), int((p > 0).sum()), shareA.sum(), int((shareA > 0).sum())


def main():
    print("q K reps | deg(sim) K^3/6 | partners/team sim formula | P(>=1 partner) sim formula |"
          " sharing pairs per agent sim formula | P(1-ball not tree) sim")
    for q, reps in ((32, 300), (48, 150), (64, 80), (96, 30), (128, 14)):
        for K in (2, 4, 6):
            s = K * q ** (-2 / 3)
            tot = np.zeros(5)
            for _ in range(reps):
                tot += np.array(one(q, K), dtype=float)
            n, ps, pp, sh, nt = tot
            f_part = 0.9 * K * K * q ** (-1 / 3) * (1 - 1 / q)
            f_any = 1 - simplex_noPartner(s, q)
            f_sh = q * q * (q - 1) * s ** 5 / 20
            print(f"{q:4d} {K} {reps:4d} | {n / (reps * q):7.3f} {K**3 / 6:7.3f} | "
                  f"{ps / n:7.4f} {f_part:7.4f} | {pp / n:6.4f} {f_any:6.4f} | "
                  f"{sh / (reps * q):8.4f} {f_sh:8.4f} | {nt / (reps * q):6.4f}")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
