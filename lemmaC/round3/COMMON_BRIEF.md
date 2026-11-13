# Common brief, round 3 (two-constants conjecture)

Every agent in this round reads this file. It is self-contained.

## Setting

Three sectors of q agents each. The pair values A[a,b], B[b,c], C[a,c] are 3q² independent
U[0,1] variables. A team is one agent from each sector and is worth
v(a,b,c) = A[a,b] + B[b,c] + C[a,c].

- T is the best total value over perfect assignments (q disjoint teams covering all 3q
  agents). M denotes the optimal assignment, which is almost surely unique.
- LP* is the optimum of the fractional relaxation: maximise the sum of v(e)·x_e over
  x ≥ 0 on the q³ teams, subject to the weights of the teams containing each agent summing
  to 1.
- D = LP* − T ≥ 0.
- The deficit of a pair is 1 minus its value. The deficit of a team is 3 − v, the sum of
  its three pair deficits. So 3q − T is the total deficit of M.
- The core of the associated market (payoffs p ≥ 0 on agents with p(e) ≥ v(e) for every
  team and total T) is the set of optimal dual solutions when LP* = T. It is nonempty
  exactly when D = 0.

## Conjecture

There is c > 0 with E D ≥ c·q^(1/3) for all large q. Equivalently, the per-agent loss
constants differ: 3q − LP* ~ a_LP·q^(1/3) and 3q − T ~ a_IP·q^(1/3) with a_LP < a_IP.

**What is enough for the theorem we want.** Anything that gives E D / q^(1/6) → ∞ proves
P(core nonempty) → 0, by fact 2 below.

## Facts from earlier rounds

Each has a one-line reason. None has been checked by a different model family, so
re-derive a fact before you build on it.

1. **Variance.** Write each pair value as X = 1 − exp(−E) with E ~ Exp(1). Almost
   everywhere dT/dE_i = 1{pair i is used by M}·(1 − X_i), so |grad T|² ≤ 3q − T. The
   Poincaré inequality for product Exp(1) measure (constant 4) gives Var T ≤ 4·E(3q − T).
   The same argument with pair loads of the LP optimum gives Var LP* ≤ 4·E(3q − T).
2. **Proposition R.** Applying the same inequality to D gives Var D ≤ 8·E(3q − T), so by
   Chebyshev P(core nonempty) = P(D = 0) ≤ 8·E(3q − T)/(E D)².
3. **Size of the loss.** E(3q − T) = Θ(q^(1/3)) (Frieze–Pegden–Tkocz, arXiv:1901.07167,
   Theorem 2; their costs are our deficits). Their greedy upper constant is about 4.87.
4. **Integer lower bound.** P(3q − T ≤ t) ≤ (q!)²·t^(3q)/(3q)! by a union bound over
   assignments, so 3q − T ≥ 2.15·q^(1/3) with probability 1 − exp(−Ω(q)).
5. **LP lower bound.** E(3q − LP*) ≥ 1.36·q^(1/3): the LP deficit is at least the sum over
   sector-1 agents of the smallest deficit among that agent's teams.
6. **Odd cycles.** If D > 0, the support of the LP optimum contains an odd Berge cycle of
   teams. Otherwise the support hypergraph is balanced and the vertex is integral
   (Fulkerson–Hoffman–Oppenheim 1974).
7. **Order alone is not enough.** If the pair values are so spread out that only their
   order matters (lexicographic maximisation), the LP optimum is always integral. A proof
   has to use how close the values are to each other, at scale q^(−2/3).
8. **Conditioning.** Given M and the values of its pairs, "M is optimal" is a decreasing
   event in the other pair values. By Harris, E[LP* | M, its values] is at most what fresh
   independent values would give. Treating the unused values as fresh overestimates LP*.
9. **Localisation.** For any subset S of M's teams, M restricted to S is optimal on the
   sub-market of S's agents, and D ≥ D(sub-market). Over a partition of M into blocks, D is
   at least the sum of the blocks' gaps.
10. **Pair-price certificates.** Let S be the sum of the three bipartite assignment optima.
    Then D ≤ S − T, and P(T = S) = 1/q! exactly.
11. **Scale.** P(a team's deficit < t) = t³/6 for t ≤ 1. Each agent is in q² teams, so it
    has O(1) teams with deficit below c·q^(−2/3).
12. **Symmetry.** P(core nonempty) = (q!)²·P(the identity assignment is LP-optimal), and
    P(the identity assignment is integer-optimal) = 1/(q!)².

## Data (floating point, evidence only)

| q | 16 | 24 | 32 | 40 | 48 | 56 | 64 |
|---|---|---|---|---|---|---|---|
| instances | 200 | 200 | 200 | 200 | 120 | 70 | 23 |
| E D | 0.039 | 0.069 | 0.108 | 0.127 | 0.156 | 0.188 | 0.212 |
| P(D > 0) | 0.645 | 0.875 | 0.965 | 0.985 | 1.000 | 1.000 | 1.000 |

- The standard deviation of D is about 0.05 at every size. A fit gives
  E D ≈ 0.11·q^(1/3) − 0.25.
- 3q − T ≈ 2.57·q^(1/3) and 3q − LP* ≈ 2.53·q^(1/3) at these sizes. The standard deviation
  of 3q − T is about 0.45.
- The LP optimum differs from M on one connected structure that touches about 80% of M's
  teams. Its weights are half-integral at small q and have large denominators at large q.

## Routes already tried, with the obstacle each one hit

- **Separating the constants with existing bounds.** The best rigorous bounds are 18%
  apart; the true constants differ by about 2%.
- **Improving structures of bounded size.** Their expected number goes to 0.
- **Summing gaps over blocks of M (fact 9).** Blocks of up to half the teams capture
  almost none of D.
- **Pair-price certificates (fact 10).** They exist with probability 1/q!.

A new idea that revives one of these is welcome if it says why the obstacle no longer
applies.

## Tools

- Interpreter with numpy, scipy (HiGHS) and ortools: `/Users/maxf/miniforge3/bin/python3`.
- `lemmaC/structure/mustruct.py` and `lemmaC/scaling/cg.py` compute LP*, T and D. Read
  them before writing your own.
- The machine is shared. Use at most 2 worker processes unless your task says otherwise,
  and keep outputs small (the disk is 96% full).

## Rules

- Write only inside your own subfolder of `lemmaC/round3/`. Put your results in a
  `REPORT.md` there.
- Do not edit `LEDGER.md` or `CLAIMS.md`. Do not commit and do not run `git add`.
- No paid API calls. No outreach, posting or submission.
- Label a step **complete** only if you would defend it to a referee. Label everything
  else a **gap**. Floating-point results are evidence, never proof. Cite only sources you
  have opened.
- If you come to believe the conjecture is false, say so and give the argument.
