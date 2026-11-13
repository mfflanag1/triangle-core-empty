# The core of a random triangle market is empty with high probability: theorem and proof map

Assembled by the orchestrating session on 11 November 2026 from the round-3 files. Each
component below was written by one agent and attacked by a second, fresh agent of the same
model family, and then by a further independent one. Codex (a different model family) read
every part at least once and found no fatal error; the gaps it listed were written out and
re-read (`codex/`, `final/`, summarised in `CLAIMS.md`). No person has read any of it. Treat
every statement as a claim with that status.

## Setting

Three sectors of q agents. The 3q² pair values are independent U[0,1]. A team is one agent
from each sector and is worth the sum of its three pair values; its deficit is 3 minus its
value. T is the best total value of a perfect assignment, LP* the optimum of the
fractional relaxation, D = LP* − T >= 0. The core of the associated three-sided assignment
game is nonempty exactly when D = 0.

## Statement (E): the tree recursion has a cheap solution

Let (xi_k) be a Poisson process of intensity xi²/2 on [0, ∞). Statement (E) says:

> there is a law of a real random variable X with X =(d) min_k (xi_k − X_k' − X_k''), the
> X', X'' being independent copies of X independent of the Poisson process, and
> 9·E X < 2.302.

(E) is now proved, with 9·E X <= 2.2519: existence by a monotonicity argument in the
dispersive order, and the bound by a computer-assisted certificate. See "Status of (E)".
Numerically 9·E X = 2.2495456.

## Main theorem

Write alpha = 9·E X <= 2.2519 for the solution given by (E). Then, as q → ∞:

1. E D >= (2.302 − alpha − o(1))·q^(1/3) >= (0.050 − o(1))·q^(1/3). In particular E D → ∞.
2. D/E D → 1 in probability, and P(core nonempty) = P(D = 0) → 0.

So the core of the random triangle market is empty with probability tending to 1. Every
component of the proof has had two independent adversarial reads by fresh agents of the
authoring model family and at least one read by Codex (table at the end of `CLAIMS.md`);
no read found an error that affects a theorem, and the corrections they asked for are in
`final/ERRATA.md`. No person has read it. All bounds are asymptotic and say nothing at any
computable size.

## Unconditional theorems

**Theorem 1 (the gap is asymptotically deterministic).** For every eps > 0, with
probability 1 − O(q^(−K)) for every K, the optimal assignment uses no pair of deficit above
C_eps·q^(−1/2+eps), and every team in the support of every LP optimum has deficit below
C_eps·q^(−1/2+eps). Consequently Var T, Var LP* and Var D are O(q^(−1/6+eps)), and
P(core nonempty) → 0 whenever E D >= q^(−1/12+eps) for all large q; in particular whenever
liminf E D > 0.
File: `final/THEOREM_1.md` (the full write-up; first version `reports/pairdeficit.md`).
Check: `v3/` (holds).

**Theorem 2 (integer lower bound).** For every beta <= 2.30224, the probability that some
perfect assignment has total deficit at most beta·q^(1/3) tends to 0. In particular,
3q − T >= 2.302·q^(1/3) with probability tending to 1.
Files: `firstmoment/THEOREM_I.md` (the model with independent team deficits, with an
interval-arithmetic certificate), `sc/TRANSFER.md` (comparison theorem),
`interface/THEOREM_I_PA.md` (the join). Checks: `v6/` (independent interval enclosure
2.302028), `v7/`, `v8/` (all hold).

**Theorem 3 (LP upper bound from the tree).** Let m be any invariant (unimodular) random
integral perfect matching on the Poisson hypertree (teams' rescaled costs at each agent a
Poisson process of intensity xi²/2), with cost c(m) per agent-triple. Then
limsup E(3q − LP*)/q^(1/3) <= c(m).
Files: `lplocal/LP_LOCALITY.md`, with its section 3 replaced by the total-variation
argument in `final/TREE_AND_LOCALITY.md` (d); it uses Theorem A of
`localconv/LOCALCONV.md` (total-variation convergence of marked balls) and nothing else
from that file. Checks: `v5/`, `v4/` (hold).

**Theorem 4 (structure on the tree).** Every invariant fractional perfect matching on the
Poisson hypertree is the edge-marginal of an invariant random integral perfect matching of
the same cost. A solution of the recursion in (E) defines an invariant perfect matching of
cost 9·E X.
Files: `final/TREE_AND_LOCALITY.md` (a)-(c) (the full write-up, including the
construction of the process on the infinite-degree tree); first versions `reports/tree.md`
§2 (check `v1/`, holds) and `rde/RDE.md` §1 (check `v9/`, holds).

## Proof of the main theorem from the pieces

By Theorem 4 and (E) there is an invariant matching m on the tree with c(m) = alpha. By Theorem 3,
E(3q − LP*) <= (alpha + o(1))·q^(1/3). By Theorem 2 and 3q − T >= 0,
E(3q − T) >= 2.302·q^(1/3)·(1 − o(1)). Subtracting, E D >= (2.302 − alpha − o(1))·q^(1/3).
By Theorem 1 and Chebyshev, P(|D − E D| > t·E D) <= Var D/(t·E D)² → 0 for every fixed
t > 0, which gives D/E D → 1 in probability and, with t = 1/2, P(D = 0) → 0.

## A consequence worth stating

Theorems 2 and 3 together show that the cheap matching on the tree cannot be read off the
tree. If an invariant perfect matching m were a measurable function of the marked tree, a
thresholded version of the construction in Theorem 3 would give integral assignments of
cost (c(m) + o(1))·q^(1/3) in the finite market, so c(m) >= 2.302 by Theorem 2. Hence every
invariant matching of cost below 2.302, in particular the one defined by the recursion in
(E), uses randomness beyond the tree and its costs: the recursive tree process is not
endogenous. No file in this folder claims otherwise. (Numerically, two copies of the
process run on the same tree settle at correlation about 0.51; `redteam2/NOTES.md`.)

## What each piece imports from the literature

- Frieze, Pegden and Tkocz (arXiv:1901.07167; Discrete Applied Mathematics 287, 2020),
  Theorem 2(b): E(3n − T_n) <= c_2·n^(1/3). Used in Theorem 1 (to turn the deficit bound
  into a variance bound) and in Theorem 3 (to price the patch). Two agents read the arXiv
  source as "uniform [0,1] entries, bound in expectation"; one checker recalled a different
  form from memory. The final write-up must quote the paper.
- The Poincaré inequality for the exponential measure with constant 4 (Bobkov–Ledoux 1997,
  inequality (1.7)), tensorised.
- Harris's inequality; Fulkerson–Hoffman–Oppenheim on balanced matrices (used only in
  remarks).

## Status of (E)

- Existence (`rdeB/ROUTES.md`; check `v9/`, holds; Codex: sound after a one-line moment
  bound): on laws with increasing failure rate the map T is monotone in the dispersive
  order, which ignores shifts; the iterates Q_n = T^n(delta_0) increase to a fixed point;
  and every fixed point has cost 9·E X <= c(Q) := 3·(E T(Q) + 2·E Q) for every IFR law Q
  below it in that order, in particular for Q = Q_n.
- Certified bound (`cert/CERTIFICATE.md`, `cert/cert.py`; check `v10/`, holds; Codex source
  audit: scheme sound): iterating T ten times from delta_0, replacing each image by a less
  dispersed piecewise-exponential law and bounding everything in interval arithmetic, gives
  9·E X <= 2.2519. The exact iterates are c_5 = 2.26786 and c_10 = 2.24964 in floating
  point, with limit 2.24954561.
- A library defect found on the way: mpmath 1.3.0's interval exponential and logarithm are
  not strictly rigorous (rare failures at relative size 1e-28 to 1e-40). Both certificates
  now use verified replacements (`cert/cert.py`; `final/certI/`).
- An earlier route (`rdeA/EXISTENCE.md`, a contraction in a weighted sup norm on quantile
  functions) was not completed and is not needed.

## Evidence that the picture is right (floating point; not part of any proof)

- Model with independent team deficits: LP cost per triple rises from 2.165 (q = 16) to
  2.200 (q = 512); exact integer optimum 2.253, 2.277, 2.285, 2.311 at q = 16, 24, 32, 40
  (`pscout/results.txt`).
- Pair-additive market: the LP constant falls from 2.53 (q <= 48) through 2.47 (160) and
  2.35 (1,000) to 2.31 (2,000), heading for the tree value (`redteam2/NOTES.md`). A
  near-tight team at cost level K shares a pair with about 0.9·K²·q^(−1/3) others, so sizes
  below about 10^4 are far from the limit, and every bound here is asymptotic only.
- P(core nonempty) measured directly: 0.35 at q = 16, 0.13 at q = 24, 0.035 at q = 32,
  0 of 120 at q = 48.

## Caveats

- All checks so far are by agents of one model family.
- Theorems 2 and 3 say nothing at any computable size: the explicit thresholds in the
  proofs are astronomically large.
- The measure-theoretic framework for marked infinite-degree trees in Theorem 3 is standard
  but not written out.
- No rate is obtained for P(core nonempty) → 0 beyond what Theorem 1 gives.
