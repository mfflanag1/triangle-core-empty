# Lower-bounding the LP/IP gap E D: proof attempt, 11 November 2026

Lane agent (Claude Code), 60-minute time box. Not reviewed, and not reproduced by another model family.
Setting and notation as in `../../THEOREM_VANISHING.md` and `../SYNTHESIS_2026-11-11.md`.
Write delta_i = 1 - X_i (iid U[0,1]) for the pair deficits, def(f) = sum of the three deficits of team f,
T = integer optimum, LP* = fractional optimum, D = LP* - T, and M = the a.s. unique integer optimum.

Files in this folder:
- `first_moment.py`: exact finite-q constants for Theorems 1-2 (output in `out/first_moment.txt`).
- `local_blocks.py`: the route-3 locality probe (output in `out/local_q{16,24,32,40}.jsonl`).

## Outcome in one paragraph

**The target (E D >= c q^(1/3), or even E D -> infinity or E D >= c > 0) is NOT proved.**
What is proved:
- **Theorem 1:** a rigorous lower bound 3q - T >= (3e^(-1/3) - o(1)) q^(1/3) = (2.1496 - o(1)) q^(1/3), with
  probability 1 - e^(-Omega(q)).
- **Theorem 2:** a rigorous LP lower bound E(3q - LP*) >= 1.3628 q^(1/3).
- **Lemma 4 (FKG direction):** conditioning on the integer optimum biases every lower-bound statistic for LP* downward,
  so route 1 cannot treat the non-M values as approximately free without a quantitative reverse comparison.
- **Lemma 5 (odd Berge cycles):** every fractional LP optimum contains an odd Berge cycle (a "hole") of teams in
  its support. So the gap can only come from cyclic structure.

What is shown numerically: the gap is not captured by blocks of bounded size, even when the blocks are chosen with
the LP optimum as an oracle. At q >= 32, blocks of 6-12 M-teams capture 0% of D. Route 2 (separate constants) is
dead with present tools. Route 1 restricted to bounded constructions is dead, by the same counting that killed the
local route. Route 3 can only work through cycles whose length grows with q, which is the cavity/objective-method
problem the synthesis already named.

---

## Route 2: separate constants

### Theorem 1 (first-moment lower bound on 3q - T). Status: complete.

For every q and t >= 0,

    P(3q - T <= t) <= (q!)^2 t^(3q) / (3q)!.

Consequently, with t*_q := ((3q)!/(q!)^2)^(1/(3q)),

    E(3q - T) >= t*_q * 3q/(3q + 1),   and   t*_q / q^(1/3) -> 3 e^(-1/3) = 2.14959...

For every beta < 3e^(-1/3) we also have P(3q - T <= beta q^(1/3)) <= C sqrt(q) (beta^3 e / 27)^q -> 0
exponentially.

*Proof.*
1. A perfect assignment is a pair of bijections (X->Y, X->Z), so there are (q!)^2 of them.
2. For a fixed assignment N, 3q - v(N) = sum of 3q distinct deficits, which are iid U[0,1].
3. For n iid U[0,1] variables and any t >= 0, P(sum <= t) <= vol{x >= 0 : sum x <= t} = t^n/n!, because the
   uniform density on the cube is at most the indicator of the orthant.
4. A union bound over the (q!)^2 assignments gives the first display.
5. Hence P(3q - T <= t) <= min(1, (t/t*)^(3q)), and E(3q - T) >= int_0^{t*} (1 - (t/t*)^(3q)) dt
   = t* 3q/(3q+1).
6. Asymptotics, by Stirling: q! <= e q^(q+1/2) e^(-q) and (3q)! >= sqrt(6 pi q) (3q/e)^(3q). These give
   (q!)^2 (beta q^(1/3))^(3q) / (3q)! <= (e^2 sqrt(q) / sqrt(6 pi)) (beta^3 e / 27)^q,
   and the base is < 1 iff beta < 3 e^(-1/3).

QED.

Finite q (`out/first_moment.txt`): the bound E(3q - T)/q^(1/3) >= 1.94, 2.03, 2.06, 2.08, 2.09, 2.10 at q = 8, 16, 24,
32, 40, 48. The data are 2.45-2.60. The first moment is 16-19% below the truth.

### Theorem 2 (a dual lower bound on 3q - LP*). Status: complete.

E(3q - LP*) >= (3/4) 6^(1/3) q^(1/3) = 1.3628 q^(1/3).

*Proof.*
1. LP duality: 3q - LP* = max { sum_i r_i : r_a + r_b + r_c <= def(f) for every team f = (a, b, c) }, with r free.
2. A feasible choice: r_a = m_a := min_{f containing a} def(f) for a in sector X, and r = 0 on the other two sectors.
3. For fixed a, m_a is a minimum over q^2 teams, each a sum of three independent U[0,1] variables, so by a union
   bound P(m_a <= t) <= q^2 t^3 / 6.
4. Therefore E m_a >= int_0^inf (1 - min(1, q^2 t^3/6)) dt = (3/4)(6/q^2)^(1/3).
5. Summing over the q agents of sector X gives the claim.

QED. The data are 2.53-2.57. This bound is not needed for route 2, which needs an *upper* bound on 3q - LP*.

### Lemma 3 (route 2 cannot close with first-moment tools). Status: complete as arithmetic; its inputs are numerical.

Route 2 needs rigorous constants alpha < beta with E(3q - LP*) <= alpha q^(1/3) and E(3q - T) >= beta q^(1/3).

- Any valid upper bound satisfies alpha >= the true a_LP. The data give a_LP of about 2.53-2.57.
- The best rigorous beta here is 2.1496 (Theorem 1), and the true a_IP is about 2.56-2.60.
- So beta would have to improve on the first moment by more than 18% and land within about 1.5% of the truth,
  because the data give (a_IP - a_LP) / a_IP = E D / (3q - T) = 0.156 / 9.35 = 1.7% at q = 48.
- On the upper side, no rigorous alpha q^(1/3) bound exists at all yet. The proved bound is
  E(3q - T) <= zeta(2) + 2 sqrt(zeta(2) q) (CONCENTRATION A2).

Bipartite calibration: for the random assignment, the first moment gives cost >= 1 against the truth zeta(2) = 1.645.
First-moment-type IP lower bounds are typically 15-40% loose. A second-moment or cavity-sharp lower bound for
3D assignment at the 1% level is a research problem in its own right. **Route 2 is dead with present tools.**

---

## Route 1: direct construction from M

### Lemma 4 (conditioning on M pushes the wrong way). Status: complete.

Fix a perfect assignment m, and condition on M = m and on the values X_m of its 3q pairs. The remaining
3q^2 - 3q values X_r then have the iid U[0,1] law conditioned on

    E_m(X_m) = { X_r : X.1_N <= X.1_m for every perfect N }.

For fixed X_m this set is a down-set in X_r, since each constraint is sum_{N \ m} X_i <= sum_{m \ N} X_i.
By Harris's inequality, for every coordinatewise nondecreasing f(X_r):

    E[ f | M = m, X_m ] <= E[ f(X_m, X'_r) ],   with X'_r a fresh iid U[0,1] vector.

LP*, and v.w for every fixed w in P, are nondecreasing in X_r. So every lower bound on LP* computed with "free"
non-M values over-estimates the conditional truth.

*Proof.* An intersection of down-sets is a down-set. Harris (FKG for product measures) gives
E[f 1_E] <= E f P(E) for f increasing and E decreasing. QED.

**Consequence.** A route-1 proof must quantify how much the optimality of M depresses the non-M values that a
fractional rotation needs. Those are exactly the values that would also make *integer* rotations profitable, so the
bias is largest precisely where the construction needs fresh randomness.

### Remark (bounded constructions give o(1); heuristic, gap)

Any construction w = w(M, X) whose support moves a bounded number k of M-teams has gain at most D(G_S) for
the k-set S. THEOREM_VANISHING F2 counts the k-sets with a positive sub-gap and finds about q^(k - 2r/3)
of them, where r is the number of new pairs used. Every bounded fractional configuration has r/k > 3/2,
so the count tends to 0.

The count is a heuristic. Making it rigorous needs a whp bound such as "every M-team has deficit <= q^(-2/3+eps)",
which is not proved. The numbers agree (N3 = 0 at q >= 20). So a bounded construction can give at most E D = o(1),
and route 1 needs constructions on Theta(q) M-teams. That needs control of the joint law of M and the non-M values
at scale q^(-2/3), which is the full problem.

---

## Route 3: the odd-cycle mechanism

### Lemma 5 (fractional optima need odd Berge cycles). Status: complete, citing a classical theorem.

Let F be a set of teams, and A_F the agents x F incidence matrix. If A_F is balanced (no square submatrix of odd order
with exactly two 1s in every row and column), then every vertex of {lambda >= 0 supported on F : A_F lambda = 1} is
integral. This is Fulkerson-Hoffman-Oppenheim 1974; see Schrijver, *Combinatorial Optimization*, Ch. 83.

Consequently, if D > 0 then the support of the LP-optimal vertex contains an odd "hole": distinct teams f_1..f_k and
agents i_1..i_k, k odd and >= 3, with i_j in f_j and f_{j+1} (indices mod k) and no other incidences among them.
In particular a Berge-acyclic support (a hypertree) is always integral.

*Proof.* The LP optimum on all teams equals the LP optimum restricted to its own support F. If A_F were balanced,
that vertex would be integral and of value <= T. QED.

### Consequence for route 3 (heuristic, gap)

- At slack c q^(-2/3), the near-tight off-M hypergraph has q-independent bounded degree (STRUCTURE section 4).
- If it is locally tree-like, like a sparse random hypergraph, the expected number of Berge cycles of length <= L is
  O_L(1), bounded in q. Then the number of agents within bounded distance of a short odd cycle is O(1), not Theta(q).
- With Lemma 5, fractional rotations supported on bounded neighbourhoods could then contribute only O(1) agents at
  q^(-2/3) each, which is o(1) in total.
- So "a positive fraction of agents in local configurations with gain q^(-2/3) each" is false in the local-limit
  picture. A q^(1/3) gap must come from odd cycles whose length grows with q (at least of order log q) and that
  thread a positive fraction of M.
- This is the mechanism of the cavity/objective method: a_LP < a_IP as a difference of fixed points on the infinite
  hypertree. It is not a sum of independent local gains.

### Numerical test of locality (`local_blocks.py`, `out/local_q*.jsonl`)

Lemma 1(ii) of THEOREM_VANISHING (complete) gives D >= sum_j D(G_{S_j}) for every partition of M into blocks.
Random partitions capture nothing (STRUCTURE section 2). This probe grows blocks greedily along:
- (lam) the support of the LP optimum itself. This is an oracle choice, the best case for any local rule.
- (tight) the near-tight hypergraph {f not in M : v_f - p_lp(f) > -q^(-2/3)}.

The captured share is sum_j D(G_{S_j}) / D, averaged over 2 partitions, on 23-29 empty-core instances out of 30 per q.

Each cell gives two numbers: the fraction of instances where any of the gap is captured, then the mean captured share.

| q | lam, s = 6 | lam, s = 12 | lam, s = q/2 | tight, s = 6 | tight, s = 12 | tight, s = q/2 |
|---|---|---|---|---|---|---|
| 16 | 0.17 / 0.112 | - | 0.30 / 0.164 (s=8) | 0.09 / 0.031 | - | 0.13 / 0.071 |
| 24 | 0.13 / 0.033 | 0.26 / 0.109 (s=12=q/2) | (= s=12) | 0.00 / 0.000 | 0.17 / 0.040 | (= s=12) |
| 32 | 0.00 / 0.000 | 0.00 / 0.000 | 0.38 / 0.099 | 0.00 / 0.000 | 0.00 / 0.000 | 0.03 / 0.008 |
| 40 | 0.00 / 0.000 | 0.03 / 0.002 | 0.21 / 0.059 | 0.00 / 0.000 | 0.00 / 0.000 | 0.03 / 0.007 |

**Reading.**
- At q >= 32, even oracle-guided blocks of 6-12 M-teams capture none of the gap.
- Halves of M grown along the LP support capture 16%, then 11%, 10% and 6%, falling with q. Halves grown along the
  near-tight hypergraph capture under 1%.
- The gap is not superadditive over any partition into local pieces. This is consistent with the
  long-odd-cycle reading of Lemma 5.
- So route 3 cannot go through Lemma 1(ii). A proof must treat the LP optimum as one global object.
- Caveat: at q <= 40 the near-tight auxiliary graph has diameter about 2, so "local" and "global" cannot be separated
  cleanly. This test rules out decomposition into pieces. It cannot test a genuine local-weak-limit statement.

---

## The exact remaining gap

With Proposition R (complete given CONCENTRATION A2, unreviewed), Conjecture V follows from

    (G)  E(LP* - T) / q^(1/4) -> infinity

(q^(1/6) would suffice if E(3q - T) = O(q^(1/3)) were proved; it is not proved here or elsewhere in the lane).
Nothing proved so far gives even E D >= c > 0 for large q (THEOREM_VANISHING gap 5 is unchanged).

The most plausible proof of (G) is in the local-weak-limit framework:
1. **Local limit.** Rescale deficits by q^(2/3). Show that the near-tight hypergraph, jointly with M, converges locally
   to a Poisson-weighted 3-partite hypertree, the analogue of Aldous's PWIT. This needs a whp bound on M-team deficits
   at scale q^(-2/3+eps), which is itself open.
2. **Fixed points.** Identify a_IP and a_LP as the per-agent costs of the optimal integral and fractional perfect
   triangle factors on that limit. Show a_LP < a_IP by exhibiting distinct recursive-distributional fixed points.
   The cavity/population-dynamics computation of synthesis next step 2 is the right first check.
3. **Transfer.** Prove E(3q - T)/q^(1/3) -> a_IP and E(3q - LP*)/q^(1/3) -> a_LP. Step 3 for the LP is the harder
   half: Wastlund-type LP arguments are the closest tool.

Each of steps 1-3 is a substantial research step. None was attempted to completion here.

**A cheaper fallback that would also suffice:** Var T = O(1) and Var LP* = O(1), plus E D -> infinity at any rate.
- The Poincare integrand is E sum_{i in M} delta_i^2. It is O(q^(-1/3)) if M uses no pair of deficit above
  q^(-2/3+eps), and then Var T = o(1).
- This turns on the same whp deficit bound as step 1. It is the most tractable piece I see and is worth attacking
  first: it improves Proposition R to "E D -> infinity" by itself.

## Confidence

| claim | confidence |
|---|---|
| Theorems 1-2, Lemmas 4-5 correct | high (about 95%); short arguments, classical citations |
| Route 2 dead with present tools | high |
| Bounded-construction route 1 dead | high as evidence; the counting step is heuristic |
| Gap is non-decomposable at q <= 40 | high as a numerical statement; floating point, 30 instances per q |
| Two-constants conjecture true (E D of order q^(1/3)) | moderate, about 60%. E D / q^(1/3) still rises at q = 48 (0.043). The alternative E D = o(q^(1/3)) but -> infinity cannot be excluded, and E D -> infinity is itself unproved. |
