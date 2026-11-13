# Concentration tools for Lemma C (random triangle market), 2026-11-11

Written by a Claude Code lane agent in a 60-minute time box. Not reviewed by a second model.
Numerics are floating point (scipy HiGHS) and are evidence only. Files: `poincare_check.py`,
`summarize.py`, `out/poin_q*.jsonl`.

Setting and notation as in `../../THEOREM_VANISHING.md`. Agents: three sectors of q. Pair entries
X = (X_i), i ranging over the n = 3q^2 cross-sector pairs, iid U[0,1]. A team e is a triple (one
agent per sector); v(e) = sum of its three pair entries; def(e) = 3 - v(e) >= 0 (team deficit).
P = axial polytope (each agent has load exactly 1). For w in P the *pair load* of pair i is
l_i(w) = sum of w_f over teams f containing i (in [0,1]). For an integer perfect assignment N write
1_N for its pair-load vector (0/1, 3q ones), d_N(w) = sum_{e in N} (1 - w_e) = q - sum_{e in N} w_e.
T = max over integer N of X.1_N, LP* = max over w in P of X.l(w); M = the (a.s. unique) integer
optimum, lambda = the (a.s. unique) LP optimum. D = LP* - T >= 0. Note sum_f w_f = q for w in P, so
sum_{f notin N} w_f = d_N(w).

**Main idea.** Write X_i = 1 - exp(-E_i) with E_i iid Exp(1). The product exponential measure
satisfies the Poincare inequality Var f(E) <= 4 E |grad f(E)|^2 for locally Lipschitz f
(Bobkov-Ledoux 1997, PTRF 107; constant 4 sharp in dimension 1, tensorizes). Since
dX_i/dE_i = 1 - X_i, every quantity below gets gradient weights (1 - X_i), which are tiny exactly on
the pairs that matter (used pairs sit near 1). This is the rigorous form of the coordinator's
"both values are determined locally / coupling of the two optima" remark: the gradient of D is
(l(lambda) - 1_M)(1 - X), which vanishes identically where lambda = M.

Summary of status:

| # | statement | status |
|---|---|---|
| A1 | Var T <= q, Var LP* <= q (Efron-Stein) and sub-Gaussian tails, variance proxy O(q) (Talagrand) | complete |
| A2 | Var T <= 4 E sum_{i in M}(1-X_i)^2 <= 4 E(3q - T) <= 4 zeta(2) + 8 sqrt(zeta(2) q); same for LP* | complete (cites Linusson-Wastlund / Nair-Prabhakar-Sharma) |
| A3 | Var D <= 8 E[ sum_{e in M} (1 - lambda_e) def(e) ] (= 0 contribution from core-nonempty instances) | complete |
| B1 | Margin coefficient lemma: every linear piece c = (l(w) - 1_N)/d_N(w) has ||c||_inf <= min(1, 1/d), ||c||_1 <= 6 | complete |
| B2 | Continuous signed surrogates G_kappa of mu*, with {G_kappa > 0} = {core empty} | complete |
| B3 | Var G_kappa <= 4 E[Q_kappa], explicit integrand; on {N* = M}: Q <= (2 barDelta - G)/d(w*) | complete |
| C  | Conditional theorem: E G_kappa >= m_q and E Q_kappa = o(m_q^2) imply P(core nonempty) -> 0 | complete given its hypotheses; hypotheses OPEN |
| P  | Positivity of the mean / uniform P(core empty) >= c | OPEN; obstacles recorded in section 4 |

Numerics (section 5): the Poincare integrand for mu* is measured at q = 6..24; it decays roughly
like q^{-2.4} (10 -> 24) but sits ~650x above the empirical variance and ~900x above (E mu*)^2 at
q = 24. So the Chebyshev route of Theorem C is *not* numerically closed at q <= 24; it would close
around q of order a few hundred if both trends continued (extrapolation, not evidence).

---

## 1. Concentration of T, LP*, D

**A1 (complete).** Var T <= q and Var LP* <= q. Moreover for t > 0,
P(|T - Med T| >= t) <= 4 exp(-t^2 / (12 q)), and the same for LP*.

*Proof.* Efron-Stein: Var f <= sum_i E[(f(X) - f(X^{(i)}))_+^2], X^{(i)} = X with X_i replaced by an
independent copy X'_i. Since M (optimal at X) is feasible at X^{(i)},
T(X) - T(X^{(i)}) <= (X_i - X'_i)_+ 1{i in M} <= (1 - X'_i) 1{i in M}, and X'_i is independent of
X, so E[(T - T^{(i)})_+^2] <= E[(1 - X'_i)^2] P(i in M) = P(i in M)/3. Summing, Var T <= (1/3) E|M|
= (1/3)(3q) = q. For LP* use lambda: LP*(X) - LP*(X^{(i)}) <= (1 - X'_i) l_i(lambda) and
l_i^2 <= l_i, sum_i l_i(lambda) = 3q. Tails: T and LP* are maxima of linear forms <a, X> with
a in {1_N} resp. {l(w)}, hence convex, and |<a, x - y>| <= ||a||_2 ||x - y||_2 with
||a||_2^2 <= sum_i a_i = 3q (a_i in [0,1]); Talagrand's inequality for convex Lipschitz functions of
independent [0,1] variables (Boucheron-Lugosi-Massart, Thm 7.12) gives the bound. QED.

These bounds are O(q) and are useless for D = O(1). A2-A3 replace them.

**A2 (complete).** With R_T := sum_{i in M} (1 - X_i)^2,
Var T <= 4 E R_T <= 4 E(3q - T) <= 4 zeta(2) + 8 sqrt(zeta(2) q),
and Var LP* <= 4 E sum_i l_i(lambda)^2 (1 - X_i)^2 <= 4 E(3q - LP*) <= the same bound.

*Proof.* f(E) := T(1 - e^{-E}) is a maximum of finitely many functions sum_{i in N} (1 - e^{-E_i}),
each 1-Lipschitz, so f is Lipschitz; off the null set where two assignments tie, its gradient is
(1{i in M} e^{-E_i})_i = (1{i in M}(1 - X_i))_i. Poincare gives Var T <= 4 E R_T, and
(1 - X_i)^2 <= 1 - X_i gives R_T <= sum_{i in M}(1 - X_i) = 3q - T. Same for LP* with gradient
l_i(lambda)(1 - X_i) (LP* is the maximum over the finitely many vertices of P of X.l(w); a.e.
unique maximizer), using l_i^2 <= l_i and sum_i l_i(lambda)(1 - X_i) = 3q - LP* <= 3q - T.
Upper bound on E(3q - T), by an explicit two-stage assignment. Convention v(a,b,c) =
A[a,b] + B[a,c] + C[b,c]. Stage 1: pi = optimal assignment a -> b for A. Uniform costs 1 - A are
entrywise <= Exp(1) costs under the coupling U = 1 - e^{-E} <= E, and the expected optimal assignment
cost with iid Exp(1) entries is sum_{k<=q} 1/k^2 < zeta(2) (Linusson-Wastlund 2004; Nair-Prabhakar-
Sharma 2005); so E(q - X.A-part) <= zeta(2). Stage 2: W[a,c] := (1 - B[a,c]) + (1 - C[pi(a),c]).
Given pi, the 2q^2 entries used are distinct and independent of A, so W is a q x q matrix of iid
copies of (1-U) + (1-U'). Then W <= 2 max(1-U, 1-U') = 2 sqrt(V), V := max(1-U,1-U')^2 ~ U[0,1]
iid. Take sigma optimal for V: E sum_a V[a,sigma(a)] <= zeta(2), and by Cauchy-Schwarz and Jensen
E sum_a W[a,sigma(a)] <= 2 E sqrt(q sum_a V[a,sigma(a)]) <= 2 sqrt(zeta(2) q). The teams
(a, pi(a), sigma(a)) have total deficit with expectation <= zeta(2) + 2 sqrt(zeta(2) q). QED.

Remark. The numerics give Var T = 0.44, 0.33, 0.24, 0.15 at q = 6, 10, 16, 24 (decreasing, as in the
bipartite random assignment where Var -> 0), while 4 E R_T = 7.1, 6.7, 6.1, 5.4 and E(3q - T) grows
like q^{1/3}. A proof that Var T -> 0 would need R_T -> 0, i.e. control of the largest pair
deficits used by M; not attempted.

**A3 (complete).** Var D <= 4 E sum_i (l_i(lambda) - 1_M(i))^2 (1 - X_i)^2
<= 8 E[ sum_{e in M} (1 - lambda_e) def(e) ].
In particular instances with nonempty core contribute nothing (lambda = M there), and an instance
contributes at most twice the total deficit of the M-teams that lambda disrupts, weighted by the
disruption.

*Proof.* D is a difference of two piecewise linear functions; at a.e. point its gradient in E is
(l(lambda) - 1_M)(1 - X). Split the sum. M-pairs i in team e in M: |l_i - 1| = 1 - l_i <= 1 - lambda_e,
so these terms are <= sum_e (1 - lambda_e) sum_{i in e}(1 - X_i)^2 <= sum_e (1 - lambda_e) def(e).
Non-M pairs: only teams f notin M contain them, l_i <= 1, so the terms are
<= sum_{f notin M} lambda_f sum_{i in f, i notin M} (1 - X_i) <= sum_{f notin M} lambda_f def(f).
Identity: v.lambda - v.M = sum_{f notin M} lambda_f v(f) - sum_{e in M}(1 - lambda_e) v(e) = D, and
sum_{f notin M} lambda_f = d_M(lambda) =: d, so sum_{f notin M} lambda_f def(f)
= 3d - [D + sum_e (1-lambda_e)(3 - def(e))] = sum_e (1 - lambda_e) def(e) - D. Total
<= 2 sum_e (1 - lambda_e) def(e) - D <= 2 sum_e (1 - lambda_e) def(e). QED.

What A3 does and does not give. It is the precise "smarter argument": the variance of the gap is
controlled by the deficits *on the disrupted part* of M only. Rigorously this is still
<= 8 E(3q - T) = O(sqrt q). To get Var D = o(1) one needs the disrupted M-teams to carry o(1) total
deficit weighted by (1 - lambda_e) in the squared form, which F2 makes doubtful (lambda disrupts
~0.8q teams, each with deficit ~q^{-2/3}, so the unsquared bound is ~q^{1/3}; the squared form
sum (1-X_i)^2 over those pairs is the relevant one and is plausibly o(1), but unproved). And as
THEOREM_VANISHING notes, concentration of D cannot by itself separate D from its atom at 0.

---

## 2. A continuous signed surrogate for mu* and its Poincare bound

**B1 (margin coefficient lemma, complete).** Let N be any integer perfect assignment and w != N a
point of P. Put d = d_N(w) > 0 and c = c_{N,w} := (l(w) - 1_N)/d in R^n, so that
v.(w - N)/d = <c, X>. Then
(a) for i in an N-team e: -min(1, 1/d) <= -(1 - w_e)/d <= c_i <= 0;
(b) for i not in N: 0 <= c_i = l_i(w)/d <= min(1, 1/d), and l_i(w) <= sum_{f notin N} w_f = d;
(c) sum over N-pairs of |c_i| <= 3, sum over non-N pairs of c_i <= 3, so ||c||_1 <= 6,
||c||_inf <= min(1, 1/d), ||c||_2^2 <= 6 min(1, 1/d).
*Proof.* A pair of an N-team e lies in e, so its load is >= w_e, and <= 1 (agent loads are 1).
A non-N pair is contained only in non-N teams. Sums: sum over N-pairs of (1 - l_i) <=
3 sum_e (1 - w_e) = 3d; sum over non-N pairs of l_i <= 3 sum_{f notin N} w_f = 3d. And
1 - w_e <= d. QED.

Consequence. Fix N = M: the function X -> mu_M(X) := max_{w != M vertex} <c_{M,w}, X> is convex and
sqrt(6)-Lipschitz in l2 *uniformly in q*, and sqrt(6/d)-Lipschitz on the region where the maximizer
has d(w) >= d. The only obstruction to dimension-free concentration of mu* is the dependence of M on X.

**B2 (surrogates, complete).** For kappa in [0, infinity) define
G_kappa(X) := min over integer perfect N of [ mu_N(X) + kappa (T(X) - X.1_N) ],
mu_N(X) := max over vertices w != N of P of <c_{N,w}, X>.
Then: (i) G_kappa is continuous and piecewise linear, G_kappa <= mu* := mu_M, and G_kappa = mu*
whenever mu* <= 0; (ii) G_kappa > 0 iff core empty (iff mu* > 0), a.s.; (iii) G_kappa increases
with kappa and tends to mu* pointwise as kappa -> infinity (a.s., where M is unique).
*Proof.* Each bracket is a sum of two maxima of finitely many linear forms, i.e. a max of finitely
many linear forms; G_kappa is a min of these. The N = M bracket equals mu*, so G_kappa <= mu*.
Assume M is the unique integer optimum (a.s.). For N != M: (v + mu 1_N) makes N LP-optimal only if
(v + mu 1_N).N >= (v + mu 1_N).M, i.e. mu (q - |N cap M|) >= v.M - v.N > 0, so mu_N > 0 (Lemma 2
applied to N: N optimal for v + mu 1_N iff mu >= mu_N), and the N-bracket is > 0. Hence if mu* <= 0
the minimum is attained at N = M and equals mu*; if mu* > 0 every bracket is > 0, so G_kappa > 0.
Lemma 2 gives mu* > 0 iff core empty. (iii): brackets increase in kappa; for N != M they tend to
+infinity. QED.

(G_0 alone is a poor surrogate: for N != M near-optimal, mu_N can be as small as about D/q, so
E G_0 may be much smaller than E mu*. kappa trades this off against the gradient size below.
mu* itself, kappa = infinity, can jump where M changes; its jumps occur only between values >= 0,
since at a tie M ~ M' each has ratio 0 against the other. Poincare needs a continuous function,
hence finite kappa.)

**B3 (Poincare bound, complete).** Var G_kappa <= 4 E[Q_kappa], where at a.e. X, with (N*, w*, N')
an active triple (G_kappa = <c_{N*,w*}, X> + kappa(X.1_{N'} - X.1_{N*}) on a neighbourhood),
Q_kappa = sum_i (c_{N*,w*,i} + kappa(1_{N'}(i) - 1_{N*}(i)))^2 (1 - X_i)^2.
N' = M always (a.e.; T is attained at M). On the event {N* = M} (which contains {core nonempty})
Q_kappa = sum_i c_i^2 (1 - X_i)^2 <= min{ 6 min(1, 1/d), (2 barDelta - G_kappa)/d },
d = d_M(w*), barDelta := sum_{e in M} (1 - w*_e) def(e) / d
(the disruption-weighted average deficit of M-teams, <= max_{e in M} def(e) <= 3).
*Proof.* G_kappa∘(1 - e^{-E}) is Lipschitz (B1: finitely many linear pieces with bounded coefficient
vectors, composed with a 1-Lipschitz map). Off the finite union of hyperplanes where two distinct
linear pieces agree, the active piece is locally constant and the gradient in E is the displayed
vector times (1 - X_i). Poincare gives the variance bound. For the estimate on {N* = M}: repeat the
proof of A3 with lambda replaced by w* and D replaced by v.(w* - M) = G d: N-pair terms
<= d^{-2} sum_e (1 - w_e)^2 sum_{i in e} (1 - X_i)^2 <= d^{-2} sum_e (1 - w_e) def(e); non-N terms
<= d^{-2} sum_{f notin M} w_f def(f) = d^{-2}[ sum_e (1 - w_e) def(e) - G d ]. Sum:
(2 barDelta - G)/d. The other bound is B1(c) with (1 - X_i) <= 1. QED.

Reading of B3. The variance of the signed margin is at most 4 E[(2 barDelta + |G|)/d(w*)]: small
when the most efficient improving structure is *mesoscopic* (d(w*) large) and the M-teams it
disrupts are *nearly perfect* (barDelta small). F2's picture (d ~ 0.3-0.4 q, deficits ~ q^{-2/3})
predicts Q ~ q^{-5/3}; the measured decay (section 5) is about q^{-2.4} between q = 10 and 24.

---

## 3. Conditional theorem

**Theorem C (complete given (H1)-(H2); (H1)-(H2) OPEN).** Suppose there are kappa_q in [0, infinity)
and m_q > 0 with
(H1) E G_{kappa_q} >= m_q, and
(H2) E Q_{kappa_q} = o(m_q^2).
Then P(core nonempty) <= 4 E Q_{kappa_q} / m_q^2 -> 0.
*Proof.* B2(ii) and Chebyshev with B3: P(G <= 0) <= P(|G - EG| >= m_q) <= Var G / m_q^2. QED.

A sufficient, more structural form of (H2): E[(2 barDelta + |G|)/d(w*); N* = M] = o(m_q^2) and
(12 + 2 kappa_q^2 * 6q) P(N* != M) = o(m_q^2) (on {N* != M}, Q_kappa <= 2*6 + 2 kappa^2 sum over
M-pairs and N*-pairs of (1 - X_i) <= 12 + 2 kappa^2 (3q - T + 3q - X.1_{N*}) <= 12 + 12 kappa^2 q). If m_q is bounded below (numerics: E mu* ~ 0.005-0.008 at
q = 16-24, not clearly decreasing), it suffices that d(w*) -> infinity in probability *and*
barDelta -> 0 in probability (both bounded quantities, so dominated convergence handles the rest);
actually d(w*) -> infinity alone suffices since barDelta <= 3 and |G| <= 6.
So, with kappa large, Lemma C (hence Conjecture V) reduces to two statements:
(P1) liminf E G_kappa > 0 (positivity of the mean of a continuous signed margin);
(P2) the most efficient improving structure is not of bounded size: P(d_M(w*) <= K) -> 0 for every
fixed K (on {N* = M}).
(P2) is a *local* first-moment statement (no improving-at-ratio structure of bounded size), the
kind of statement F2's heuristic count q^{k - 2r/3} -> 0 suggests; it still needs control of the
conditional law of M near a bounded set of teams, which is not available. (P1) is the genuine
difficulty and is untouched.

Caveat on (H2) as d -> infinity only: the rate matters if m_q -> 0. And with kappa finite one must
also control pieces with N* != M (gradient picks up kappa (1_{N'} - 1_{N*})(1 - X), supported on the
symmetric difference of two near-optimal assignments).

---

## 4. Positivity of the mean / uniform lower bound: what was tried (all OPEN)

Nothing in this section is proved. Recorded so the next attempt does not repeat it.

1. *E D -> infinity* is probably false: numerically D has median 0.05-0.07, roughly flat in q
   (THEOREM_VANISHING F2/F3). Not a viable target.
2. *Paley-Zygmund on D* (P(D > 0) >= (E D)^2/E D^2) needs E D^2 = O((E D)^2). The only rigorous upper
   bound on D available is D <= 3q - T, with E(3q - T) = O(sqrt q) (A2); no O(1) upper bound on E D^2
   is known. A rigorous E D >= c > 0 is also not available.
3. *Planting a bounded block* (Proposition 3 style) fails uniformly: for M to use a planted block the
   block's internal pairs must be within O(q^{-2/3}) of 1 (otherwise its agents are absorbed into
   near-perfect outside teams), so the probability is polynomially small in q at best, and making the
   block's cross pairs low costs exp(-cq). Neither gives a uniform c.
4. *Two-round exposure in exponential coordinates* (most promising, not carried out): write
   E_i = min(E'_i, E''_i) with E', E'' iid Exp(rate 1/2) - this reproduces iid U[0,1] values
   exactly. Round 1 (E') fixes a near-optimal structure; in round 2, each pair independently gets
   a fresh deficit at density 1/2 near 0. The scale problem: an improving half-integral structure
   must improve on M by exploiting ~q^{-2/3} deficits spread over Theta(q) teams (F2), so round 2
   must create a *mesoscopic* odd structure, and one needs a second-moment count of mesoscopic
   structures. This is the same obstruction as (P1).
5. *Comparison with the non-bipartite random matching LP* (where odd cycles make LP > IP whp): the
   analogy is heuristic. Here every team's value is a sum of three bipartite pair values; the
   fractional structures that win (weights 1/2 on 2k teams over k M-teams) are the analogue of odd
   cycles, but no monotone coupling between the two models was found.
6. *Exponential lower bound*: a Proposition-3 variant with cross pairs <= theta gives at best
   exp(-Cq) and its proof does not close cleanly (the rest-market optimum is not 3(q-3)); not pursued.

Conclusion: no uniform P(core empty) >= c > 0 was obtained; limsup P(core nonempty) < 1 remains open.

---

## 5. Numerics (evidence only)

`poincare_check.py`: gate-1 instances (`instance('triangle','unif',q,idx)`), MILP for M, signed mu*
by 45-step bisection on mu in [-6, 1] (the vertex just below mu* is w*), c = (l(w*) - 1_M)/d,
Q = sum_i c_i^2 (1 - X_i)^2, R_T = sum_{i in M}(1 - X_i)^2. Q is the exact Poincare integrand of
G_kappa on core-nonempty instances; on core-empty ones it is the integrand of the M-piece (G_kappa
may sit on another piece there). Run `summarize.py`.

| q | n | mean mu* | Var mu* | 4 E Q | 4EQ / Var | Var T | 4 E R_T | E(3q - T) | E[1/d] | max ||c||_1 | 4EQ/(E mu*)^2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 6 | 120 | -0.0294 | 1.63e-3 | 0.763 | 468 | 0.437 | 7.07 | 4.38 | 0.353 | 5.60 | 885 |
| 10 | 120 | -0.0090 | 4.60e-4 | 0.292 | 636 | 0.332 | 6.68 | 5.44 | 0.251 | 6.00 | 3600 |
| 16 | 100 | +0.0056 | 1.56e-4 | 0.104 | 662 | 0.241 | 6.05 | 6.45 | 0.197 | 5.43 | 3329 |
| 24 | 40 | +0.0065 | 5.66e-5 | 0.0366 | 646 | 0.228 | 5.28 | 7.43 | 0.139 | 5.35 | 878 |

Readings. (a) B1's bound ||c||_1 <= 6 is attained (6.000 at q = 10) and never exceeded.
(b) The Poincare bound is valid but loose by a factor ~500-700, stable in q (so it captures the
right *rate*, not the constant). (c) 4 E Q decays like about q^{-2} to q^{-2.4} over q = 10..24.
(d) Var T decreases with q although the rigorous bound 4 E R_T only decreases slowly.
(e) 4 E Q/(E mu*)^2 is far above 1 at q <= 24: Theorem C's Chebyshev step is not numerically
closed in this range. A Chebyshev-free route (e.g. the empirical variance itself, ~600x smaller,
gives Var/(mean)^2 ~ 1.3 at q = 24) would need a sharper variance inequality than Poincare.
(f) Check of the analytic bound Q <= (2 barDelta - G)/d (rerun with `--tag b`, `out/poinb_q*.jsonl`,
q = 6, 10, 16, n = 40, 40, 24): 0 violations; E Q = 0.162, 0.075, 0.021 against E bound = 0.55,
0.31, 0.19; barDelta = 0.74, 0.61, 0.44 (the disrupted M-teams have deficit of the same order as an
average M-team, ~2.5 q^{-2/3}); median d(w*) = 3.0, 4.4, 5.3 (growing slowly).

---

## 6. Strict list of what is proved and what is open

Proved (complete, would defend): A1, A2, A3, B1, B2, B3, Theorem C as an implication.
Cited, not re-proved: Bobkov-Ledoux Poincare constant 4 for Exp(1) products; Talagrand convex
concentration (BLM Thm 7.12); E[optimal Exp(1) assignment] = sum_{k<=q} 1/k^2.
Open: (P1) positivity of E G_kappa; (P2)/(H2) mesoscopic maximizer and the needed rate; Var T -> 0;
any uniform lower bound on P(core empty); Conjecture V.
Confidence: high for A1-A3, B1-B3 (short proofs, numerically consistent: ||c||_1 <= 6 and
Var <= 4 E Q hold in every sample); the conditional theorem is only as useful as (P1)-(P2).
