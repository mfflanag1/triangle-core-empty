# Adversarial check of Proposition R inputs, 11 November 2026

Written by a Claude Code subagent (same model family as the authors of the checked notes, so this
is not the cross-family review asked for in SYNTHESIS next step 1). Time box 40 minutes. No API
spend, no LEDGER edit. Checked: `../SYNTHESIS_2026-11-11.md` (Proposition R) and
`../concentration/CONCENTRATION.md` (A1-A3).

Sources read in full:
- FPT = A. Frieze, W. Pegden, T. Tkocz, "On random multi-dimensional assignment problems",
  arXiv:1901.07167 (v1, 22 Jan 2019, the only arXiv version; LaTeX source read), published as
  Discrete Appl. Math. 287 (2020) 1-9, doi:10.1016/j.dam.2020.07.013. The DAM version was not
  seen (paywalled); theorem numbers below are those of arXiv v1.
- BL = S. Bobkov, M. Ledoux, "Poincaré's inequalities and Talagrand's concentration phenomenon for
  the exponential distribution", Probab. Theory Relat. Fields 107 (1997) 383-400 (author PDF,
  pages 383-388 read).

## Verdicts

| Item | Verdict |
|---|---|
| 1. FPT gives E(3q - T) = Theta(q^(1/3)) for our model | **Holds.** Model matches exactly; Thm 2(b) is the needed upper bound. Two caveats on the paper, neither affecting the order. |
| 2a. U = 1 - e^{-E} is U[0,1] | **Holds.** |
| 2b. Gradient of T, \|grad T\|^2 <= 3q - T | **Holds.** |
| 2c. Gradient of LP*, \|grad LP*\|^2 <= 3q - LP* | **Holds with fix.** The LP optimum need not be unique in team weights, even a.s. Uniqueness is needed only for the optimal pair-load vector, and that is a.s. unique. |
| 2d. Poincaré, Var f <= 4 E\|grad f\|^2 for product Exp(1) | **Holds.** Citation should be BL (1.7) / Lemma 2.1 plus tensorisation; conditions met (f bounded and globally Lipschitz). |
| 2e. Var T, Var LP* <= 4 E(3q - T) | **Holds.** |
| 2f. Chebyshev step and Proposition R thresholds (q^(1/4); q^(1/6) under FPT) | **Holds.** Constant 16 can be 8. |

## 1. Frieze-Pegden-Tkocz

**Model (FPT eq. (costs), Section 1).** Minimise sum_i W_{i,sigma(i),tau(i)} over pairs of
permutations sigma, tau of [n], where W_{i,j,k} = a_{i,j} + b_{i,k} + c_{j,k} and all 3n^2 of the
a, b, c are independent U[0,1]. The objects are costs and the problem is a minimisation; Z_n is the
optimum value. It is three-dimensional with pairwise-decomposable costs. The intro calls it
"(planar)", but the formulation (two permutations, each index used once) is the axial 3-index
problem. That is our problem, and the word is a slip.

**Our model maps onto it exactly.** In our model T = max over perfect assignments of
sum over teams of [A + B + C], with 3q^2 iid U[0,1] pair values. Then
3q - T = min over assignments of sum over teams of [(1-A) + (1-B) + (1-C)], the entries 1 - U are
again iid U[0,1], and the sector-to-index map (sector 1 -> i, sector 2 -> j, sector 3 -> k; our
X-Y pair -> a_{ij}, X-Z -> b_{ik}, Y-Z -> c_{jk}) is a bijection. So 3q - T has the law of Z_q,
identically for each q, and E(3q - T) = E Z_q.

**Statements (arXiv v1 numbering; a shared counter, so Frieze-Sorkin is Theorem 1).**
- **Theorem 2.** There are constants c_1, c_2 with (a) E(Z_n) >= c_1 n^{1/3}, and (b) GREEDY(n - n^{1/4})
  with an arbitrary completion finds a solution of expected value <= c_2 n^{1/3}. Because the
  greedy output is feasible, (b) gives E Z_n <= c_2 n^{1/3}.
- **Lemma 3.** E(W_n) ~ c_1 n^{-2/3}, where W_n is the minimum weight in a fixed plane.
  (a) follows from E Z_n >= n E W_n.
- **Lemma 4.** For m <= n - n^{1/4}, E(W_m) <~ c_1 (n - m + 1)^{-2/3} for the weight added at greedy
  round m. Summing (eq. (greedycost)) gives about 3 c_1 n^{1/3} plus at most 3n^{1/4} for the
  completion.
- **Theorem 5** (sketch only). In d dimensions with (d-1)-wise decomposable costs the order is
  n^{1/d}, with c_d = (d!)^{1/d} Gamma(1 + 1/d). The pairwise 4-dimensional model is "not
  considered".

**Constants.** The lower-bound constant is asymptotic, c_1 = 6^{1/3} Gamma(4/3) = 1.6227. The
greedy constant is asymptotic, c_2 ~ 3 c_1 = 4.868. The paper's final comments say it is a
(3 + o(1))-approximation whp.
*Typo:* Lemma 3's statement prints c_1 = 6^{1/3} Gamma(3/4) (= 2.227), but its own proof computes
int_0^inf e^{-x^3/6} dx = 6^{1/3} Gamma(4/3). The correct value is Gamma(4/3). Anyone quoting the
constant should use 1.6227.
Our data (3q - T ~ 2.55 q^{1/3}) lie between 1.62 and 4.87, as they should.

**Rigour caveats (adversarial).**
- (a) is robust. An elementary non-asymptotic version needs nothing from the paper: in plane i,
  P(min_{j,k} W_{ijk} <= t) <= q^2 t^3 / 6 for t <= 1, so E min >= (3/4)(6/q^2)^{1/3}. Hence
  E(3q - T) >= (3/4) 6^{1/3} q^{1/3} = 1.363 q^{1/3} for all q >= 3. Proposition R does not use
  the lower bound at all.
- (b) is the bound Proposition R uses, and its proof (Lemma 4) is terse:
  - typo p_m = L/nu_m, which should be L/nu_m^{2/3};
  - the conditioning in eq. (one) is handled in one line;
  - the claim that the c_{jk} outside the exposed sets X_l are "unconditioned by the history" is
    asserted rather than proved.

  These are presentation gaps in a refereed paper, not evident errors. Cite it as Theorem 2(b) of
  the DAM paper. Have a human check the DAM text for renumbering, and note that the
  q^{1/6} version of Proposition R rests on this lemma.
- **Independence from FPT.** Proposition R's q^{1/4} version needs only A2's own bound,
  E(3q - T) <= zeta(2) + 2 sqrt(zeta(2) q). I re-checked that two-stage construction: given pi,
  the 2q^2 entries used in stage 2 are distinct and independent of A; max(1-U, 1-U')^2 ~ U[0,1];
  Cauchy-Schwarz then gives the bound. That version stands without FPT.

**What FPT does not cover.** The LP relaxation, the gap LP* - T, duals, concentration (the paper
contains no variance statement), other value distributions, and K >= 4 pairwise models. This
agrees with GATE0_REPORT section 2(c). That report's theorem numbers (2, 5) are correct for arXiv
v1, and its consequence "E[3q - T] = Theta(q^{1/3})" is correct. The SYNTHESIS hedge "still needs
checking against the paper" can be discharged: E(3q - T) = E Z_q = Theta(q^{1/3}), with the upper
bound by FPT Thm 2(b).

## 2. The variance argument, step by step

**(a) Uniform representation.** For x in [0,1), P(1 - e^{-E} <= x) = P(E <= -ln(1 - x)) =
1 - e^{ln(1-x)} = x. So 1 - e^{-E} ~ U[0,1], and X = 1 - e^{-E} coordinatewise gives iid U[0,1]
pair values. The map E -> X is a C^infinity diffeomorphism of (0,inf)^n onto (0,1)^n, so it maps
null sets to null sets. **Holds.**

**(b) T.** f_T(E) = max over the finitely many perfect assignments N of
L_N(E) = sum_{i in 1_N} (1 - e^{-E_i}).
- *Pair count.* The 3q pairs of N are distinct, because two disjoint teams share no agent and hence
  no pair. The map N -> 1_N is injective, since the X-Y and X-Z pairs recover each team. So
  distinct assignments give distinct linear forms in X.
- *Differentiability.* Each L_N is C^1 with partial derivatives in [0,1], so f_T is globally
  Lipschitz (constant sqrt(3q)). Off the set where two distinct forms tie, which is a finite union
  of hyperplanes {<1_N - 1_N', X> = 0} with nonzero normal (null in X, hence in E), the maximiser
  M is unique and locally constant. So f_T = L_M near E, and
  dT/dE_i = 1{i in M} e^{-E_i} = 1{i in M}(1 - X_i).
- *Gradient bound.* \|grad\|^2 = sum_{i in M}(1 - X_i)^2 <= sum_{i in M}(1 - X_i) = 3q - T, using
  0 <= 1 - X_i <= 1 and T = sum_{i in M} X_i over exactly 3q pairs.

**Holds.**

**(c) LP\*.** LP* = max over w in P of sum_i X_i l_i(w), where l_i(w) = sum of w_f over teams f
containing pair i.
- *Loads are in [0,1].* The teams containing pair (a,b) are a subset of the teams containing agent
  a, whose total weight is 1, so l_i in [0,1]. Each team contains 3 pairs, so
  sum_i l_i(w) = 3 sum_f w_f = 3q. (The parenthetical "pair weight = sum of team weights
  containing the pair, at most 1" is right.)
- *Finite maximum.* P is a bounded polytope, so LP* = max over its finitely many vertices. That is
  a maximum of finitely many linear forms <l(w^v), X>, and f_LP is Lipschitz exactly as in (b).
- **Hole (needs a fix).** CONCENTRATION.md says "lambda = the (a.s. unique) LP optimum". In team
  weights this is false in general: w -> l(w) is not injective on P, and distinct vertices can
  share a load vector, in which case they tie for every X. Example, q = 2 (exact, by hand):
  - w1 = 1/2 on {111, 122, 212, 221} and w2 = 1/2 on {112, 121, 211, 222}.
  - Each is the unique solution of the agent constraints on its support, so each is a vertex of P.
  - Both have every one of the 12 pair loads equal to 1/2.
- **Fix.** Index the maximum by the finite set of distinct load vectors {l(w^v)}. Ties between
  distinct load vectors lie on finitely many null hyperplanes. Off them the optimal load vector
  l* is unique and locally constant, and dLP*/dE_i = l*_i (1 - X_i). Read lambda_i in SYNTHESIS
  as l*_i. Then
  \|grad\|^2 = sum l*_i^2 (1-X_i)^2 <= sum l*_i (1-X_i) = 3q - LP* <= 3q - T.
- *Knock-on.* A3 and B3 use team weights lambda_e and w*_e. Their bounds hold for any choice of
  optimal team weights with the given load vector, so they survive. Only the words "a.s. unique"
  should be changed.

**Holds with fix.**

**(d) Poincaré.** BL Lemma 2.1 (p. 387) proves Var_nu(f) <= 4 int f'^2 dnu for the two-sided
exponential nu_1. BL say that everything holds "with the obvious modifications" for the one-sided
exponential nu~_1, and they display that one-sided inequality as (1.7) (p. 386), attributing it to
Talagrand [Ta1] and Klaassen [Kl]. A two-line proof for Exp(1), which is BL's proof with sgn
dropped: for g = f - f(0), integration by parts gives E g^2 = 2 E[g g'], which is at most
2 (E g^2)^{1/2} (E g'^2)^{1/2}, so Var f <= E g^2 <= 4 E f'^2.
- *Conditions.* BL's class L_n asks for continuous, a.e. differentiable f with f and \|grad f\|
  integrable and e^{-x_i} f -> 0 at infinity. In fact the integration by parts needs f absolutely
  continuous on lines, since a.e. differentiability alone admits Cantor-type counterexamples. Our
  f_T and f_LP are bounded (in [0,3q]) and globally Lipschitz, so all of these hold.
- *Constant.* 4 is sharp in dimension 1 (BL remark p. 385-386).
- *Tensorisation.* Var_{mu^n} f <= sum_i E[Var_i f] (Efron-Stein; Boucheron-Lugosi-Massart,
  Concentration Inequalities (2013), Thm 3.1). Apply (1.7) to each coordinate section, which is
  Lipschitz and so absolutely continuous, with derivative equal a.e. (Fubini) to the a.e. partial
  derivative. This gives Var f <= 4 E \|grad f\|^2. The CONCENTRATION.md attribution
  "Bobkov-Ledoux 1997, PTRF 107; constant 4 sharp in dimension 1, tensorizes" is accurate. The
  precise pointer is BL (1.7) / Lemma 2.1 plus Efron-Stein.

**Holds.**

**(e)** Var T <= 4 E(3q - T), and Var LP* <= 4 E(3q - LP*) <= 4 E(3q - T). **Holds.**

**(f) Chebyshev / Proposition R.**
- *Core criterion.* D = LP* - T >= 0 and {core nonempty} = {D = 0}. This LP criterion is taken from
  PROOF_NOTES Thm 1 and THEOREM_VANISHING Lemma 2 and was not re-checked here.
- *Chebyshev.* If E D > 0, then P(D <= 0) <= P(\|D - ED\| >= ED) <= Var D / (ED)^2.
- *Variance of D.* Var D <= 2 Var T + 2 Var LP* <= 16 E(3q - T) holds. Poincaré applied to D
  directly gives 8: \|grad D\|^2 = sum (l*_i - 1_M(i))^2 (1-X_i)^2, and (l - m)^2 <= l + m on
  [0,1], so this is at most 2(3q - T). That is A3's first line.
- *Thresholds.* P(core nonempty) <= 8 E(3q - T) / (E D)^2.
  - With A2's O(sqrt q) bound, E D / q^{1/4} -> infinity suffices.
  - With FPT Thm 2(b), E D / q^{1/6} -> infinity suffices.

**Holds as stated.**

**Not holes, but worth recording.**
1. Proposition R's hypothesis is far from the data. E D = 0.156 at q = 48, while q^{1/6} = 1.9.
   D <= S - T = O(q^{1/3}), so the hypothesis requires E D between q^{1/6} and q^{1/3}. Only the
   "two constants" conjecture E D >= c q^{1/3} would deliver it.
2. CONCENTRATION.md section 4 item 1 says "E D -> infinity is probably false (median of D flat)".
   That is in tension with SYNTHESIS, where E D is still rising as D/q^{1/3} at q <= 48. Median
   and mean can differ, but the two notes should be reconciled before Proposition R is presented
   as the route.
3. The rigorous variance bound (~4 x 2.55 q^{1/3}, about 37 at q = 48) is about 40-150x the
   measured Var T (0.15-0.44 at q = 6..24). Proposition R is therefore far from sharp. This is not
   an error.

## Bottom line

- **Check 1: holds.** FPT covers our model exactly after U -> 1 - U, and E(3q - T) = E Z_q =
  Theta(q^{1/3}). The upper bound is FPT Thm 2(b), via greedy, constant ~ 3 x 1.623 asymptotically.
  The lower bound is Thm 2(a) / Lemma 3 with constant 1.623; the printed Gamma(3/4) is a typo for
  Gamma(4/3). An elementary bound of 1.363 q^{1/3} holds for all q >= 3.
- **Check 2: holds with one fix.** Replace "a.s. unique LP optimum" by "a.s. unique optimal
  pair-load vector". Optimal team weights can be non-unique, as the q = 2 example shows. With that
  fix, Proposition R is correct, and its constant can be improved from 16 to 8.
