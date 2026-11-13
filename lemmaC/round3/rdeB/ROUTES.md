# Step II-a by a second route: monotonicity in the dispersive order

Round 3, line `rdeB`. 60-minute time box. Nothing here has been read by another agent or model
family. Labels: **complete** = I would defend it to a referee; **gap** = not proved. Floating point
is evidence only.

## 0. Result

| Item | Status |
|---|---|
| Theorem M: T is monotone for the dispersive order on IFR laws (§2) | **complete** (pencil proof; 60 random pairs and the iterates agree) |
| Corollary A: every fixed point has cost c ≤ c_n := 3·Γ(T^n δ_0), for every n (§3) | **complete** |
| Corollary B: a fixed point exists (monotone limit below an explicit super-solution) (§4) | **complete**, proof written in one pass, not re-read by anyone |
| c_5 < 2.30 (§5) | **gap**: floating point gives c_5 = 2.2678, c_6 = 2.2558. Needs one validated evaluation of an explicit 5-fold iterated integral; one-sided rounding is free by Theorem M |
| Routes (a) inequalities only, (b) variational, (c) positive temperature, capped recursion (§6) | assessed; none closes alone |

So the step no longer needs Newton–Kantorovich, an approximate inverse, or operator norms. What is
left is a one-sided numerical bound on an explicit number with margin 0.032 (n = 5) or 0.044 (n = 6).

## 1. Setting

E, E_1, E_2 are iid Exp(1). A law with survival exp(−G), G convex nondecreasing (an **IFR** law), is
the law of Q(E) with Q = G^(−1) : (0,∞) → R **concave nondecreasing**. Let 𝒬 be the set of such Q for which Q(E) has moments of every order.

- Order: Q_0 ⪯ Q_1 iff Q_1 − Q_0 is nondecreasing (the dispersive order). It ignores additive
  constants. 0 ⪯ Q for every Q (0 is δ_0). 𝒬 is convex and Q_0 ⪯ Q_1 implies Q_0 ⪯ Q_τ ⪯ Q_1 for
  Q_τ = Q_0 + τ(Q_1 − Q_0).
- The map: g_Q(t) = E(t + Q(E_1) + Q(E_2))_+^3 / 6 and T(Q) = g_Q^(−1) on (0,∞). This is RDE.md's T
  written on quantile functions. g_Q is convex, so **T(Q) ∈ 𝒬 for every input** (moments: RDE.md 1.2).
- Shift law T(Q + a) = T(Q) − 2a. Scaling law (new, exact): **T(aQ)(θ) = a·T(Q)(θ/a³)**, because
  g_{aQ}(t) = a³·g_Q(t/a).
- Slope at level: s_Q(θ) = g_Q'(g_Q^(−1)(θ)), so T(Q)'(θ) = 1/s_Q(θ). Hence
  T(Q_0) ⪯ T(Q_1) iff s_{Q_1} ≤ s_{Q_0} on (0,∞).
- Γ(Q) = E[T(Q)(E)] + 2·E[Q(E)], invariant under Q → Q + a. At a fixed point Q* = T(Q*):
  Γ(Q*) = 3·E X, so the cost of RDE.md §1.5 is **c = 9·E X = 3·Γ(Q*)**.

## 2. Theorem M (complete)

**Lemma 1.** For X = Q(E), Q ∈ 𝒬, put M_k(u) = E(u + X)_+^k / k!. Then M_1² ≥ M_0·M_2, so M_2/M_1
is nondecreasing on {M_1 > 0}.
Proof. M_0(u) = P(X > −u) = exp(−G(−u)) is log-concave because G is convex. A primitive
∫_{−∞}^u of a log-concave function is log-concave (Prékopa), so M_1 and M_2 = ∫M_1 are
log-concave; log-concavity of M_2 with M_2' = M_1, M_2'' = M_0 a.e. is M_1² ≥ M_0 M_2, and
(M_2/M_1)' = (M_1² − M_0M_2)/M_1². ∎

**Theorem M.** If Q_0, Q_1 ∈ 𝒬 and Q_0 ⪯ Q_1 then T(Q_0) ⪯ T(Q_1).
Proof. Δ = Q_1 − Q_0 is nondecreasing, Q_τ = Q_0 + τΔ ∈ 𝒬. Write X^i = Q_τ(E_i), Y = X^1 + X^2,
g = g_{Q_τ}. By dominated convergence and symmetry in (E_1, E_2),

  ∂_τ g(t) = E[(t+Y)_+² Δ(E_1)],  ∂_τ g'(t) = 2·E[(t+Y)_+ Δ(E_1)].

Fix θ > 0 and t = g^(−1)(θ) (g'(t) > 0). The implicit function theorem gives ∂_τ t = −∂_τ g/g', so

  d/dτ s_{Q_τ}(θ) = ∂_τ g'(t) − g''(t)·∂_τ g(t)/g'(t)
                 = 2·E(t+Y)_+ · ( E^(1)[Δ(E_1)] − E^(2)[Δ(E_1)] ),

where E^(k) is the expectation under the probability proportional to (t+Y)_+^k dP. Under P^(k) the
law of E_1 has density proportional to e^(−e)·M_k(t + Q_τ(e)), with M_k as in Lemma 1 for the law of
X^2. The ratio of the P^(2)-density to the P^(1)-density is proportional to
(M_2/M_1)(t + Q_τ(e)), nondecreasing in e by Lemma 1. So the P^(2)-law of E_1 dominates the
P^(1)-law (likelihood ratio, hence stochastically), and Δ is nondecreasing, so the bracket is ≤ 0.
Hence s_{Q_1}(θ) ≤ s_{Q_0}(θ) for every θ. ∎

Remarks. (i) This is the monotonicity RDE.md §3 item 3 guessed and dropped ("convolution preserves
the dispersive order only against log-concave laws"). No such preservation is needed: only the
slope-at-level function of the output has to be compared, and IFR is exactly what makes M_2
log-concave. (ii) IFR costs nothing: every output of T is IFR. (iii) It agrees with RDE.md's
spectrum: the leading non-shift multiplier of DT is −2·(−0.175) = +0.35 > 0.

Checks (`disp_thmM.py`, `disp_check.py`): 60 random IFR pairs Q_0 ⪯ Q_0 + λQ_b, 200 levels θ in
[1e-6, 20]: largest (s_1 − s_0)/s_0 = −1.1e-3 (never positive). Along T^n δ_0, n = 1..9, θ in
[1e-9, 40]: s decreases at every level at every step.

## 3. Corollary A: enclosure for every fixed point (complete)

**Lemma 3.** Q_0 ⪯ Q_1 in 𝒬 implies Γ(Q_0) ≥ Γ(Q_1).
Proof. Centre both (Γ is shift-invariant). Then Q_1(E) = Q_0(E) + Δ(E), Δ nondecreasing, E Δ(E) = 0.
For convex φ and a constant a, E φ(Q_1(E) + a) ≥ E φ(Q_0(E) + a) + Cov(φ'(Q_0(E) + a), Δ(E)) and
the covariance of two nondecreasing functions of E is ≥ 0. Apply this twice (once in each of
E_1, E_2) with φ(y) = (t + y)_+³/6: g_{Q_1} ≥ g_{Q_0} pointwise, so T(Q_1) ≤ T(Q_0) pointwise and
Γ(Q_1) ≤ Γ(Q_0). ∎

**Corollary A.** Let Q_n = T^n(0) (so Q_1(θ) = (6θ)^(1/3), the law with survival exp(−x³/6)) and
c_n = 3·Γ(Q_n). Then Q_0 ⪯ Q_1 ⪯ Q_2 ⪯ …, c_0 ≥ c_1 ≥ …, and **every proper solution of
G(x) = E(x + X' + X'')_+³/6 has c = 9·E X ≤ c_n for every n.**
Proof. A solution is an output of T, hence Q* ∈ 𝒬 with T(Q*) = Q*. 0 ⪯ Q*, so by Theorem M and
induction Q_n ⪯ Q*; Lemma 3 gives 3Γ(Q*) ≤ 3Γ(Q_n). Monotonicity of Q_n: 0 ⪯ Q_1, apply T. ∎

Floating-point values (`disp_iter.py`, grid h = 0.01 on [−16, 16]):

| n | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| c_n | 4.86795 | 3.41521 | 2.66013 | 2.39798 | 2.30163 | 2.26777 | 2.25584 | 2.25169 | 2.25027 | 2.24963 |

c_0 = 3·6^(1/3)·Γ(4/3) is Frieze–Pegden–Tkocz's greedy constant 4.87, as it should be (one step of
the recursion from "partners have price 0"). The ratio of successive differences tends to 0.35.
c_4 is above 2.30; **n ≥ 5 is needed.**

## 4. Corollary B: existence (complete; written once, please re-read)

**Lemma 4 (super-solution).** Let R(θ) = ln θ for θ ≤ 1 and 3(θ^(1/3) − 1) for θ ≥ 1 (concave,
increasing; the law has P(X ≤ x) ≈ e^x below and survival exp(−(1 + x/3)³) above). There is a < ∞
with T(aR) ⪯ aR.
Proof. Write ℓ_Q(θ) = θ·Q'(θ). By the scaling law, T(aR) ⪯ aR iff ℓ_{TR}(ϑ) ≤ ℓ_R(a³ϑ) for all
ϑ > 0, and ℓ_R(θ) = max(1, θ^(1/3)). With Y = R(E_1) + R(E_2) and M_k for Y,
ℓ_{TR}(ϑ) = g/g' at the point t where g(t) = ϑ, and g/g' = M_3/M_2 ≤ M_2/M_1 ≤ M_1/M_0: each M_k is
log-concave in u as in Lemma 1, applied to Y, which is IFR because a sum of independent IFR
variables is IFR (Barlow–Marshall–Proschan 1963, cited from memory, not opened). M_1/M_0 is the
mean residual life of Y above −t, at most 1/(hazard rate of Y at −t).
(i) ϑ → 0 means t → −∞. Y is IFR and has all exponential moments (cubic upper tail), so its hazard
rate tends to infinity and its mean residual life tends to 0. Hence ℓ_{TR}(ϑ) → 0; pick ϑ_1 with
ℓ_{TR} ≤ 1 on (0, ϑ_1].
(ii) ϑ → ∞: for t ≥ t_0 with P(Y > −t_0) ≥ 1/2, g/g' ≤ M_1/M_0 ≤ 2(t + E|Y|), and
ϑ = g(t) ≥ (t + E Y)_+³/6 (Jensen), so ℓ_{TR}(ϑ) ≤ C'(1 + ϑ^(1/3)). With continuity, A := sup_{ϑ ≥ ϑ_1} ℓ_{TR}(ϑ)/ϑ^(1/3) < ∞.
Take a = max(1, A). For ϑ ≤ ϑ_1: ℓ_{TR}(ϑ) ≤ 1 ≤ ℓ_R(a³ϑ). For ϑ ≥ ϑ_1:
ℓ_{TR}(ϑ) ≤ a·ϑ^(1/3) = (a³ϑ)^(1/3) ≤ ℓ_R(a³ϑ). ∎

**Corollary B.** T has a fixed point, namely the least one: Q_n = T^n(0) increases to Q_∞ with
T(Q_∞) = Q_∞ + k, and Q_∞ + k/3 is a fixed point of T with cost lim c_n ≤ c_n for all n.
Proof. 0 ⪯ aR and T(aR) ⪯ aR give Q_n ⪯ aR for all n (Theorem M). The slopes Q_n'(θ) increase in
n and are bounded by aR'(θ); let Q_∞' be the limit (nonincreasing in θ), and normalise
Q_n(1) = Q_∞(1) = 0. Then |Q_n(E)| ≤ a|R(E)|, which has all moments, so by dominated convergence
g_{Q_n} → g_{Q_∞} and g'_{Q_n} → g'_{Q_∞} locally uniformly, hence s_{Q_n}(θ) → s_{Q_∞}(θ) > 0 and
T(Q_∞)' = lim Q_{n+1}' = Q_∞'. So T(Q_∞) = Q_∞ + k and, by the shift law,
T(Q_∞ + k/3) = Q_∞ + k/3 (RDE.md's quotient lemma). ∎

With RDE.md §1.3–1.5 (tree process, matching, cost identities, all conditional only on a fixed
point) this gives: **an invariant random integral perfect matching exists on the Poisson hypertree,
with cost ≤ c_n for every n.** Numerically a decreasing iteration from several super-solutions
(`disp_tests.py`: G^0(x) = A·E(x + σN)_+³, min s_out/s_in between 1.09 and 1.44) converges to the
same cost 2.2495 from below (2.2486 after 8 steps), so the fixed point is very likely unique; I do
not need or claim uniqueness.

## 5. The remaining gap: a certified upper bound on c_5 (or c_6)

c_5 is an explicit number: five applications of an explicit integral transform to δ_0, then one mean.
Theorem M and Lemma 3 make every rounding safe as long as it goes one way:

- replace an iterate Q_n by any Q̃ ∈ 𝒬 with Q̃ ⪯ Q_n; then T(Q̃) ⪯ Q_{n+1} and the final 3Γ only
  goes up;
- a valid Q̃: on a mesh θ_0 < … < θ_N take slope 1/ŝ_j on (θ_{j−1}, θ_j] with ŝ_j ≥ s(θ_j) a
  validated upper bound (s is nondecreasing, so this is below the true slope and stays concave),
  slope 1/ŝ_0 below θ_0, slope 0 above θ_N (an atom at the top). Tails need no estimate at all.
- Q̃ is then a finite piecewise-exponential law; g, g' for the sum of two such variables are finite
  sums of exponential-polynomials, so ŝ_j can be computed in interval arithmetic.

Floating-point trial of exactly this scheme (`disp_rounded.py`; crude first-order rounding,
log-spaced cells, rounding applied at every stage):

| mesh | stage 5 | stage 6 | stage 8 |
|---|---|---|---|
| 400 cells on [1e-6, 30] | 2.3827 | 2.3724 | 2.3677 |
| 800 cells on [1e-4, 20] | 2.3094 | 2.2981 | 2.2928 |
| 1600 cells on [1e-4, 20] | 2.2887 | 2.2771 | 2.2717 |
| 3200 cells on [1e-4, 20] | 2.2783 | 2.2666 | 2.2610 |

The loss is O(1/N). With 1600 cells the rounded bound is below 2.30 from stage 5 on, with margin
0.02 to 0.03 left for the interval-arithmetic slack. This is a finite validated computation of
modest size (six stages, each a table of a few thousand values of an explicit function); it was not
carried out. **Until it is, "alpha < 2.30" is not proved.** The trial itself is floating point on an
x-grid and smears the top atom, so its numbers are indicative only.

## 6. The brief's other routes

- (a) Inequalities only. Jensen in g gives c ≤ 3·E W = 4.868 for every fixed point (this is c_0), and
  Lyapunov gives c ≥ 1.36. Corollary A is the systematic version: c ≤ c_n. Without a computation the
  best pencil value is c_1 = 3.415 (a single explicit double integral), not 2.30.
- (b) Variational. Not pursued. The monotone structure makes it unnecessary.
- (c) Positive temperature. The shift law is exact at every β: exp(−βX) = Σ_k exp(−β(ξ_k − X' − X''))
  sends a shift c to −2c. So no contraction in the unbounded message, at any temperature, and the
  energy would still need a certified solve. Not better than the zero-temperature route.
- Capped recursion (new). X = min(κ, min_k(ξ_k − X' − X'')) maps laws on [−2κ, κ] into themselves,
  so Schauder–Tychonoff gives a fixed point for every κ with no shift problem; it is the
  maximum-weight matching with reward 3κ per team, leaving a fraction exp(−E(κ + Y)_+³/6) of agents
  unmatched (about 1e-14 at κ = 5). A κ → ∞ limit gives an invariant perfect matching, but with cost
  equal to a liminf of unlocated numbers; the cap breaks the shift-invariance Theorem M lives on.
  Superseded by Corollary B. Gap, not needed.

## 7. Honest list of gaps and risks

1. c_5 < 2.30 certified (§5). Not done.
2. Corollaries A and B rest on RDE.md §1.2–1.5 (moments, tree process, c = 9·E X), which I read but
   did not re-derive line by line.
3. Theorem M and Lemma 4 were proved and written in under an hour and have had no second reader.
   Points a referee should check: differentiation under the expectation along Q_τ (needs third
   moments; all moments are available); Lemma 1 when the law has an atom at its upper end
   (M_0 has a jump; the a.e. form of the inequality is what is used); in Lemma 4, "hazard rate of Y
   tends to infinity" (IFR plus finite exponential moments of every order).
4. The matching built from the fixed point uses extra randomness (it is not claimed to be a factor
   of the tree); LP_LOCALITY's hypothesis (H0) allows that.

5. No research-handoff or ledger entry was made from this line (writes were confined to `rdeB/`).

## 8. Files

`disp_iter.py` (iterates and c_n), `disp_check.py` (slope monotonicity along the iteration),
`disp_thmM.py` (random IFR pairs and a non-IFR control), `disp_tests.py` (super-solutions and the
decreasing iteration), `disp_rounded.py` (one-sided rounding trial). Each runs in under a minute with
`/Users/maxf/miniforge3/bin/python3`.

## 9. Confidence

Theorem M 90%. Corollary A given M 95%. Corollary B (existence) 85%. A certified c_5 or c_6 below
2.30 by §5: 90% that it works, perhaps one to two days of interval arithmetic, much less delicate
than Newton–Kantorovich because nothing has to be inverted and every error has a safe direction.
