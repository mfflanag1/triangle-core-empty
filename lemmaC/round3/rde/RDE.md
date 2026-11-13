# Step (II-a): the tree recursion X = min_k(xi_k − X_k' − X_k'')

Round 2, line `rde`. 55-minute time box. Floating point is evidence only. Nothing here has
been read by another model family. Labels: **complete** = I would defend it to a referee;
**gap** = not proved.

## 0. Result in one table

| Item | Status |
|---|---|
| Equation for G, consistency of the matching, cost formulas c = E S_+^4/8 = 9·E X, normalisation E S_+^3 = 6 (§1) | **complete**, given a fixed point |
| Shift law, quotient lemma, universal a priori bounds, linearisation DT = −2P with P a reversible Markov kernel (§2) | **complete** |
| Existence of a fixed point (§3) | **gap**. Reduced to one checkable invariance statement; not verified |
| Cost of the numerical fixed point (§4) | c = 2.2495456100, E X = 0.2499495122 (three formulas agree to 1e-10; fourth-order grid convergence observed). Not an enclosure |
| Certified enclosure (§5) | **gap**. Plan and constants given; not carried out |
| Simpler invariant matching below 2.30 (§6) | none found |

The constant is not 9/4: c = 2.24954561, E X = 0.24994951.

## 1. Derivations (complete, conditional on a fixed point)

Notation: X, X', X'', X_1, X_2, X_3 iid with law mu; Y = X' + X''; S = X_1 + X_2 + X_3;
Fbar(x) = P(X > x); G = −ln Fbar.

**1.1 Equation for G.** The points (xi_k, Y_k) form a Poisson process on [0,∞) × R with
intensity (xi^2/2) dxi × law(Y). The number of points with xi_k − Y_k ≤ x is Poisson with mean
∫ (xi^2/2) P(Y ≥ xi − x) dxi = E[(x + Y)_+^3]/6. So X = min_k(xi_k − Y_k) has

  G(x) = E[(x + X' + X'')_+^3] / 6.   (1)

A proper law solves (1) only if E Y_+^3 < ∞; then G is finite everywhere, C^3, with
G' = E(x+Y)_+^2/2, G'' = E(x+Y)_+, G''' = P(Y > −x), and G'''' = density of −Y. The minimum is
attained and a.s. unique (finitely many points below any level; continuous law).

**1.2 Moments are automatic.** For a solution: G(x) ≥ (x − m)^3·P(Y > −m)/6 for x > m, so X_+
has all moments; and P(X ≤ −t) ≤ G(−t) = E(Y − t)_+^3/6, so E X_− ≤ E Y_+^4/24 < ∞. Hence
E|X| < ∞ and E S_+^4 < ∞.

**1.3 Recursive tree process.** Fix a solution mu. For each n put iid mu-messages on the
agents at depth n (directed to their parent teams) and compute all messages
X_{v→e} = min over teams k ∋ v, k ≠ e, of (xi_k − X_{b_k→k} − X_{c_k→k}) inside the n-ball, in
both directions (outward messages are functions of inward ones, by induction on depth). By
(1) and independence of subtrees the laws for different n are consistent, so Kolmogorov
gives a process on the infinite tree in which the messages entering any finite connected
sub-hypertree from outside are iid mu. Unimodularity: under the measure that counts a team e
at the root (Mecke), e has cost density xi^2/2 and its three agents carry iid copies of the
same planted subtree-with-messages (the root's message to e is the minimum over its other
teams, which is again mu and independent of the rest). This description is symmetric in the
three agents, which is involution invariance.

**1.4 The matching.** Define e = (a,b,c) ∈ M iff xi_e < X_{a→e} + X_{b→e} + X_{c→e}.
The condition is symmetric in a, b, c. At an agent v write Z_k = xi_k − X_{b_k→k} − X_{c_k→k}
for its teams. Then X_{v→e} = min_{k≠e} Z_k, so e ∈ M iff Z_e < min_{k≠e} Z_k, i.e. iff e is
the argmin of Z at v. The argmin exists and is unique a.s., so every agent is in exactly one
team of M. M is an invariant random perfect matching. (It is a function of tree and
messages; whether it is a function of the tree alone is not claimed. Given the multiplier
−2 below I expect it is not.)

**1.5 Cost and normalisation.** By Mecke at the root, for f ≥ 0,
E Σ_{e∋o} f(xi_e, X_{o→e}, X_{b→e}, X_{c→e}) = ∫ (xi^2/2) E f(xi, X_1, X_2, X_3) dxi. Hence
- f = 1{xi < S}: 1 = E #{e ∋ o, e ∈ M} = **E S_+^3 / 6**.
- f = xi·1{xi < S}: c := E[cost of the root's team] = **E S_+^4 / 8**.
- f = (X_2 + X_3)·1{xi < S} (integrable by 1.2): E[sum of the two partner messages into the
  root's team] = E[(X_2+X_3)·S_+^3]/6 = (2/3)·E[S·S_+^3]/6 = E S_+^4 / 9.
- The root's own minimum is X_o = xi_{e*} − (two partner messages), and X_o has law mu. So
  E X = E S_+^4/8 − E S_+^4/9 = E S_+^4/72, i.e. **c = 9·E X**.
The cavity note's "Bethe form" 3·E X + E S_+^4/12 equals E S_+^4·(1/24 + 1/12) = E S_+^4/8, so
all three formulas are the same identity. Also E[message of an agent into its own team] =
E S_+^4/18 = 4c/9.

**1.6 Normalisation of the market.** c is the expected rescaled cost xi = q^(2/3)·deficit of
the team containing a given agent. A perfect assignment has q teams and each of the q agents
of sector 1 is in exactly one, so total deficit = q·c·q^(−2/3) = **c·q^(1/3)**. (Per agent
the cost is c/3; tree.md's "per team" convention is this c.)

## 2. Structure of the map T: mu ↦ law with G_out(x) = E(x+Y)_+^3/6 (complete)

**2.1 Antitone, shift law.** If mu ≤_st nu then T(mu) ≥_st T(nu). T(mu∗δ_c) = T(mu)∗δ_{−2c}
exactly (replace x + Y by x + Y + 2c).

**2.2 Quotient lemma.** If T(mu) = mu∗δ_c then mu∗δ_{c/3} is a fixed point:
T(mu∗δ_{c/3}) = mu∗δ_{c − 2c/3}. So T descends to a map Tbar on laws modulo shift, and fixed
points of Tbar are exactly shift classes containing one fixed point of T.

**2.3 Universal bounds on any output.** For every input: (a) the density G'e^{−G} is at most
sup_g (6g)^{2/3} e^{−g}/2 = 0.647 (Lyapunov: E u_+^2 ≤ (E u_+^3)^{2/3}); (b) if κ is the point
with G(κ) = θ, then for x ≥ 0, θ + x·G'(κ) + x^2 G''(κ)/2 + x^3 G'''(κ)/6 ≤ G(κ + x) ≤
(x + (6θ)^{1/3})^3/6 (lower: G'''' ≥ 0; upper: Minkowski in L^3). (c) If the input has
P(X > 0) = p then κ ≤ (6θ/p^2)^{1/3}.

**2.4 Linearisation.** Move each atom of mu by εφ(X), φ bounded. Then
δG(x) = ε·E[(x+Y)_+^2·φ(X')] = 2ε·G'(x)·(Pφ)(x), and the output atoms move by εψ with

  ψ = −2·Pφ,   P(x, da) = f(a)·g_1(x + a) da / G'(x),   g_1(t) = E(t + X)_+^2 / 2.

P(x,·) is the conditional law of a partner message of the argmin team given X = x. It is a
Markov kernel (P1 = 1, which is the shift law) and, because g_1(x + a) is symmetric, it is
reversible with respect to π ∝ f·G'. So DT = −2P is self-adjoint on L^2(π) with real
spectrum in [−2, 2]; eigenvalue −2 is the shift.

**Numerical spectrum of P at the fixed point** (`rde_spectrum.py`, grids h = 0.04 and 0.02
agree to 6 digits; finite-difference check of ψ = −2Pφ on three test fields agrees to 1e-5):

  1, −0.1750, 0.0417, −0.0123, 0.0042, −0.0016, 0.0006, …

So on shapes (modulo shift) Tbar is a local contraction with factor 2·0.175 = **0.350**, and
‖(I − DT)^{−1}‖ on L^2(π) is 1/min|1 + 2p| = 1/0.650 = **1.54**. The fixed point is isolated
and non-degenerate; the only expanding direction is the shift. The observed rate of the
iteration f ← (2f + Tf)/3 is 0.783 = (2/3)(1 + 0.175), as predicted.

## 3. Existence: gap

**What is proved.** By 2.2 it is enough to find a fixed point of R = (re-centre)∘T on a
slice. Take the slice Σ_p = {mu: mu(0,∞) = p} (convex under mixtures) and re-centre the
output at the point κ with G_out(κ) = −ln p.

**Proposition E (complete).** Let K ⊂ Σ_p be nonempty, convex, weakly closed, tight, with
X_+^3 uniformly integrable over K and densities bounded by 0.647. If R(K) ⊂ K then T has a
fixed point mu* = nu∗δ_{κ(T nu)/3} for some nu ∈ K, the process of §1 exists, and
c = 9·E X* ≤ 9·sup_{nu ∈ K}(E_nu X + κ(T nu)/3).
Proof: on K, T is weakly continuous (uniform integrability), κ is continuous (G_out strictly
increasing where positive), mu(0,∞) is weakly continuous (bounded densities), so R is a
continuous self-map of a compact convex subset of a locally convex space; Schauder–Tychonoff;
then 2.2 and 1.5. ∎

**What is not proved: an explicit K with R(K) ⊂ K.** What I tried:
1. Order intervals [L, U] in stochastic order. **Dead (complete negative lemma):** if
   T[L,U] ⊂ [L,U] and the interval contains shifts mu*∗δ_c exactly for c ∈ [c−, c+], then
   −2c+ ≥ c− and −2c− ≤ c+, which forces c− = c+ = 0. No interval with room around the fixed
   point is invariant; pinned envelopes through (0, p) fail the same way, because order
   bounds do not pinch at the pin. T^2 is monotone but multiplies shifts by 4, so
   Tarski/Knaster does not apply either; shift-compensated variants of T or T^2 are not
   monotone.
2. A priori envelopes from 2.3 (upper tail exp(−θ − a·x^3/6), or exp(−θ − b·x)). The output
   constant is far below the input constant in every version I evaluated by hand (rough
   example: θ = 1, a = 0.5 returns only a ≥ 0.012). The inequalities lose the correlation between the
   re-centring and the tail. Not closed; I do not claim it cannot be closed.
3. Dispersive-order intervals (bounds on the quantile density, which are shift-free and convex
   in quantile coordinates). Tbar looks monotone in this order heuristically, but the
   convolution step preserves the dispersive order only against log-concave laws, and
   G'e^{−G} is not log-concave for general inputs. Not pursued.

**Route I recommend.** §2.4 says Tbar is a contraction with factor about 0.35 near the
solution and I − DT has inverse norm about 1.5. That is the textbook setting for a
computer-assisted Newton–Kantorovich proof (§5), which gives existence, local uniqueness and
the cost enclosure together. A pencil proof of existence alone would not give the cost.

## 4. Numerical fixed point (evidence)

`rde_solve.py`: density on a uniform grid on [−14, 14]; Y by FFT convolution (trapezoid on
the line, spectrally accurate); G, G' by exact partial moment sums (error O(h^4) from the
third-derivative kink); iteration f ← (2f + Tf)/3 to L1 residual below 1e-13.

| h | E X | 9·E X | E S_+^4/8 | E S_+^3/6 | sd X |
|---|---|---|---|---|---|
| 0.04 | 0.249949519622 | 2.2495456766 | 2.2495456103 | 1.0000000010 | 0.863244141741 |
| 0.02 | 0.249949512687 | 2.2495456142 | 2.2495456101 | 1.0000000001 | 0.863244152949 |
| 0.01 | 0.249949512254 | 2.2495456103 | 2.2495456100 | 1.0000000000 | 0.863244153649 |
| 0.005 | 0.249949512226 | 2.2495456100 | 2.2495456100 | 1.0000000000 | 0.863244153693 |

Observed order 4.00; Richardson value **c = 2.2495456100(1)**, E X = 0.2499495122.
Cross-checks: cavity note 2.24954 (E X = 0.24995, sd 0.8632); red-team script 2.2496. Three
independent codes, two different discretisations, three cost formulas.

Heuristic a posteriori reading: with ‖(I − DT)^{−1}‖ ≈ 1.5 and discretisation defect about
1e-9, the error of E X is about 1e-9. Reaching 2.30 needs an error of 5.6e-3 in E X, six
orders of magnitude more.

## 5. What a rigorous enclosure needs (gap; not started)

1. Represent F_0 as an explicit piecewise polynomial on [−L, L] (rational coefficients) with
   explicit tails. Then G_0(x) = E(x + Y_0)_+^3/6 is piecewise polynomial and computable
   exactly; the defect ‖exp(−G_0) − Fbar_0‖ needs interval arithmetic only for exp.
2. Newton-type operator N = I − A(T − I) with A a finite-rank approximate inverse of DT − I
   (Galerkin block plus identity on the tail). DT is an integral operator with the explicit
   smooth kernel of 2.4, so the projection error is O(h^2) with explicit constants.
3. Radii-polynomial bounds in a weighted L^1 or sup norm on Fbar: Y_0 (defect), Z_0, Z_1
   (linear), Z_2 (second derivative of a quadratic-then-exponential map, explicit).
4. Tails: upper by 2.3(b); lower by P(X ≤ −t) ≤ E(Y − t)_+^3/6.
5. Cost: |E X* − E X_0| ≤ ‖Fbar* − Fbar_0‖_{L^1} ≤ r. The margin allows r = 5.6e-3, so a
   coarse discretisation (a few hundred cells) should be enough.
Estimated effort: two to four days with interval arithmetic (mpmath `iv` or Arb), mostly the
tail and Z_2 estimates. None of it was carried out in this time box.

## 6. Simpler invariant matchings below 2.30 (task 4): none

- Greedy or mutually-cheapest rules: Frieze–Pegden–Tkocz's greedy constant is 4.87; no local
  greedy rule is within 2.3% of the optimum.
- Finite-depth dynamic programming with a boundary law that is only approximately fixed: the
  messages at depth j from the boundary have law T^j(nu), which leaves the solution at rate
  2^j, so the three agents of a team disagree with probability not tending to 0. A packing
  with slack → 0 needs the exact fixed point (or boundary accuracy 2^(−r), which is existence
  again).
- Entropic (Sinkhorn-type, positive temperature) scalings x_e = u_a u_b u_c w(xi_e): the
  message recursion has the same shift multiplier −2 and the same existence problem, and a
  cost at or above the zero-temperature value.
- Local weak limits of finite LP optima are invariant fractional matchings, but their cost is
  the unknown we are bounding.
- No closed form: the bipartite analogue is logistic; here E X = 0.249949512, not 1/4.
The fractional fixed point (2.197) has twice the margin but a function-valued message; its
existence is harder, not easier.

## 7. Files

- `rde_solve.py` (solver; writes `fixed_point_h005.npy`, 180 KB: x, f, G, G').
- `rde_spectrum.py` (spectrum of P), `rde_fdcheck.py` (check of DT = −2P and iteration rate).
Run: `/Users/maxf/miniforge3/bin/python3 rde_solve.py` (6 s).

## 8. Confidence

- §1 and §2.1–2.4 as mathematics: 95%.
- A fixed point exists with c = 2.24954561: 99% (as a fact), unproved.
- A certified alpha < 2.30 is attainable by §5: about 85–90%; the risk is effort in the tail
  estimates, not the margin.
- A pencil-only existence proof with cost bound: below 30%.
