# Certificate: c* < 2.30 (rigorous upper bound on the cost of every fixed point)

Round 3, line `cert`. Written and run by one agent in one sitting; no second reader yet.
Labels: **complete** = I would defend it to a referee; **assumed** = taken from elsewhere.

## 0. Result

**Certified: every fixed point Q* = T(Q*) has c* = 9·E X ≤ 2.251887 < 2.30** (run B, stage 10).
A coarser run of under a minute (run A) already gives c* ≤ 2.257381 at stage 8 and 2.274755 at stage 5.

| stage n | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| certified, run A (h = 2^-10, 4000 levels) | 4.867956 | 3.417530 | 2.665068 | 2.404165 | 2.308352 | 2.274755 | 2.262954 | 2.258825 | 2.257381 | |
| certified, run B (h = 2^-12, 12000 levels) | 4.867955 | 3.415996 | 2.661760 | 2.400021 | 2.303815 | 2.270053 | 2.258185 | 2.254029 | 2.252574 | 2.251887 |
| reference c_n without rounding (referee's code, floating) | 4.86795 | 3.41528 | 2.66022 | 2.39809 | 2.30170 | 2.26786 | 2.25596 | | | limit 2.24954561 |
| floating c_n in ROUTES.md (low by about 1e-4 according to the referee) | 4.86795 | 3.41521 | 2.66013 | 2.39798 | 2.30163 | 2.26777 | 2.25584 | 2.25169 | 2.25027 | 2.24963 |

Each certified entry is an upper bound on 3·Γ(Q̃_n) for an explicit law Q̃_n, hence on c*. Every
entry from stage 5 on is below 2.30 in both runs.

What the certificate rests on:

1. **assumed**: Theorem M and Lemma 3 of `rdeB/ROUTES.md` (T is monotone for the dispersive order on
   𝒬; Γ is antitone), applied to laws in 𝒬 that have an atom at their upper end. Not re-proved here.
   The round's referee reports that Theorem M, Corollary A and Corollary B hold (no error found).
2. **complete**: sections 1 to 4 below (the scheme, its safety, the enclosures).
3. Trust in the four basic operations and the comparisons of `mpmath.iv` (outward-rounded interval
   arithmetic; its `exp` is **not** used, see section 5) and in IEEE-754 binary64 round-to-nearest
   for `+` and `×` in numpy/BLAS.
4. The code `cert.py` (about 245 lines) does what sections 1 to 5 say. It has been tested (section 6) but
   not read by anyone else.

## 1. The rounded scheme

Notation as in ROUTES.md: E ~ Exp(1); a law in 𝒬 is the law of Q(E), Q concave nondecreasing;
g_Q(t) = E(t + X' + X'')_+³/6 with X', X'' iid ~ Q(E); T(Q) = g_Q^(−1) on (0,∞);
Γ(Q) = E T(Q)(E) + 2·E Q(E); s_Q(θ) = g_Q'(g_Q^(−1)(θ)), so T(Q)'(θ) = 1/s_Q(θ).

**Family.** A *rounded law* is given by finitely many double-precision numbers
0 < θ_0 < … < θ_N and 0 < ŝ_0 ≤ … ≤ ŝ_N (read as exact rationals). Its quantile function Q̃ has
slope 1/ŝ_0 on (0, θ_0], slope 1/ŝ_i on (θ_{i−1}, θ_i], slope 0 on (θ_N, ∞), and is normalised by
Q̃(θ_N) = 0. Equivalently X = Q̃(E) has P(X > x) = exp(−G̃(x)) for x < 0, where G̃ is piecewise
linear and convex: G̃(x_i) = θ_i, slope ŝ_i on (x_{i−1}, x_i], x_N = 0,
x_{i−1} = x_i − (θ_i − θ_{i−1})/ŝ_i, G̃ = 0 at and below x_min = x_0 − θ_0/ŝ_0; and an atom of mass
exp(−θ_N) at 0. The support is [x_min, 0]. Q̃ is concave, nondecreasing and bounded, so Q̃ ∈ 𝒬.

**Stage map.** Q̃_0 = 0 (the law δ_0). Given Q̃_n with g = g_{Q̃_n}:
choose points 0 < t_0 < … < t_N (any choice; plain floating point), compute rigorous bounds
 θ_i ≤ g(t_i) and ŝ_i ≥ g'(t_i),
replace ŝ by its running maximum, check θ_i strictly increasing, and let Q̃_{n+1} be the rounded law
with parameters (θ_i, ŝ_i).

The slopes must be nondecreasing in i, otherwise Q̃_{n+1} is not concave and Theorem M does not
apply (the referee's required fix). The program enforces this with the running maximum
max_{k≤i} ŝ_k, which is still ≥ g'(t_i), and asserts the result. The referee's variant, the suffix
minimum min_{k≥i} ŝ_k, is also valid (g'(t_k) ≥ g'(t_i) for k ≥ i) and never larger. In practice
neither does anything: in run A the raw bounds were already strictly increasing at every stage, and
both operations changed 0 entries (`check_monotone.log`).

**Lemma R1 (each stage is safe; complete).** Q̃_{n+1} ⪯ T(Q̃_n), and Q̃_{n+1} ∈ 𝒬.
Proof. Y = X' + X'' ≤ 0 with P(Y > −t) > 0 for every t > 0, so g = 0 on (−∞, 0], g is C², convex and
strictly increasing on (0,∞), and g' is nondecreasing; hence s(θ) = g'(g^(−1)(θ)) is nondecreasing
in θ. Since g(t_i) ≥ θ_i, g^(−1)(θ_i) ≤ t_i and s(θ_i) ≤ g'(t_i) ≤ ŝ_i. For θ in (θ_{i−1}, θ_i]
(or (0, θ_0] for i = 0): s(θ) ≤ s(θ_i) ≤ ŝ_i, so T(Q̃_n)'(θ) = 1/s(θ) ≥ 1/ŝ_i = Q̃_{n+1}'(θ). For
θ > θ_N: Q̃_{n+1}' = 0 ≤ T(Q̃_n)'. Both functions are locally absolutely continuous, so
T(Q̃_n) − Q̃_{n+1} is nondecreasing on (0,∞). ∎
No tail estimate is used anywhere: below θ_0 the slope bound at θ_0 is valid for all smaller levels,
and above θ_N the rounded slope is 0.

**Lemma R2 (the chain).** Q̃_n ⪯ Q* for every n and every fixed point Q* = T(Q*) ∈ 𝒬, hence
c* = 3·Γ(Q*) ≤ 3·Γ(Q̃_n) = 3·(E T(Q̃_n)(E) + 2·E Q̃_n(E)).
Proof. 0 ⪯ Q*. If Q̃_n ⪯ Q* then Q̃_{n+1} ⪯ T(Q̃_n) ⪯ T(Q*) = Q* by Lemma R1 and Theorem M
(**assumed**). Then Lemma 3 of ROUTES.md (**assumed**). ∎

## 2. Enclosures for g and g' of a rounded law (complete)

Fix h = 2^(−HEXP) and lattice cells C_k = [kh, (k+1)h], k = K_min, …, −1, K_min = ⌊x_min/h⌋. Let μ be
the law of X = Q̃(E), S(x) = P(X > x) = exp(−G̃(x)) for x < 0, S_k = S(kh) for k < 0, S_0 = S(0−) =
exp(−θ_N), J_k = (1/h)∫_{C_k} S(x) dx. For the continuous part of μ on C_k:

 p_k := μ(C_k) = S_k − S_{k+1},  D_k := ∫_{C_k} (x − kh) dμ = h·(J_k − S_{k+1}) ≥ 0

(integration by parts), and the atom has p_0 = exp(−θ_N), D_0 = 0. Write φ_j(u) = u_+^j/j!.

**Lemma R3 (lower bounds, tangent planes).** For every real t,
 g(t) ≥ Σ_m (p∗p)_m·φ_3(t + mh) + 2·Σ_m (D∗p)_m·φ_2(t + mh),
 g'(t) ≥ Σ_m (p∗p)_m·φ_2(t + mh) + 2·Σ_m (D∗p)_m·φ_1(t + mh).
Proof. φ_3 is convex and C¹ with φ_3' = φ_2 ≥ 0. For x ∈ C_k, x' ∈ C_k':
φ_3(t + x + x') ≥ φ_3(t + (k+k')h) + φ_2(t + (k+k')h)·((x − kh) + (x' − k'h)). Integrate over
C_k × C_k' against μ⊗μ and sum over k, k' (atom included as the point cell k = 0). Same for φ_2 with
φ_2' = φ_1. ∎ All terms are nonnegative, so the inequalities survive replacing p, D by termwise
lower bounds p^lo, D^lo ≥ 0.

**Lemma R4 (upper bound, chords).** With ν_j := J_{j−1} − J_j for K_min ≤ j ≤ 0 (J_{K_min −1} := 1,
J_0 := 0), g'(t) ≤ Σ_m (ν∗ν)_m·φ_2(t + mh) for every t.
Proof. φ_2 is convex, so for x ∈ C_k and w = (x − kh)/h,
φ_2(c + x) ≤ (1 − w)·φ_2(c + kh) + w·φ_2(c + (k+1)h). Integrating over C_k gives the weights
p_k − D_k/h = S_k − J_k at kh and D_k/h = J_k − S_{k+1} at (k+1)h; adding the two contributions at a
lattice point j gives J_{j−1} − J_j (at j = K_min: 1 − J_{K_min}; at j = 0: J_{−1} − S_0 + atom =
J_{−1}). Apply this in x for fixed x'', then in x''. ∎ φ_2 ≥ 0, so ν may be replaced by termwise
upper bounds ν^hi.

Both bounds are second order in h. The final functional:

**Lemma R5.** E Q̃(E) = K_min·h + h·Σ_k J_k (because X ≤ 0 and S = 1 below x_min), and, since
T(Q̃)(E) has survival function exp(−g(t)) on t > 0 and is positive,
 E T(Q̃)(E) = ∫_0^∞ e^(−g) ≤ t_0 + Σ_{i<N} e^(−θ_i)·(1 − e^(−a_i·d_i))/a_i + e^(−θ_N)/a_N,
for any 0 < t_0 < … < t_N, d_i = t_{i+1} − t_i, θ_i ≤ g(t_i), 0 < a_i ≤ g'(t_i).
Proof. On [0, t_0] use e^(−g) ≤ 1. On [t_i, t_{i+1}] and on [t_N, ∞) use convexity:
g(t) ≥ g(t_i) + g'(t_i)(t − t_i) ≥ θ_i + a_i(t − t_i). ∎

## 3. What the program computes at stage n (cert.py)

1. `law_to_lattice` (interval arithmetic, 120 bits): encloses x_i, x_min; walks the merged sequence
   of breakpoints and lattice points (every order comparison must be decided by the enclosures,
   otherwise the program stops); encloses S at each lattice point (one `rexp`), ∫_{C_k} S as a
   sum of (S(a) − S(b))/ŝ over the linear pieces of G̃ in the cell; outputs doubles
   p^lo_k ≤ p_k, D^lo_k ≤ D_k, ν^hi_j ≥ ν_j and an upper bound on E Q̃(E).
2. `conv_nonneg` (float64): p^lo∗p^lo, D^lo∗p^lo, ν^hi∗ν^hi; deflated/inflated (section 5).
3. `choose_mesh` (plain floating point, no effect on validity): points t_i, multiples of 2^−20,
   placed so that g(t_i) is about log-spaced on [1e-6, 32] with extra points for levels in [0.01, 8].
4. `evalsums` (float64): the three sums of Lemmas R3, R4 at every t_i, giving θ_i = g^lo(t_i),
   a_i = g'^lo(t_i), ŝ_i = g'^hi(t_i).
5. `mean_T_hi` (interval arithmetic): the bound of Lemma R5; then
   3·(E T upper + 2·E Q̃ upper) in interval arithmetic, rounded up. This is the printed number.
6. Next law: parameters θ_i and running maximum of ŝ_i.

Stage 0 is the same code with the one-point law δ_0.

## 4. Size of the computation

Run A: lattice ≤ 9862 cells, 4000 levels: under 1 min for eight stages. Run B: ≤ 39539 cells,
≈ 12000 levels: about 4 min for ten stages. One process each.

## 5. Arithmetic: what is rigorous, what is floating point

- **Interval (mpmath.iv, prec = 120)**: all of step 1 and step 5. Doubles enter as exact rationals.
  Interval → double conversion uses `nextafter` in the safe direction and is asserted exactly
  (`f_hi`, `f_lo`). Only `+ − × ÷`, `abs` and comparisons of `mpmath.iv` are used (in mpmath these
  are computed from exact integer arithmetic with directed rounding; sampled test below).
- **The exponential is `rexp` in cert.py, not `iv.exp`.** exp(x) = exp(r)^(2^k) with r = x/2^k,
  |r| ≤ 2^−10 (asserted in interval arithmetic), exp(r) ∈ Σ_{j<12} r^j/j! + [−ρ, ρ],
  ρ = 2·2^−120/12! ≥ |r|^12·e^|r|/12! (Lagrange remainder), Horner and the k squarings in interval
  arithmetic. Reason: **mpmath 1.3.0's `iv.exp` is not a rigorous enclosure.** At prec 120 its upper
  endpoint was below the true value (400-bit reference) for 5 of 20000 random arguments, by a
  relative 1e−40 (about 2e−4 ulp). The defect is far too small to change any digit here (the first
  version of this certificate used `iv.exp`; its logs are in `old_ivexp/` and are identical), but a
  certificate should not rest on it. `test_rexp.py`: `rexp` encloses the reference for 40000 of
  40000 arguments (relative width ≤ 5e−31), and interval `+ − × ÷` enclose the exact rational
  result in 160000 of 160000 random cases.
- **float64 with a proved error bound**: steps 2 and 4. Every quantity is a sum of at most 2^17
  products of nonnegative numbers; the arguments t_i + m·h are exact (both are multiples of 2^−20 of
  size < 32; asserted); each term involves at most 3 rounded multiplications. For such a sum,
  evaluated in any order with or without FMA, |computed − exact| ≤ ((1 + 2^−53)^(n+3) − 1)·exact
  + n·η, where n ≤ 2^17 is the number of terms and η ≤ 2^−1022 covers underflow (even under
  flush-to-zero). So the relative error is below 3e−11 and the absolute part below 3e−303. Every
  lower bound is multiplied by (1 − 1e−9) and reduced by 1e−290 (clamped at 0); every upper bound is
  multiplied by (1 + 1e−9) and increased by 1e−290. The margin factor 30 also covers the rounding
  of these adjustments and of the divisions by 2 and 6.
- **Plain floating point**: only the choice of the t_i (and the spot checks' random indices).
- Extra evidence, not part of the proof: at every stage `spot_check` recomputes two convolution
  entries and the three sums at two mesh points in interval arithmetic and asserts that the float64
  bounds are on the correct side.

## 6. Checks (floating point; evidence, not proof)

- `check_float.py` (independent code path: non-dyadic cells of width 0.0007, mass at cell midpoints,
  FFT convolution, trapezoid rule, no directed rounding), on the saved parameters of run A
  (`check_float_A.log`):
  (B) certified laws re-evaluated: float c(Q̃_n) = 4.867955, 3.417527, 2.665064, 2.404162, 2.308348,
  2.274751, 2.262951, 2.258822, 2.257377 for n = 0..8; the certified bounds exceed these by at most
  4e−6. On every mesh, θ_cert/g_float ∈ [0.999982, 1.000000] and ŝ_cert/g'_float ∈
  [1.000000, 1.000003]: the certified parameters are on the safe side of the float values.
  (A) the scheme iterated from scratch in floats on the same meshes: 2.274746 at stage 5,
  2.257372 at stage 8.
- `test_small.py` (`test_small.log`): on a 20-piece law and lattices h = 1/8, 1/16, 1/32 the bounds
  are compared with tensor Gauss–Legendre quadrature of the exact law. g^lo/g ∈ [0.47, 0.9984],
  [0.77, 0.9996], [0.92, 0.9999]; g'^hi/g' ∈ [1.00016, 1.54], [1.00004, 1.13], [1.00001, 1.03]:
  one-sided as required and second order in h. Cell masses, cell moments and ν agree with a
  trapezoid reference to 2e−10.
- For n = 1..6 the certified values exceed the referee's unrounded reference c_n by 0.0007 to
  0.0022 (run B) and by 0.0023 to 0.0070 (run A): always above, and the loss shrinks with the mesh,
  as expected from rounding that only loses dispersion. The stage-10 bound of run B is 0.0023
  above the referee's limit 2.24954561.
- `check_monotone.py` (`check_monotone.log`): the saved parameters of run A are reproduced exactly
  from the saved laws, and the monotonisation of the slopes was inactive at every stage.

## 7. Gaps, risks, remarks for the referee of ROUTES.md

1. Theorem M is used for laws with an atom at the top (Q̃ constant above θ_N) and for δ_0. In
   Lemma 1 of ROUTES.md the survival function is then log-concave with a jump to 0; I see no
   problem, but the referee should confirm the proof covers Q that is not strictly increasing.
2. I found no error in Lemma 3 or Corollary A while using them.
3. **For other lines of this round**: `firstmoment/certify.py` uses `iv.exp` and `iv.log` (prec 80).
   By the defect above its enclosures may be too narrow by a relative amount of order 1e−28. This
   cannot matter numerically unless a margin there is of that size, but for a strict certificate
   each `iv.exp`/`iv.log` result should be widened by a few ulps (or replaced as here).
4. `cert.py` has not been read by a second person or model. The mesh choice is deterministic but
   may differ across BLAS builds; any mesh is valid, and run A's parameters are saved in `stages/`.
5. No LEDGER, CLAIMS or research-handoff entry was made from this line (writes confined to `cert/`).

## 8. Files and reproduction

`cert.py` (certificate), `run_A_h10.log`, `run_B_h12.log` (logs), `stages/A_stage{1..8}.npz`
(parameters θ, ŝ, t of run A), `check_float.py`, `check_float_A.log`, `test_small.py`,
`test_small.log`, `check_monotone.py`, `check_monotone.log`, `test_rexp.py` (checks of `rexp` and of interval `+ − × ÷`), `old_ivexp/` (logs
of the first version, which used `iv.exp`; superseded, identical numbers).

```
P=/Users/maxf/miniforge3/bin/python3
$P cert.py 10 1500 2500 8 stages/A   # run A, under 1 min
$P cert.py 12 4000 8000 10           # run B, about 4 min
$P check_float.py stages/A 8         # 7 s
$P test_small.py 4                   # 5 s
$P test_rexp.py                      # 15 s
```

## 9. Confidence

Given Theorem M and Lemma 3: 93% that the statement "c* ≤ 2.2575" (run A) is correctly certified,
and the same for 2.2519 (run B); the residual is an unnoticed slip in the code or in
Lemmas R3 to R5, which the two floating checks make unlikely (they would have to miss an error at
least ten thousand times larger than the 4e−6 agreement to matter: the margin to 2.30 is 0.043).
