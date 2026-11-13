# Existence of a fixed point of X = min_k(xi_k − X_k' − X_k''): state after rdeA

Round 3, line `rdeA`, 60-minute box. **The existence theorem is not proved.** What exists:
an exact identity that turns the shape map into an average of Markov kernels (complete), a
reduction of existence plus cost enclosure to a finite list of inequalities on an explicit
reference law (bulk part written out; tail part specified, not evaluated), and a float64
measurement of every constant in that list. No inequality has been verified in interval or
exact arithmetic. Labels: **complete** = I would defend it to a referee; **gap** = not proved.

| Item | Status |
|---|---|
| Lemma 1, exact displacement identity (§1) | **complete** |
| Lemma 2, cell-wise kernel envelope (§2) | **complete** |
| Theorem A, bulk invariance from (H1)–(H4) (§3) | proof written; **gap** until the tail term (H5) and the Schauder wrapper are written in full |
| (H1)–(H4) for the reference law of §4 | **float64 only** (measured, with margin) |
| (H5) tail envelope | **gap**: not evaluated; reference-only value measured |
| Cost enclosure | float, conditional: c ∈ [2.2465, 2.2526] |

## 0. Setting

Laws are represented by quantile functions Q on (0,1) (nondecreasing). U, U1, U2 are
independent uniforms. T(Q) is the quantile function of the law with survival function
exp(−G_Q(x)), G_Q(x) = E(x + Q(U1) + Q(U2))_+^3/6. T(Q + c) = T(Q) − 2c.

Function space: quantile functions with the weighted sup norm sup_u |φ(u)|/W(u) on a bulk
quantile range, shifts removed by an explicit linear functional. This is a weighted W_∞
metric on laws modulo translation. (L²(du), i.e. W_2, does not work: the measured norm of the
centred linearisation there is 0.87–0.94 before any envelope loss; the weighted sup norm
gives 0.351, equal to the spectral value.)

## 1. Lemma 1 (exact displacement identity) — complete

Let Q_0, Q be quantile functions with E Q(U)_+^3 < ∞, u ∈ (0,1), x_0 = T(Q_0)(u),
x = T(Q)(u), ξ = x − x_0. Put B = x_0 + Q_0(U1) + Q_0(U2), A = x + Q(U1) + Q(U2) and

  κ = ∫_0^1 3·(B + t(A − B))_+^2 dt ≥ 0,  so that A_+^3 − B_+^3 = (A − B)·κ.

Both G-values equal −ln(1−u), so E[A_+^3 − B_+^3] = 0. For any event E symmetric in (U1,U2),
with φ = Q − Q_0 and τ_E = E[A_+^3 − B_+^3; E^c]:

  ξ·E[κ; E] + 2·E[κ·φ(U1); E] + τ_E = 0.

Hence, with E = {U1 ∈ 𝔅, U2 ∈ 𝔅} for a set 𝔅 ⊂ (0,1) and
p̃(u,v) = E[κ | U1 = v; U2 ∈ 𝔅]/E[κ; E] (a probability density on 𝔅 for each u):

  **T(Q)(u) − T(Q_0)(u) = −2·∫_𝔅 p̃(u,v)·φ(v) dv − τ_E/E[κ; E].**   (1)

Proof: the display is the definition of κ plus symmetry of κ in (U1, U2). ∎
With 𝔅 = (0,1) and φ → 0 this is DT = −2P of RDE.md §2.4. The identity is nonlinear and
exact: the kernel depends on Q only through the argument of κ, and since s ↦ 3s_+^2 is
nondecreasing,

  3(B − δ)_+^2 ≤ κ ≤ 3(B + δ)_+^2 whenever |ξ + φ(U1) + φ(U2)| ≤ δ.   (2)

## 2. Reference law and Lemma 2 (kernel envelope) — complete

Reference μ_0: density constant on cells [Jh, (J+1)h), masses p_J, support [A_0, B_0];
a = its quantile function; b = T(a). Bulk: cells with [Jh,(J+1)h) ⊂ [α, β]; 𝔅 = a^{-1}[α, β].
g_B(t) = E[(t + X_0)_+^2; X_0 ∈ [α, β]]/2 is nondecreasing, and at grid points
g_B(Mh) = h²·Σ_J p_J·c3(M + J), c3(s) = (3s² + 3s + 1)/6 for s ≥ 0, else 0. Likewise
G_0(Ih) = h³·Σ_K (p∗p)_K·c5(I + K), G_0'(Ih) = h²·Σ_K (p∗p)_K·c4(I + K), with
c5(s) = (20s³ + 60s² + 70s + 30)/120, c4(s) = (12s² + 24s + 14)/24 for s ≥ 0,
c5(−1) = 1/120, c4(−1) = 1/24, zero below. If the p_J are dyadic rationals all these grid
values are exact rationals.

**Lemma 2.** Suppose a(u) ∈ [Ih, (I+1)h], |b(u) − a(u)| ≤ d, |φ| ≤ ρ on 𝔅, |ξ| ≤ ξ⁺ and
d + ξ⁺ + 2ρ ≤ h. Then for every bulk cell J and every v with a(v) in cell J (pointwise in v,
which is what the proof of Theorem A uses, since φ may vary inside a cell),
P⁻_{IJ} ≤ p_J·p̃(u,v) ≤ P⁺_{IJ}, where N⁻_{IJ} = p_J·g_B((I + J − 1)h),
N⁺_{IJ} = p_J·g_B((I + J + 3)h), P⁻_{IJ} = N⁻_{IJ}/Σ_J' N⁺_{IJ'}, P⁺_{IJ} = N⁺_{IJ}/Σ_J' N⁻_{IJ'}.
Proof: by (2) with δ = ξ⁺ + 2ρ, E[κ | U1 = v; U2 ∈ 𝔅] lies between 6·g_B(x_0 + a(v) ∓ δ);
x_0 + a(v) ranges in [(I+J)h − d, (I+J+2)h + d]; g_B is monotone. ∎

## 3. Theorem A (reduction) — bulk part written; tail part a gap

Data: μ_0 as above; weights W_J ≥ 1 and a probability vector q_J on bulk cells (W, q extended
as step functions of u through a); numbers m_0, r > 0; per-cell numbers Y_I, t_I.
K_0 = {Q nondecreasing: |Q − a| ≤ r·W on 𝔅; tail envelopes (K2), (K3) below}.
R(Q) = T(Q) − m'(Q), m'(Q) = m_0 − 2·Σ_J q_J·(mean of Q − a over cell J).

Hypotheses.
- (H1) defect: |b(u) − a(u) − m_0| ≤ Y_I for a(u) in bulk cell I.
- (H2) invariance: Y_I + r·Σ_J 2·max(P⁺_{IJ} − q_J, q_J − P⁻_{IJ})·W_J + t_I ≤ r·W_I for all I.
- (H3) envelope room: 6·r·W_max + 2·t_max + sup_𝔅|b − a| ≤ h.
- (H4) q is a probability vector, W ≥ 1 (trivial once fixed).
- (H5) tails: for every Q ∈ K_0 and u ∈ 𝔅 in cell I, |τ_E|/E[κ; E] ≤ t_I; and R maps the
  tail envelopes into themselves. Proposed envelopes: (K2) for u above 𝔅,
  Q(u) ≤ c_U + (G_0^{BB})^{-1}(−ln(1−u)) with G_0^{BB}(x) = E[(x + a(U1) + a(U2))_+^3; E]/6 and
  c_U ≥ 2rW_max + sup|m'|, which is preserved because G_Q ≥ G_0^{BB}(· − 2rW_max);
  (K3) for u below 𝔅, P(Q < −s) ≤ E(Q⁺(U1) + Q⁺(U2) − s)_+^3/6 with Q⁺ the pointwise upper
  bound of K_0. Then |τ_E| ≤ max(E[A_+^3; E^c], E[B_+^3; E^c]) is bounded by explicit one- and
  two-dimensional integrals of the envelopes. **Not evaluated.**

Conclusion (given H1–H5): R(K_0) ⊂ K_0; K_0 is convex and compact in L¹(du) and R is
continuous (Helly plus the uniform tail envelopes; **written only at this level of detail:
gap**), so by Schauder R has a fixed point Q, T(Q) = Q + c with c = m'(Q), and Q + c/3 is a
fixed point of T (RDE.md 2.2) with

  |E X* − E X_0 − m_0/3| ≤ r·(Σ_J p_J W_J + (2/3)·Σ_J q_J W_J) + (tail integrals of K2, K3).

Proof of the bulk invariance. For Q ∈ K_0 and u ∈ 𝔅, (1) with Q_0 = a gives
T(Q)(u) − a(u) − m'(Q) = (b − a − m_0)(u) − 2∫_𝔅 (p̃(u,v) − q(v))·φ(v) dv − τ_E/E[κ;E].
The middle term is at most 2r·Σ_J max(P⁺ − q, q − P⁻)·W_J by Lemma 2, provided |ξ| ≤ ξ⁺ with
d + ξ⁺ + 2rW_max ≤ h. By (1), p̃ being a probability kernel, |ξ| ≤ 2rW_max + t_max =: ξ⁺ as soon as (H5) holds for
the actual ξ, so 4rW_max + t_max + d ≤ h is enough here. Evaluating (H5) will itself need an
a priori bound on ξ; the intended argument is continuity along Q_s = (1−s)a + sQ ∈ K_0
(ξ_0 = 0, ξ_s continuous, assume |ξ_s| ≤ 2ξ⁺ and recover |ξ_s| ≤ ξ⁺), which is why (H3) is
stated with the factor 6. With q(v) = q_J/p_J on cell J, the pointwise form of Lemma 2 gives
∫_{cell J} |p̃ − q| dv ≤ max(P⁺_{IJ} − q_J, q_J − P⁻_{IJ}). Then (H1), (H2) give |R(Q) − a| ≤ rW on 𝔅. ∎

## 4. Measured constants (float64, `measure.py`, log in `measure_variants.log`)

Reference: cell masses of the saved solution `rde/fixed_point_h005.npy` on h = 0.01,
support [−6, 7], m_0 = 0. q = Perron left vector of the midpoint kernel.
W = (1 − A_env/θ)^{-1}·1 with A_env the envelope matrix of (H2), so (H2) holds with factor
≤ θ by construction once the spectral radius of A_env is below θ.

| quantity | bulk [−4.2, 4.4], θ = 0.6 | bulk [−4.5, 4.5], θ = 0.7 |
|---|---|---|
| spectral radius, midpoint kernel 2·abs(P − q) | 0.3509 | 0.3509 |
| spectral radius, envelope matrix A_env | 0.4050 | 0.4050 |
| Z = max_I Σ_J A_env·W_J / W_I | 0.5713 | 0.6531 |
| W range | [1.50, 20.9] | [1.34, 14.9] |
| defect sup of abs(b − a) (3 samples per cell) | 1.88e-4 | 1.95e-4 |
| max_I Y_I/W_I | 3.44e-5 | 4.96e-5 |
| reference tail term (G_0 − G_0^{BB})/G_0^{BB}' , max | 2.29e-6 | 1.68e-6 |
| radius r (tail term ×4 as placeholder for H5) | 6.76e-5 | 1.04e-4 |
| (H3) left side, must be ≤ 0.01 | 8.58e-3 | 9.42e-3 |
| cost half-width 9·r·(ΣpW + (2/3)ΣqW) | 3.05e-3 | 3.63e-3 |
| c enclosure (conditional, float) | [2.24649, 2.25260] | [2.24591, 2.25318] |

9·E X_0 = 2.24954561. Target 2.30 is 0.050 above; the measured half-width is 0.003.

## 5. What is and is not verified

Rigorous: Lemma 1, Lemma 2 (pencil). Nothing numerical.
Float64 only: every number in §4. Floating point enters everywhere in `measure.py`.
Not done:
1. (H1) rigorously. The sampled defect is second order in h inside a cell. A first-order
   bound (density of T(μ_0) against p_I/h on the cell) would give about 1.6e-3 at the bulk
   edges, which breaks (H3) at h = 0.01. Needs either a second-order within-cell bound or a
   reference on h ≤ 0.0025 for the defect with merged cells for the kernel sums.
2. (H2) in exact rational or interval arithmetic: one 860 × 860 sum per row, on exact grid
   values of g_B; needs dyadic p_J.
3. (H5): evaluate the envelope integrals; the measured reference value 2e-6 against
   r·W ≈ 1e-4 leaves a factor of about 40.
4. Compactness and continuity for Schauder, written in full; tail integrals in the cost bound.
5. (H3) has only 6–14% room at h = 0.01; the factor 6 is crude (v-dependent shifts would
   give room), or shrink r via item 1.

## 6. Files

`measure.py` (run: `python3 measure.py α β θ`), `measure.log`, `measure_variants.log`,
`explore1.py` (sup-norm versus Hilbert–Schmidt constants).
