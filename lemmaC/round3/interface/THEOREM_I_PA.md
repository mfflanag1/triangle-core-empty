# Theorem I-pa: the integer lower bound in the pair-additive market

Round 3, line `interface`. Written 2026-11-11. This file joins `../firstmoment/THEOREM_I.md`
(Theorem I) and `../sc/TRANSFER.md` (Theorem T, Proposition 5). Nothing here has been read by
a second agent or a different model family. Every part is labelled **complete** or **gap**.
Floating point (`toy_check.py`) is evidence only.

## 0. Result and status

**Theorem I-pa.** In the pair-additive market (3q² independent U[0,1] pair values), for
every β ≤ 2.30224, and in particular for β = 2.30 and β = 2.302,

  P( some perfect assignment has total deficit ≤ β·q^(1/3) ) → 0,

equivalently 3q − T > β·q^(1/3) with probability tending to 1. No rate is obtained.

| Step | Content | Label |
|---|---|---|
| D | Dictionary between the two files | complete |
| F | The first moment of Theorem I depends on the stars only through ν_b and Φ̂_s(u) | complete |
| S | Stability of both functionals under changing k teams | complete |
| L | Box Legendre identity H_A(μ) = H(μ') − (A/2)·‖μ − μ'‖₁ | complete |
| E | Entropy modulus | complete |
| Q | Quantitative replacement for Lemma 3.2 of Theorem I (explicit A, explicit error) | complete |
| I-a* | The interface, with explicit constants | complete, given the imports of section 7 |
| A | Assembly of Theorem I-pa | complete, given the imports of section 7 |

The certificate did not have to be changed or re-run with new parameters: the perturbation
is absorbed analytically by Lemma Q, using only the certified number Ψ_0 ≤ −3.84751328.
(A reproduction run of the unchanged script is in `out_certify_rerun.txt`.)

**Errors found in the two input files:** none in the parts used. Three things had to be
changed at the join, none of them an error:
1. Theorem I's typical set 𝒯_δ (distance with weight 1 + |d|) is not of the
   bounded-increment form that TRANSFER.md asks for: one team at an agent of degree |d| can
   move that distance by up to (3 + 2|d|)/q. It is replaced by the two functionals the
   proof actually uses (section 3), which have bounded one-team increments.
2. Lemma 3.2 of Theorem I is a compactness argument with no explicit A or δ. Lemma Q
   replaces it by an explicit bound, which is what makes the stability parameter explicit.
3. Theorem T is stated for the threshold (e_1* − ε)·q. Its proof uses the threshold only
   through the event A_N and through E'[Z] = exp(O(q)), so it holds verbatim for any fixed
   threshold β·q. The constant e_1* of TRANSFER.md is not used here.

## 1. Dictionary (complete)

| Object | THEOREM_I.md | TRANSFER.md | Here |
|---|---|---|---|
| independent-teams model | "the model", P | P', E' | P', E' |
| pair-additive market | not present | P, E | P, E |
| truncation, bins | K, B equal bins I_b of width K/B | K, bins κ_0 < … < κ_L, minimal width δ_b | K = 5, B = L = 50, δ_b = 1/10 |
| bin masses | λ_b | Λ_β = q²·p_β | λ_b |
| reference star law | ν = product of Poisson(λ_b) | π_K | π_K |
| binned near-tight hypergraph | G | H (under P), H' (under P') | H, H' |
| exponential tilt | theta | not present | ϑ = 1.6712 |
| stability fraction | not present | θ | τ |
| pair-sharing exponent | not present | η | η_T |
| count of cheap assignments | Z (X(N) < βq) | Z_ε (Σξ ≤ (e_1* − ε)q) | Z_β (X(N) ≤ βq) |

The one-team law is the same in both files: the deficit of a team is a sum of three
independent U[0,1] variables (law F), and teams are independent under P'. The labelled low
hypergraph G of Theorem I and the binned hypergraph H' of Theorem T are the same random
object. Z_β counts all perfect assignments with X(N) ≤ βq, whatever mix of teams below and
above K they use; since 1{X ≤ βq} ≤ exp(ϑ(βq − X)), section 2 of THEOREM_I.md applies to
it unchanged. A *configuration* h is a set of teams with bin labels; D(h) is its degree
array, D^s_{v,b}(h) = number of teams of h of bin b containing agent v of sector s.

## 2. Statements

Notation from THEOREM_I.md: Φ(u) = E ln(e^{u_inf} + Σ_b D_b·e^{u_b}) for the reference star,
Φ̂_s(u; D) = (1/q)·Σ_{v∈V_s} ln(e^{u_inf} + Σ_b D^s_{v,b}·e^{u_b}), ν_b(D) = n_b/q,
ℓ_b = ln(κ_b/λ_b), ℓ_inf = ln J, Ψ_0 = Ψ_0(K, B, ϑ, u_0).

New notation. U_A = {u ∈ R^(B+1): max u − min u ≤ A}. ℓ_max = max(max_b |ℓ_b|, |ℓ_inf|).
λ_max = max_b λ_b. For δ_1, δ_2, A > 0:

  𝒯(δ_1, δ_2; A) = { degree arrays D : Σ_b |ν_b(D) − λ_b| ≤ δ_1 and
            sup over u ∈ U_A of |Φ̂_s(u; D) − Φ(u)| ≤ δ_2 for s = 1, 2, 3 },
  𝒢_q(δ_1, δ_2; A) = { configurations h : D(h) ∈ 𝒯(δ_1, δ_2; A) }.

𝒢^(θ) is as in TRANSFER.md: configurations within θq teams (symmetric difference, labels
included) of a member of 𝒢.

**Interface (I-a*), precise form.** Let K, B, ϑ, u_0, β satisfy m := −(ϑβ + Ψ_0) > 0.
Let c' > 0, A = (8c' + 4 + 2ℓ_max)/3, and δ_1, δ_2, τ > 0 be such that

  c' ≥ 1 + ln^+(λ_max + δ_1 + 2τ),
  Err := 3·(δ_2 + 2τ·(A + ln 2)) + 4c'·(δ_1 + 2τ) + (6B + 2)·exp(−1 − c') < m.

Put 𝒢_q = 𝒢_q(δ_1, δ_2; A). Then
 (i) P(H ∈ 𝒢_q) → 1 in the pair-additive market;
 (ii) 𝒢_q^(2τ) ⊂ 𝒢_q(δ_1 + 2τ, δ_2 + 2τ(A + ln 2); A);
 (iii) for q^(2/3) ≥ K, E'[Z_β; H' ∈ 𝒢_q^(2τ)] ≤ C·q^(3B+2)·exp(−(m − Err)·q), with C
   depending only on K, B.

**Certified constants** (`constants.py`, interval arithmetic, output `out_constants.txt`;
K = 5, B = 50, ϑ = 1.6712, ℓ_max = ϑK = 8.356, λ_max = 7351/6000, Ψ_0 ≤ −3.84751328):

| β | m ≥ | c' | A | δ_1 | δ_2 | τ | Err ≤ | m − Err ≥ |
|---|---|---|---|---|---|---|---|---|
| 2.30 | 0.0037532 | 14 | 44.2374 | 10⁻⁵ | 2·10⁻⁴ | 5·10⁻⁶ | 0.0031604 | 0.00059 |
| 2.302 | 0.0004108 | 17 | 52.2374 | 10⁻⁶ | 2·10⁻⁵ | 2.5·10⁻⁷ | 0.0002460 | 0.00016 |

For any β < −Ψ_0/ϑ (certified ≥ 2.30224586) such constants exist: take c' large, then
δ_1, δ_2, τ small.

## 3. What the first moment uses (complete)

THEOREM_I.md (2.2) holds for an arbitrary set 𝒯 of degree arrays. With the bounds of its
sections 3.1 and 3.2, for every 𝒯 and q^(2/3) ≥ K,

  **(3.1)** E'[Z_β; D(H') ∈ 𝒯] ≤ P_q · exp( q · sup over D ∈ 𝒯 and profiles μ of Ξ(μ; D) ),

  Ξ(μ; D) = ϑβ + Σ_s inf_u [Φ̂_s(u; D) − ⟨μ, u⟩] − 2·Σ_b g(μ_b, ν_b)
       + Σ_b μ_b·ℓ_b + ε·ℓ_inf + 2ε·ln ε − 2ε,

where μ_b = m_b/q ≤ ν_b, ε = 1 − Σ_b μ_b ≥ 0, g(x, y) = y·h(x/y), and
P_q = exp(K⁶/(36q(1−p)))·(1−p)^(−q)·e²·(q+1)^(B+2)·Π_b (n_b + 1)², which is at most
C·q^(3B+2)·(1 + max_b ν_b)^(2B).

So the star profile enters only through
 (F1) the bin densities ν_b = n_b/q (in g and in the constraint μ_b ≤ ν_b), and
 (F2) the three functions u ↦ Φ̂_s(u; D), through inf_u.
Both are averages over agents of a function of the binned star: d ↦ d_b and
d ↦ ln(e^{u_inf} + Σ_b d_b·e^{u_b}). Nothing else about the stars is used. (Theorem I's set
𝒯_δ is a convenient sufficient condition for closeness of F1 and F2; by its section 3.3,
𝒯_δ ⊂ 𝒯(Bδ, (A+1)δ; A).)

## 4. Lemmas

**Lemma S (stability; complete).** Let h, h' be configurations with |h Δ h'| ≤ k. Then
Σ_b |n_b(h) − n_b(h')| ≤ k, and for every u ∈ U_A and every s,
|Φ̂_s(u; h) − Φ̂_s(u; h')| ≤ (k/q)·ln(1 + e^A) ≤ (k/q)·(A + ln 2). Hence (I-a*)(ii).

Proof. Pass from h to h' by removing the teams of h − h' and then adding those of h' − h,
one at a time. Each step changes one n_b by 1 and, in each sector, one star in one
coordinate by 1. If W = e^{u_inf} + Σ_b d_b·e^{u_b} then adding a team of bin b changes ln W
by ln(1 + e^{u_b}/W) ≤ ln(1 + e^{u_b − u_inf}) ≤ ln(1 + e^A), because W ≥ e^{u_inf}; removing
one is the same step read backwards. QED

**Lemma L (box Legendre identity; complete).** Let μ ∈ R^(B+1) with Σμ = 1 and A > 0. The
function u ↦ Φ(u) − ⟨μ, u⟩ attains its minimum over U_A at some u*, and with
μ' = grad Φ(u*) ∈ 𝒟 (the set of Lemma 4.3 of THEOREM_I.md; 𝒟 ⊂ P),

  Φ(u*) − ⟨μ, u*⟩ = H(μ') − (A/2)·‖μ − μ'‖₁.

Proof. The function f(u) = Φ(u) − ⟨μ, u⟩ is convex, differentiable and invariant under
u → u + c·1; U_A is convex and invariant under the same shifts, and U_A ∩ {Σu = 0} is
compact. So a minimiser u* exists, and ⟨grad f(u*), v − u*⟩ ≥ 0 for all v ∈ U_A. Here
grad f(u*) = μ' − μ =: z with Σz = 0. Thus ⟨z, u*⟩ = min over v ∈ U_A of ⟨z, v⟩. For
v ∈ U_A and c = (max v + min v)/2, ⟨z, v⟩ = ⟨z, v − c·1⟩ ≥ −(A/2)·‖z‖₁, with equality at
v = −(A/2)·sign(z). So ⟨μ' − μ, u*⟩ = −(A/2)·‖μ − μ'‖₁. Finally u* is a critical point of
the convex function u ↦ Φ(u) − ⟨μ', u⟩, so H(μ') = Φ(u*) − ⟨μ', u*⟩, and
f(u*) = H(μ') + ⟨μ' − μ, u*⟩. QED

**Lemma E (entropy modulus; complete).** Let φ(x) = −x·ln x on [0, ∞), φ(0) = 0. For
a, b ∈ [0, M] and c' ≥ 1 + ln^+ M: φ(a) − φ(b) ≤ c'·|a − b| + exp(−1 − c').

Proof. If a ≥ b: φ is concave with φ(0) = 0, hence subadditive, so
φ(a) − φ(b) ≤ φ(a − b), and sup over z ≥ 0 of [φ(z) − c'z] = exp(−1 − c'). If a < b:
φ(a) − φ(b) = ∫_a^b (1 + ln x) dx ≤ (b − a)·(1 + ln^+ b) ≤ c'·(b − a). QED

**Lemma Q (quantitative form of Lemma 3.2 of THEOREM_I.md; complete).** Let c' > 0,
A = (8c' + 4 + 2ℓ_max)/3, and D ∈ 𝒯(δ_1', δ_2'; A) with c' ≥ 1 + ln^+(λ_max + δ_1'). Then for
every profile μ (μ_b ≤ ν_b, ε = 1 − Σ_b μ_b ≥ 0),

  Ξ(μ; D) ≤ ϑβ + Ψ_0 + 3δ_2' + 4c'·δ_1' + (6B + 2)·exp(−1 − c').

Proof. Take u*, μ' from Lemma L for this μ. The same u* is admissible in all three infima,
and u* ∈ U_A, so by the definition of 𝒯 and Lemma L

  Σ_s inf_u [Φ̂_s(u) − ⟨μ, u⟩] ≤ 3·[Φ(u*) − ⟨μ, u*⟩] + 3δ_2' = 3H(μ') − (3A/2)·‖μ − μ'‖₁ + 3δ_2'.

μ' ∈ P, so Proposition 4.5 of THEOREM_I.md gives 3H(μ') ≤ Ψ_0 + 2S_0(μ') − ⟨ℓ, μ'⟩, where
S_0(μ') = Σ_b λ_b·h(μ'_b/λ_b) + φ(ε') + ε'. Write S̃(μ, ν) = Σ_b g(μ_b, ν_b) + φ(ε) + ε, so
that the entropy part of Ξ is −2S̃. Then

  Ξ(μ; D) − ϑβ − Ψ_0 − 3δ_2' ≤ 2·[S_0(μ') − S̃(μ, ν)] + ⟨ℓ, μ − μ'⟩ − (3A/2)·‖μ − μ'‖₁.   (4.1)

Since g(x, y) = φ(x) + φ(y − x) − φ(y),

  g(μ'_b, λ_b) − g(μ_b, ν_b) = [φ(μ'_b) − φ(μ_b)] + [φ(λ_b − μ'_b) − φ(ν_b − μ_b)] + [φ(ν_b) − φ(λ_b)].

All arguments lie in [0, λ_max + δ_1']. With z_b = |μ_b − μ'_b|, d_b = |ν_b − λ_b| and Lemma E,
this is at most c'·(2z_b + 2d_b) + 3·exp(−1 − c'). For the high coordinate, with
z_inf = |ε − ε'| and ε, ε' ∈ [0, 1], φ(ε') + ε' − φ(ε) − ε ≤ (c' + 1)·z_inf + exp(−1 − c').
Also ⟨ℓ, μ − μ'⟩ ≤ ℓ_max·‖μ − μ'‖₁. So the right side of (4.1) is at most

  (4c' + 2 + ℓ_max − 3A/2)·‖μ − μ'‖₁ + 4c'·Σ_b d_b + (6B + 2)·exp(−1 − c'),

and the first coefficient is 0 by the choice of A. QED

Remarks. (1) Lemma Q needs no compactness and no continuity of H at the boundary of its
domain. (2) The error is continuous in (δ_1', δ_2') and tends to 0 as c' → ∞, δ' → 0: this
is the statement "the exponent moves continuously" with explicit modulus. (3) Floating-point
check in toy models with 1 to 3 bins (`toy_check.py`, `out_toy_check.txt`): the identity of
Lemma L holds to 3·10⁻⁶ (optimiser tolerance) in 200 random cases, the inequality of Lemma Q
holds in all of them (largest left minus right: −0.12), Lemma E in 200,000 samples.

## 5. Proof of (I-a*) (complete, given section 7)

(ii) is Lemma S with k = 2τq.

(iii) By (ii), {H' ∈ 𝒢_q^(2τ)} ⊂ {D(H') ∈ 𝒯(δ_1', δ_2'; A)} with δ_1' = δ_1 + 2τ,
δ_2' = δ_2 + 2τ(A + ln 2). Apply (3.1) to this set and Lemma Q:
E'[Z_β; H' ∈ 𝒢_q^(2τ)] ≤ P_q·exp(q·(ϑβ + Ψ_0 + Err)) = P_q·exp(−(m − Err)·q), and
ν_b ≤ λ_max + δ_1' on this set, so P_q ≤ C·q^(3B+2).

(i) Normalising max u = 0, the supremum over U_A in the definition of 𝒯 is a supremum over
the cube [−A, 0]^(B+1) (both Φ̂_s and Φ shift by c under u → u + c·1). For u in the cube
put g_u(d) = ln(e^{u_inf} + Σ_b d_b·e^{u_b}). Then −A ≤ g_u(d) ≤ ln(1 + |d|), so
|g_u| ≤ A + degree; E_{π_K} g_u = Φ(u); and |g_u(d) − g_{u'}(d)| ≤ ‖u − u'‖_∞ for every d,
because the gradient of g_u in u is a probability vector. Hence Φ̂_s and Φ are 1-Lipschitz
in the sup norm. Take a finite net 𝒩 of the cube of mesh δ_2/4 (it does not depend on q).
By Proposition 5(d) of TRANSFER.md, applied to the finitely many functions g_u, u ∈ 𝒩, and
d ↦ d_b, b ≤ B, in each of the three sectors: with probability tending to 1 under P,
|Φ̂_s(u; D(H)) − Φ(u)| ≤ δ_2/2 for all u ∈ 𝒩 and s, and |ν_b − λ_b| ≤ δ_1/B for all b. On
this event sup over the cube of |Φ̂_s − Φ| ≤ δ_2/2 + 2·δ_2/4 = δ_2 and Σ_b |ν_b − λ_b| ≤ δ_1,
that is H ∈ 𝒢_q. QED

## 6. Assembly of Theorem I-pa (complete, given section 7)

**Order of choices.**
1. K = 5, B = 50 equal bins, ϑ = 1.6712, u_0 = the weights of `../firstmoment/w_K5_e10.json`.
   This fixes Ψ_0 ≤ −3.84751328 (certificate) and the transfer constants C1 = 36K⁵/δ_b⁶ =
   1.125·10¹¹, W = C1·q^(2/3).
2. β < −Ψ_0/ϑ, giving the margin m = −(ϑβ + Ψ_0) > 0.
3. c', then A = (8c' + 4 + 2ℓ_max)/3.
4. δ_1, δ_2, τ with Err < m (table of section 2). This fixes 𝒢_q = 𝒢_q(δ_1, δ_2; A).
5. η_T = 1/6, s_0 = q^(5/6).
6. q → ∞. (TRANSFER's ε is not chosen: the threshold is βq.)

**Chain.** 3q − T = q^(−2/3)·min_N X(N), so {3q − T ≤ β·q^(1/3)} = {Z_β ≥ 1}. The event
{H ∈ 𝒢_q, s(H) ≤ s_0, |H| ≤ K³q} is a function of H, so

  P(Z_β ≥ 1) ≤ P(H ∉ 𝒢_q) + P(s(H) > s_0) + P(|H| > K³q)
         + E[Z_β; H ∈ 𝒢_q, s(H) ≤ s_0, |H| ≤ K³q].

| Term | Bound | Source |
|---|---|---|
| P(H ∉ 𝒢_q) | o(1), no rate | (I-a*)(i), which uses Proposition 5(d) |
| P(s(H) > q^(5/6)) | (3/40)·K⁵·q^(−1/6) | Proposition 5(b) |
| P(\|H\| > K³q) | O(1/q) | Proposition 5(a) |
| E[Z_β; …] | exp(R_q)·( E'[Z_β; H' ∈ 𝒢_q^(2τ)] + 3·E'[Z_β]·exp(−(τ/4)·q·ln q) ) | Theorem T with θ = τ, η = 1/6, threshold βq |

with R_q = Δ + 2K⁶ + s_0·ln W + K³W/2 ≤ (3/40)K⁵q^(2/3) + 2K⁶ + q^(5/6)·ln(C1·q^(2/3))
+ (K³C1/2)·q^(2/3) = O(q^(5/6)·ln q). The two terms in the bracket:

- E'[Z_β; H' ∈ 𝒢_q^(2τ)] ≤ C·q^(3B+2)·exp(−(m − Err)·q) by (I-a*)(iii). The factors
  exp(o(q)) inside this bound are: the number of profiles (q+1)^B; the binomial and
  factorial corrections e²(q+1)²·Π(n_b+1)²; the simplicity factor exp(K⁶/(36q(1−p))); and
  (1−p)^(−q). All are in P_q.
- E'[Z_β] ≤ (q!)²·(β·q^(1/3))^(3q)/(3q)! ≤ e²(q+1)²·(β³e/27)^q ≤ e²(q+1)²·exp(0.22q) for
  β ≤ 2.31 (fact 4 of the brief; the 3q pair deficits of one assignment are independent
  under P, and the team deficits have the same joint law under P').

Hence
  E[Z_β; …] ≤ exp(R_q + O(ln q) − (m − Err)·q) + exp(R_q + 0.22q + O(ln q) − (τ/4)·q·ln q) → 0,
because R_q = o(q), m − Err > 0 and τ > 0. All four terms tend to 0. QED

Every excluded event: {H ∉ 𝒢_q}, {s(H) > s_0}, {|H| > K³q} under P (table above), and
{X_N > τq} under P' inside Lemma 4 of TRANSFER.md (the term with exp(−(τ/4) q ln q)).
Theorem I's own exceptional event {D(G) ∉ 𝒯_δ} under P' does not appear: typicality is
needed only under P.

Size of q. The statement is asymptotic only. The term (K³C1/2)·q^(2/3) = 7.0·10¹²·q^(2/3) of
R_q is below (m − Err)·q = 5.9·10⁻⁴·q (β = 2.30) only for q > 1.6·10⁴⁸. The tail term of
Lemma 4 of TRANSFER.md needs more: its base e·K³·C1/(2τ·q^(1/3)) is below 1 only for
q^(1/3) > 3.8·10¹⁸, that is q > 5.5·10⁵⁵ (τ = 5·10⁻⁶).

**The checker's three constraints on (I-a*)** (received from the lead after the check of
TRANSFER.md), and where each is met:
- (a) *The saving must be proved on the 2τ-enlargement with τ small.* It is: τ = 5·10⁻⁶
  (β = 2.30) or 2.5·10⁻⁷ (β = 2.302), and (I-a*)(iii) is a bound on 𝒢_q^(2τ) itself. Lemma Q
  makes the limit visible: the cost of the enlargement in the exponent is
  2τ·(3(A + ln 2) + 4c') ≈ 382·τ at β = 2.30, so the method tolerates τ up to about 10⁻⁵
  and nothing near 2τ = 1.
- (b) *𝒢 may use only statistics with bounded one-team increments that are typical by
  Proposition 5.* It uses two families of per-sector averages: d ↦ d_b (increment 1) and
  d ↦ g_u(d), u ∈ U_A (increment at most A + ln 2), both bounded by A + degree, so
  Proposition 5(d) applies. There is no maximum-degree condition and no truncation of
  stars anywhere in the argument.
- (c) *All of Z; exact one-team law; exactly q³ Bernoulli teams.* Z_β counts every
  assignment with X(N) ≤ βq through the high option of THEOREM_I.md (2.1). The model of
  THEOREM_I.md section 1.1 is q³ independent teams with deficit law U + U' + U''; its
  Lemma 2.2 uses the exact Bernoulli probability Π p_b^{n_b}·(1 − p)^{q³ − n}. No Poisson
  number of teams is introduced; "Poisson" refers only to the limit law of a star.

## 7. Hypotheses imported (for a referee)

From `../firstmoment/THEOREM_I.md` (independent-teams model, law F):
- T1. (F1), (F2) and inequality (2.1): tilt, split into low and high teams, (r!)² count for
  the teams above K.
- T2. Lemma 2.1 (exact configuration count) and Lemma 2.2 (simplicity conditioning), giving
  (2.2) for an arbitrary set of degree arrays.
- T3. Section 3.1 (choice-function bound, any weight vector) and section 3.2 (binomial and
  factorial bounds), giving the function Ξ.
- T4. Lemmas 4.1 to 4.3 and Propositions 4.4, 4.5: 3H − 2S_0 + ⟨ℓ, ·⟩ ≤ Ψ_0 on P, for every
  u_0. I re-derived the algebra of 4.1, 4.2, 4.4 and 4.5 while writing this and found no
  error; this is a reading by the same model family, not an independent check.
- T5. Section 5: Ψ_0 ≤ −3.84751328 for the stated parameters (interval arithmetic, one
  author). Reproduced here by running the unchanged script: `out_certify_rerun.txt`.
- Not used: Lemma 3.1 (typical stars under P'), section 3.3, Lemma 3.2 (compactness).

From `../sc/TRANSFER.md`:
- X1. Theorem T (Lemmas 1 to 4), read with the threshold βq in place of (e_1* − ε)q.
  I re-read the proofs of Lemmas 1, 3, 4 and section 8 and found no error.
- X2. Proposition 5(a), (b).
- X3. Proposition 5(d) for the functions d ↦ d_b and d ↦ g_u(d), u in a finite net. Both are
  bounded by A + degree. This is the only place where local convergence of stars in the
  pair-additive market is used.
- Not used: Proposition 5(c), the constant e_1*.

From the brief: fact 4 (plain first moment), re-derived in section 6.

New here: Lemmas S, L, E, Q and the constants of `constants.py`.

## 8. Gaps and confidence

- **Gap (verification, not mathematics):** Lemmas S, L, E, Q and the assembly have one
  author and no second reader. The imports T1 to T5 and X1 to X3 are being checked by
  other lines; Theorem I-pa is exactly as strong as they are.
- **No rate.** P(H ∉ 𝒢_q) → 0 comes from Proposition 5(d), which has no rate.
- **Asymptotic only**, q above about 10⁴⁸ (section 6).
- The explicit τ is small because Lemma E is used with the same constant in every bin; it
  has not been optimised and its size plays no role.

Confidence: the interface (sections 3 to 5) is correct, 92%. Theorem I-pa follows from the
imports as stated, 90%. Theorem I-pa is true with constant 2.302, 85% (the difference is the
risk in the imports).

## Files

- `THEOREM_I_PA.md`: this file.
- `constants.py`, `out_constants.txt`: interval-arithmetic error budget (section 2 table).
- `toy_check.py`, `out_toy_check.txt`: floating-point check of Lemmas L, Q, E (evidence).
- `out_certify_rerun.txt`: reproduction of `../firstmoment/certify.py` on `w_K5_e10.json`.
