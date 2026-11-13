# Step (I-b): transfer of the star-conditioned first moment to the pair-additive market

Round 2, line `sc`. Written 2026-11-11. Nothing here has been read by another agent.
Floating point (section 9) is evidence only. Sources: see section 4.

## 0. Result and status

**Proved here (complete, sections 3-8):** Theorem T, a comparison inequality between the
pair-additive market and the independent-teams ("Poissonised") model, and Proposition 5
(the good event is typical in the pair-additive market).

**Conditional conclusion:** the target statement follows from Theorem T, Proposition 5 and
the Poissonised first moment in the form (I-a*) of section 2. (I-a*) is the other agent's
lemma; it is the only input not proved here.

| Step | Content | Status |
|---|---|---|
| L1 | Cluster bound: factor C1·q^(2/3) per unit of pair-sharing excess, uniform in the values of N's pairs | complete |
| L2 | Janson-type bound for "no near-tight team in a family", proved directly | complete |
| L3 | Planted pointwise domination of the law of H − N given all pair values of N | complete |
| L4 | Removing the weight from near-tight teams that touch N (product structure of the Poissonised model) | complete |
| P5 | |H|, pair-sharing excess and empirical star law are typical in the pair-additive market | complete |
| T | Theorem T | complete |
| I | Interface: (I-a) must be delivered in the form (I-a*) (same one-team law; typical event stable under changing θq teams; exponentially small bound) | gap (owed by line `firstmoment`) |

**The orchestrator's sketch needed two repairs.** (a) A domination of the law of H alone
is not enough, because a cheap assignment uses teams above K and their costs are not a
function of H; the domination has to be for the law of H − N *given* the pair values of the
candidate assignment N. (b) In that planted form the factor is not exp(o(q)) on the good
event: it contains W^X with W = C1·q^(2/3) and X the number of near-tight teams through
pairs of N, and X is not controlled by any event of H alone. It is removed in expectation
on the Poissonised side (L4), at the price of a stability requirement on the typical
event. Everything else in the sketch survives.

## 1. Setting and notation

Pair deficits d_π, π ranging over the 3q² pairs, are independent U[0,1]. A team e = (a,b,c)
has the three pairs AB(a,b), BC(b,c), AC(a,c), deficit D_e = sum of its three pair
deficits, and rescaled cost ξ_e = q^(2/3)·D_e. F is the law of D_e; F(t) = t³/6 on [0,1].

- **P** is the pair-additive market. **P'** is the independent-teams model: the q³ team
  deficits are independent with law F. (P' has no pair values.) E, E' are the expectations.
- Fix K ≥ 1 and bins 0 = κ_0 < κ_1 < … < κ_L = K, of minimal width δ_b. Take q so large
  that K·q^(−2/3) ≤ 1. Put p_β = P(ξ_e in bin β) = (κ_β³ − κ_(β−1)³)/(6q²) ≥ δ_b³/(6q²),
  p_K = K³/(6q²), Λ_β = q²·p_β, and t = K·q^(−2/3).
- **H** is the binned near-tight hypergraph: the set of teams with ξ ≤ K, each labelled
  with its bin. H' is the same object under P'. A *configuration* h is a set of teams with
  bin labels. pairs(h) is the set of pairs lying in a team of h, and
  s(h) = 3|h| − |pairs(h)| = Σ_π (mult_h(π) − 1)^+ is its **pair-sharing excess**; s is
  monotone in h.
- N is a perfect assignment. Its q teams are disjoint, so they use 3q distinct pairs, and
  (ξ_e)_{e∈N} =: ξ_N is an i.i.d. family with law F under P and under P'. d_N is the vector
  of the 3q pair deficits of N. H⁻ = H − N is the binned set of near-tight teams outside N;
  H is a function of (ξ_N, H⁻).
- **Geometry of N.** (G1) Two distinct teams share at most one pair. (G2) A team f ∉ N
  contains at most one pair of N: if AB(a,b) and BC(b,c) are both pairs of N then N's team
  at b is (a,b,c) = f; the other two cases are the same. (G3) So the family 𝒯_N of teams
  outside N containing a pair of N has exactly 3q(q−1) members, and each pair of N lies in
  q − 1 of them.
- For h⁻ disjoint from N: s(h⁻;N) = 3|h⁻| − |pairs(h⁻) − pairs(N)|, the excess when N's
  pairs count as already used. With h_touch = h⁻ ∩ 𝒯_N, by (G2),
  **s(h⁻;N) = s(h⁻) + |pairs(h⁻) ∩ pairs(N)| ≤ s(h⁻) + |h_touch|.**
- Z_ε is the number of perfect assignments N with Σ_{e∈N} ξ_e ≤ (e_1* − ε)·q, i.e. total
  deficit at most (e_1* − ε)·q^(1/3).
- For a set 𝒢 of configurations and θ > 0, 𝒢^(θ) is the set of configurations that differ
  from some member of 𝒢 in at most θq teams (symmetric difference, labels included).
- Constants: C1 = 36·K⁵/δ_b⁶, W = C1·q^(2/3) ≥ 1, Δ = 3q³(q−1)·t⁵/40 ≤ (3/40)·K⁵·q^(2/3).

## 2. Statements

**Theorem T (complete).** Fix K, the bins, ε > 0, θ > 0 and η in (0, 1/3). Let
s_0 = q^(2/3+η) and let 𝒢 = 𝒢_q be any sets of configurations. For all q ≥ q_0(K, bins, θ),

  E[ Z_ε ; H ∈ 𝒢, s(H) ≤ s_0, |H| ≤ K³q ]
   ≤ exp(R_q) · ( E'[ Z_ε ; H' ∈ 𝒢^(2θ) ] + 3·E'[Z_ε]·exp(−(θ/4)·q·ln q) ),

  R_q = Δ + 2K⁶ + s_0·ln W + K³W/2 = O(q^(2/3+η)·ln q) = o(q).

Consequently

  P(Z_ε ≥ 1) ≤ P(H ∉ 𝒢) + P(s(H) > s_0) + P(|H| > K³q)
        + exp(R_q)·E'[Z_ε; H' ∈ 𝒢^(2θ)] + exp(−(θ/5)·q·ln q).

(E'[Z_ε] = E[Z_ε] ≤ (q!)²·(e_1*·q^(1/3))^(3q)/(3q)! = exp(O(q)) is the plain first moment,
fact 4 of the brief.)

**Proposition 5 (complete).** Under P, for fixed K, bins, δ > 0, η > 0, D:
(a) P(|H| > K³q) = O(1/q); (b) P(s(H) > q^(2/3+η)) ≤ (3/40)·K⁵·q^(−η);
(c) for each sector σ, the empirical law L_σ(H) of the binned stars of the q agents of σ
converges in total variation, in probability, to π_K, the law of independent
Poisson(Λ_β) bin counts; (d) for every function g of the binned star with
|g| ≤ c·(1 + degree)^m, (1/q)·Σ_{v∈σ} g(star_v) → E_{π_K} g in probability.

**Interface (I-a*) (gap; owed by the Poissonised first moment).** For every ε > 0 there are
K, bins, θ > 0, c > 0 and sets 𝒢_q with (i) P(H ∈ 𝒢_q) → 1 in the pair-additive market and
(ii) E'[Z_ε; H' ∈ 𝒢_q^(2θ)] ≤ exp(−c·q) for large q, where P' has one-team law F.
Remarks on what this asks:
- Z_ε counts all cheap assignments, whatever mix of teams below and above K they use. The
  sum over profiles, the binning error and the truncation at K stay inside (I-a).
- (i) holds by Proposition 5 for any 𝒢 of the form
  {d_TV(L_σ(H), π_K) ≤ δ for all σ} ∩ {|(1/q)Σ_{v∈σ} g_i(star_v) − E_π g_i| ≤ δ, i ≤ r},
  with g_i polynomially bounded. Changing one team changes one star per sector, so if
  each g_i changes by at most Lip_i < ∞ when one team is added to or removed from a star,
  𝒢^(2θ) is contained in the same event with δ replaced by δ + 2θ·max(1, max_i Lip_i).
  (For d_TV the change is at most 1/q per team. Lip ≤ D + 1 for g = degree·1{degree > D};
  Lip ≤ ln(1 + max w/w_∞) for g = ln(w_∞ + Σ_{star} w), the function in the checker's
  truncated entropy. An unbounded-increment statistic such as degree² is typical by
  Proposition 5(d) but is not stable in this sense.) So (I-a*) is the usual statement "the star-conditioned first
  moment is exponentially small on a δ-neighbourhood of the typical star law, for small
  δ", and a proof of (I-a) by conditioning on stars and summing over typical star
  sequences gives it as long as its error terms are continuous in δ.
- If (I-a) is written for a model with a Poisson number of teams, the ratio of the two
  probabilities of {H' = h} is exp(O(K³|h|/q²)) = exp(O(1)); this is not an issue.

**Corollary (target statement, conditional on (I-a*)).** For every ε > 0, whp the
pair-additive market has no perfect assignment of total deficit below
(e_1* − ε)·q^(1/3). Proof: Theorem T with η = 1/6, (I-a*) and Proposition 5(a),(b).

## 3. Lemma 1: cluster bound (complete)

**Lemma 1.** Let N be a perfect assignment or N = ∅, and let h⁻ be a configuration
disjoint from N. Then, almost surely in d_N,

  P( ξ_f in bin β_f for all f ∈ h⁻ | d_N ) ≤ W^{s(h⁻;N)} · Π_{f∈h⁻} p_{β_f}.

Proof. Call f ∈ h⁻ *isolated* if none of its pairs is a pair of N or of another team of
h⁻. Let I be the isolated teams and R the rest. The 3|I| pairs of I are distinct and are
disjoint from pairs(R) and pairs(N), so the events of I are independent of each other, of
d_N and of the events of R, and each has probability exactly p_{β_f}. Let Q_R be the pairs
of R outside pairs(N). The event for R forces d_π ≤ t for every π ∈ Q_R, which has
probability t^{|Q_R|} and does not involve d_N. Hence the left side is at most
Π_{f∈I} p_{β_f} · t^{|Q_R|}.

Write s = s(h⁻;N) = 3|R| − |Q_R| (isolated teams contribute 0). s equals the number of
slots (f,π) with π ∈ pairs(N) plus Σ_{π∈Q_R} (mult(π) − 1). Every team of R has a slot in
pairs(N) or a free pair of multiplicity at least 2; the former are at most the number of
such slots and the latter at most Σ_{mult ≥ 2} mult(π) ≤ 2·Σ (mult(π) − 1). So |R| ≤ 2s.
Using p_β ≥ δ_b³/(6q²),

  t^{|Q_R|} / Π_{f∈R} p_{β_f} ≤ (K q^(−2/3))^{3|R|−s} · (6q²/δ_b³)^{|R|}
   = (6K³/δ_b³)^{|R|} · K^{−s} · q^{2s/3} ≤ (36K⁵/δ_b⁶ · q^(2/3))^s = W^s,

since 6K³/δ_b³ ≥ 1 and |R| ≤ 2s. QED

Remarks. (1) The bound is uniform in d_N; this is what allows conditioning on the exact
values of N's pairs. (2) The lowest bin costs nothing extra: only δ_b enters. (3) The true
factor for a single sharing, without bins, is P(both ≤ K)/p_K² = (t⁵/20)/(t³/6)² =
1.8·q^(2/3)/K, so the q^(2/3) per unit of excess is the right order and cannot be improved.

## 4. Lemma 2: Janson-type bound (complete, proved directly)

**Lemma 2.** Let ℱ be any family of teams. Then

  P( ξ_e > K for all e ∈ ℱ ) ≤ exp( −|ℱ|·p_K + Δ ).

Proof. Let B_e = {D_e ≤ t}. It is a down-set in the deficits (an up-set in the pair
values). Write e ~ f if e ≠ f share a pair. B_e is independent of the σ-field generated by
{B_f : f ≁ e, f ≠ e}, because those teams use none of e's three pairs. Order ℱ as
e_1, …, e_n. For fixed i let 𝒟 = ∩_{j<i, e_j ~ e_i} B_j^c and ℐ = ∩_{j<i, e_j ≁ e_i} B_j^c.
Then

  P(B_i ∩ 𝒟 ∩ ℐ) ≥ P(B_i ∩ ℐ) − Σ_{j<i, e_j~e_i} P(B_i ∩ B_j ∩ ℐ)
          ≥ P(B_i)·P(ℐ) − Σ_{j<i, e_j~e_i} P(B_i ∩ B_j)·P(ℐ).

The first term is independence. The second is Harris's inequality in the product space of
pair deficits: B_i ∩ B_j is decreasing and ℐ is increasing. Since P(𝒟 ∩ ℐ) ≤ P(ℐ),

  P(B_i | 𝒟 ∩ ℐ) ≥ P(B_i) − Σ_{j<i, e_j~e_i} P(B_i ∩ B_j)

(if P(𝒟 ∩ ℐ) = 0 the product below is 0 and there is nothing to prove). So
P(B_i^c | ∩_{j<i} B_j^c) ≤ exp(−p_K + Σ_{j<i, e_j~e_i} P(B_i ∩ B_j)), and multiplying over i
gives exp(−|ℱ|p_K + Δ_ℱ) with Δ_ℱ the sum of P(B_e ∩ B_f) over unordered sharing pairs in
ℱ. For the full family the number of unordered sharing pairs is 3q³(q−1)/2 by (G1), and
for a sharing pair P(B_e ∩ B_f) = ∫_0^t ((t−u)²/2)² du = t⁵/20. So Δ_ℱ ≤ Δ. QED

Note that no conditioning on B_i is used, so no correlation inequality is needed for the
law on the simplex {D_e ≤ t}, which is not positively associated.

Source status. I did not open Riordan-Warnke myself. A fetch helper read the HTML of
arXiv:1203.1024v2 (Random Structures and Algorithms 46 (2015) 391-395) and reported
Theorem 1 as P(X = 0) ≤ exp(−μ + Δ/2) with Δ summed over *ordered* dependent pairs, for
events in a family closed under ∩ and ∪ with pairwise positive correlation. That is the
same inequality as Lemma 2. The proof above is self-contained and does not rely on it.

## 5. Lemma 3: planted pointwise domination (complete)

**Lemma 3.** For every perfect assignment N and configuration h⁻ disjoint from N, almost
surely in d_N,

  P( H⁻ = h⁻ | d_N ) ≤ exp( Δ + K³(|h⁻|+q)/(2q) + K⁶/(36q) ) · W^{s(h⁻;N)} · P'( H'⁻ = h⁻ ).

Proof. Let ℱ be the teams with no pair in pairs(h⁻) ∪ pairs(N). Then
{H⁻ = h⁻} ⊂ E_1 ∩ E_2 with E_1 = {ξ_f in bin β_f for all f ∈ h⁻} and
E_2 = {ξ_f > K for all f ∈ ℱ}. E_1 is determined by the deficits of pairs(h⁻) ∪ pairs(N)
and E_2 by the other pairs, so P(E_1 ∩ E_2 | d_N) = P(E_1 | d_N)·P(E_2). Lemma 1 bounds the
first factor and Lemma 2 the second. At most 3(|h⁻| + q) pairs are excluded and each lies
in q teams, so |ℱ|·p_K ≥ q³p_K − K³(|h⁻|+q)/(2q). On the other side
P'(H'⁻ = h⁻) = Π p_{β_f}·(1 − p_K)^{q³−q−|h⁻|} ≥ Π p_{β_f}·exp(−q³p_K − q³p_K²), and
q³p_K² = K⁶/(36q). QED

**Consequence (*).** Under P', ξ'_N and H'⁻ are independent and ξ'_N has the law of ξ_N.
So for every measurable G ≥ 0, with ρ(h⁻,N) the factor in front of P' in Lemma 3,

  E[ G(ξ_N, H⁻) ] = E[ Σ_{h⁻} G(ξ_N, h⁻)·P(H⁻ = h⁻ | d_N) ] ≤ E'[ G(ξ'_N, H'⁻)·ρ(H'⁻, N) ].

This holds for all h⁻ with an explicit factor, so nothing has to be restricted to "h with
no pair-sharing", and it is a pointwise statement about the conditional law, which is
stronger than a comparison in expectation.

Why the planted form is needed. A cheap assignment has about 2% of its teams above K = 5,
and K cannot grow with q (Δ is of order K⁵q^(2/3)). The number of cheap assignments is
therefore not a function of H, and P(H = h) ≤ exp(o(q))·P(H' = h) alone does not bound
its expectation. In (*) the costs of all teams of N, above and below K, are carried along
exactly, with the same law in both models. No separate union bound for teams above K is
needed in the transfer.

## 6. Lemma 4: removing the weight of teams that touch N (complete)

Apply (*) to G = 1{Σ_{e∈N} ξ_e ≤ (e_1*−ε)q}·1{H ∈ 𝒢}·1{s(H) ≤ s_0}·1{|H| ≤ K³q}.
On the support of G, |h⁻| ≤ K³q and s(h⁻;N) ≤ s(h⁻) + |h_touch| ≤ s_0 + |h_touch|, so

  ρ ≤ exp(Δ + 2K⁶ + s_0·ln W)·W^X,  X = |H'_touch|

(the middle terms of Lemma 3 are K⁶/2 + K³/2 + K⁶/(36q) ≤ 2K⁶).

X is not small on any event of H alone: under P, given that a team of N is near-tight, each
of the q − 1 other teams through one of its pairs is near-tight with probability of order
q^(−4/3), not q^(−2). The inflation is real, and X is of order q^(2/3) for a cheap N.

**Lemma 4.** Let A_N = {Σ_{e∈N} ξ'_e ≤ (e_1*−ε)q}. For q ≥ q_0(K, bins, θ),

  E'[ 1_{A_N}·1{H' ∈ 𝒢}·W^X ] ≤ exp(K³W/2)·( E'[1_{A_N}·1{H' ∈ 𝒢^(2θ)}] + 3·P'(A_N)·exp(−(θ/4)·q·ln q) ).

Proof. Under P' the three objects ξ'_N, H'_touch (near-tight teams in 𝒯_N) and H'_free
(near-tight teams outside N ∪ 𝒯_N) are independent, since N and 𝒯_N are fixed families.
X is Binomial(3q(q−1), p_K), and H'_0 := H' − H'_touch is a function of (ξ'_N, H'_free).
If H' ∈ 𝒢 and X ≤ θq then H'_0 ∈ 𝒢^(θ); if H'_0 ∈ 𝒢^(θ) and X ≤ θq then H' ∈ 𝒢^(2θ). So

  E'[1_A 1{H'∈𝒢} W^X] ≤ E'[1_A 1{H'_0 ∈ 𝒢^(θ)}]·E'[W^X] + P'(A)·E'[W^X; X > θq],
  E'[1_A 1{H'_0 ∈ 𝒢^(θ)}] ≤ E'[1_A 1{H' ∈ 𝒢^(2θ)}] + P'(A)·P'(X > θq).

E'[W^X] ≤ (1 + W p_K)^{3q²} ≤ exp(K³W/2). For x_0 = θq,
E'[W^X; X > x_0] ≤ Σ_{x > x_0} (K³W/2)^x/x! ≤ 2·(e·K³·C1/(2θ·q^(1/3)))^{θq}
≤ exp(−(θ/4)·q·ln q) for large q, and P'(X > θq) is smaller. Since exp(K³W/2) ≥ 1 the
three error terms combine as stated. QED

## 7. Proposition 5: the good event in the pair-additive market (complete)

(a) E|H| = K³q/6 =: μ. Covariances of the indicators vanish except for sharing pairs, so
Var|H| ≤ μ + 2Δ and Chebyshev gives P(|H| > 6μ) ≤ (μ + 2Δ)/(25μ²) = O(1/q).

(b) s(H) = Σ_π (mult(π) − 1)^+ ≤ Σ_π C(mult(π), 2) = Y, the number of unordered sharing
pairs of teams inside H, and E Y = Δ exactly. Markov.

(c), (d) Fix an agent a of sector A; its teams are (a,b,c), (b,c) ∈ [q]². Let X_β count
those in bin β. For k = (k_β) with |k| = Σ k_β, the joint factorial moment
E Π_β (X_β)_{k_β} is the sum, over ordered |k|-tuples of distinct teams at a with assigned
bins, of the probability that each is in its bin. The BC pairs of the tuple are distinct,
so its excess is s = 2|k| − #{b_i} − #{c_i}.
- s = 0: there are (q)_{|k|}² such tuples, the teams are pair-disjoint, and each tuple has
  probability exactly Π p_{β_i}. Total: Π Λ_β^{k_β}·(1 + O(|k|²/q)).
- s ≥ 1: at most |k|^{2|k|}·q^{2|k|−s} tuples, each of probability at most W^s·Π p_{β_i} by
  Lemma 1 with N = ∅. Total at most |k|^{2|k|}·Π Λ_{β_i}·Σ_{s≥1} (W/q)^s = O(q^(−1/3)).
So all joint factorial moments converge to those of independent Poisson(Λ_β) variables,
a law determined by its moments; the binned star of a converges in distribution to π_K.

For two agents a ≠ a' of sector A the same computation applies to the mixed factorial
moments. A team at a and a team at a' share a pair only if they have the same (b,c).
Write the excess of a mixed tuple as s_1 + s_2 + s_3: s_1, s_2 the excesses of the two
single-agent tuples and s_3 the number of BC coincidences between them. If the s_3
coincident teams at a' use r values of b and c values of c, those r + c values are forced,
so the number of tuples is at most const·q^{2|k|−s_1}·q^{2|k'|−s_2−(r+c)}, and
r + c ≥ max(2, 2s_3 − s_2) when s_3 ≥ 1 (s_3 distinct teams on r + c values have excess
2s_3 − r − c ≤ s_2). With the factor W^{s_1+s_2+s_3} of Lemma 1 the net power of q is
−(s_1+s_2)/3 − (r+c) + 2s_3/3 ≤ −1/3 whenever the total excess is positive (it is at most
−4/3 when s_3 ≥ 1). Tuples of excess 0 number (q)_{|k|}²·(q)_{|k'|}²·(1 − O(q^(−2))) and
have product probability; the rest contribute O(q^(−1/3)). Hence (star_a, star_a') converges in distribution to π_K ⊗ π_K, with
all moments of the degrees bounded uniformly in q. For g with |g| ≤ c(1 + degree)^m this
gives E g(star_a) → E_π g and Cov(g(star_a), g(star_a')) → 0 (uniform integrability), so
Var((1/q)·Σ_{v∈A} g(star_v)) ≤ Var(g(star_a))/q + Cov(g(star_a), g(star_a')) → 0 by
exchangeability. Chebyshev gives (d). For (c) take g the indicator of a star type: the
empirical law converges pointwise in probability on a countable set to the probability
measure π_K; for a finite set S_0 of types with π_K(S_0) ≥ 1 − η', the total variation
distance is at most Σ_{S_0} |L(𝔰) − π(𝔰)| + L(S_0^c) + η' → at most 2η' in probability.
Sectors B and C are the same by symmetry. QED

This is all the local convergence the transfer needs: depth-1 stars, law of large numbers
over agents, no rate. The `localconv` line is not used.

## 8. Proof of Theorem T

By symmetry over N and (*), with the bound on ρ from section 6,

  E[Z_ε; H ∈ 𝒢, s(H) ≤ s_0, |H| ≤ K³q] = Σ_N E[G_N(ξ_N, H⁻)]
   ≤ exp(Δ + 2K⁶ + s_0 ln W)·Σ_N E'[1_{A_N}·1{H' ∈ 𝒢}·W^{X_N}],

and Lemma 4 bounds each term; Σ_N E'[1_{A_N} 1{H' ∈ 𝒢^(2θ)}] = E'[Z_ε; H' ∈ 𝒢^(2θ)] and
Σ_N P'(A_N) = E'[Z_ε]. The good event is a function of H alone, so
P(Z_ε ≥ 1) ≤ P(H not good) + E[Z_ε; H good]; exp(R_q)·3·E'[Z_ε]·exp(−(θ/4) q ln q) ≤
exp(−(θ/5) q ln q) for large q because R_q = o(q) and E'[Z_ε] = exp(O(q)). QED

Size of the lower-order factors, all exp(o(q)): Janson correction exp(Δ),
Δ ≤ 0.075·K⁵·q^(2/3); pair-sharings inside H, exp(s_0·ln W) = exp(O(q^(2/3+η) ln q));
teams sharing a pair with the candidate assignment, exp(K³W/2) = exp(O(q^(2/3)));
excluded teams in the void probability, exp(2K⁶). The number of profiles and the binning
error do not appear in the transfer.

## 9. Numerics (evidence only)

`cmp.py`, output `out_cmp.txt`: q = 16, 24, 32, 48 and K = 0.8, 1.0, 1.3, 2.0; 40,000,
40,000, 30,000 and 15,000 markets; pair-additive (PA) against independent teams (IT),
no bins. s is the pair-sharing excess, Y the number of unordered sharing pairs in H.

| Check | Prediction | Observed |
|---|---|---|
| E[Y] under PA | Δ = 3q³(q−1)t⁵/40 exactly | agrees in all 16 cells within sampling error (e.g. q=48, K=2: 31.19 against 31.04) |
| P_PA(H empty) | at most exp(−μ+Δ) (Lemma 2); at least the IT value (Harris) | holds in all 11 cells with nonzero counts; the bound is nearly attained (q=32, K=0.8: 0.0822 against 0.0829; q=24, K=1: 0.0300 against 0.0333; IT 0.0652 and 0.0183) |
| P_PA(|H|=n, s=0)/P_IT(same), max over n with at least 30 events on each side | at most exp(Δ + K³n/2q + K⁶/36q) (Lemma 3 with N = ∅) | 1.13 to 4.24 in the 12 cells with K ≤ 1.3. One cell exceeded its bound in the first run: q=32, K=1, n=0 gave 2.22 against exp(Δ) = 2.08, from about 290 and 130 events. `cmp2.py` (`out_cmp2.txt`) reran q=32, K=1 with 150,000 PA and 600,000 IT markets: the ratio is 1.71 ± 0.06 at n=0, falls with n, and is below the bound for every n from 0 to 12. The exceedance was sampling noise |
| P_PA(|H|=n, s=1)/P_IT(same), max over n | at most the above times the exact one-sharing factor 1.8·q^(2/3)/K | 6.0 to 19.6 for K ≤ 1.3; the one-sharing factor alone is 8.8 to 29.7; two cells (K=1.3, q=16 and 24) exceed the bare factor, by ratios 1.23 and 1.14, and are inside the full bound (exp(Δ) alone is 5.2 and 9.2 there) |
| P_PA(s=1)/P_IT(s=1), all n | of order q^(2/3) when Δ is small | 8.5, 13.3, 15.2, 17.4 at K=0.8, q=16..48 |
| star size of one agent | Poisson(K³/6) with error O(K⁵q^(−1/3)) | TV 0.0090, 0.0083, 0.0078, 0.0069 (K=0.8); 0.151, 0.144, 0.136, 0.126 (K=2), q=16..48 |
| covariance of 1{star empty} for two agents | tends to 0 | at most 0.0033 in absolute value in all cells |

Reading. Nothing contradicts Lemmas 1-3 or Proposition 5. The convergence in Proposition
5 is slow at K = 2 and will be slower at K = 5: Δ/μ = 0.45·K²·q^(−1/3), which at K = 5 is
below 0.1 only for q above about 10⁶. The transfer is an asymptotic statement and says
nothing about the constants seen at q ≤ 160; this agrees with item 5 of the round-1
synthesis. For K = 2 the classes s ≤ 1 are too rare under PA to compare.

## 10. Remaining gaps and confidence

1. **(I-a*)**: the Poissonised first moment has to be delivered for the one-team law F,
   for Z_ε (all cheap assignments), on a typical event that contains the 2θq-neighbourhood
   of an event typical in the pair-additive market, with an exponentially small bound. I
   expect the star-conditioned computation to give this (its exponent is continuous in the
   star law), but I have not seen its write-up. In particular, if its typical event uses
   statistics other than the total-variation distance and polynomially bounded averages,
   Proposition 5(d) and the stability remark have to be re-checked for them.
2. Theorem T, Lemmas 1-4 and Proposition 5 have not been read by anyone else.
3. e_1* = 2.3025 is taken as given.

Confidence: Theorem T correct as stated, 90%. Target statement true, 85%. Target statement
provable this way once (I-a*) is written, 80%; the residual risk is in the exact form of
the typical event in (I-a), not in the transfer.
