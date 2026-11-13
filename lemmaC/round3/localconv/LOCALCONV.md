# Local weak convergence of the near-tight hypergraph (round 3, lane `localconv`)

Status labels: **complete** = I would defend it to a referee; **gap** = not proved here.
Floating point is evidence only. Nothing here has been checked by a second reader.
Scripts: `pairshare.py` (output `out_pairshare.txt`), `tail.py` (output `out_tail.txt`).

## 0. Setting and notation

Pair deficits δ_p (p a pair ab, bc or ac) are 3q² iid U[0,1]. Team deficit
d(e) = sum of its three pair deficits, ξ_e = q^(2/3)·d(e). Fix K and put

    s = K·q^(−2/3),   λ = K³/6,   S_r = Σ_{j=0}^{r−1} (2λ)^j.

Throughout K ≤ q^(2/3), i.e. s ≤ 1, so that P(d(e) ≤ t) = t³/6 and
P(δ+δ' ≤ t) = t²/2 hold exactly for t ≤ s.

- H_q(K): vertices = the 3q agents (each carries its sector label), edges = teams with
  ξ_e ≤ K, marked by ξ_e.
- Distance between agents = least number of teams in a chain joining them. B_r(o) = the
  rooted marked hypergraph formed by all edges of H_q(K) containing an agent at distance
  ≤ r−1 from o, with their agents and sector labels, up to root-preserving isomorphism.
  Agents at distance ≤ r−1 are *interior*.
- T_K: the Poisson hypertree. Every agent has child teams whose marks form a Poisson
  process of intensity ξ²/2 on [0,K] (so Poisson(λ) many), independently over agents;
  each child team has two new agents, one in each of the two other sectors. ν = ν_{K,r,X}
  is the law of B_r(root) when the root has sector X. μ = μ_{q,K,r,X} is the law of B_r(o)
  for an agent o of sector X (it does not depend on which o).
- For a finite tree t: |E| = number of teams, V = 1+2|E| = number of agents, m = number
  of interior agents, P2 = number of unordered pairs of teams sharing an agent.

A marked tree with distinct marks has a canonical form: order the child teams of each
interior agent by increasing mark; the two children of a team are told apart by sector.
With respect to (counting measure on shapes) × (Lebesgue measure on the ordered marks),
ν has density

    g(t,ξ) = e^(−mλ) · Π_{e∈t} ξ_e²/2.                                   (0.1)

(Ordered Poisson points: no factorial.) Under ν, E|E| = λ·S_r and E m = S_r.

## 1. One-team facts (all **complete**, elementary integrals)

(1a) Conditional on ξ_e ≤ K the three pair deficits (x₁,x₂,x₃) of e are uniform on the
solid simplex {x_i ≥ 0, Σx_i ≤ s}; the density of x₁ is 3(s−x)²/s³.

(1b) A team α sharing exactly the pair p with e is near-tight with conditional probability
(s−x_p)²/2 given e's deficits, and E[(s−x_p)²/2 | e near-tight] = 3s²/10. Given e's
deficits, the 3(q−1) teams sharing a pair with e are independent (their other pairs are
pairwise distinct).

(1c) **Pair-sharing rate.** For a near-tight team e,

    E[# near-tight teams sharing a pair with e | e near-tight]
        = 3(q−1)·3s²/10 = (9/10)·K²·q^(−1/3)·(1 − 1/q).

    P(e has no such partner | e near-tight)
        = E Π_{i=1..3} (1 − (s−x_i)²/2)^(q−1),  x uniform on the solid simplex.

So the working guess is right, with constant 9/10: pair-sharing per near-tight team is of
order K²q^(−1/3), and the hypertree picture needs q^(1/3) ≫ K².

(1d) Two teams sharing a pair are both near-tight with probability
∫₀^s ((s−x)²/2)² dx = s⁵/20. Hence for a fixed agent o, the number X_o of unordered pairs
of near-tight teams through o sharing a pair that contains o has

    E X_o = q²(q−1)·s⁵/20 = (K⁵/20)·q^(−1/3)·(1 − 1/q).

(The sister report's q³s⁵/40 is this count for one choice of the second agent's sector;
both sectors together give 1/20. Confirmed numerically in section 6.)

## 2. Theorem A: one root (**complete**)

**Theorem A.** For K ≤ q^(2/3), r ≥ 1 and every root sector X,

    d_TV(μ_{q,K,r}, ν_{K,r}) ≤ ε₁(q,K,r)
      := (9/10)·K²·q^(−1/3)·λ·S_r            [pair-sharing with a tree team]
       + (2λ+λ²)·S_r·K·q^(−2/3)              [a team closing a triangle on two tree teams]
       + Q_r/q                               [collisions of agents in the embedding]
       + λ²·S_r/q²,                          [Bernoulli versus Poisson]

where Q_r = E_ν V² ≤ (r+1)²(2r+1)·max(1,2λ)^(2r). In particular
P(B_r(o) is not a hypertree, or two of its teams share a pair) ≤ ε₁.

Whenever the bound is below 1 the first term dominates, and it equals
(9/10)K²q^(−1/3) × (expected number of teams in the limiting r-ball). For K ≥ 2 it is at
most 0.72·K²·(K³/3)^r·q^(−1/3). So radius-r convergence starts when
q^(1/3) ≫ K^(3r+2)/3^r; per team the requirement is q^(1/3) ≫ K².

**Proof.**

*Step 1 (exact density of μ on trees).* Fix a canonical marked tree (t,ξ) with distinct
marks and root sector X. Let φ range over injective sector-respecting maps of the agents
of t into the market with φ(root) = o. Then

    {B_r(o) ≅ (t,ξ)} = disjoint union over φ of
        {φ(e) has mark ξ_e for every e ∈ t} ∩ {W_φ = 0},

where W_φ counts near-tight teams that are not tree teams and contain an interior agent
of φ(t). (If the tree teams are near-tight and no other near-tight team meets the
interior, induction on distance shows the agents at distance j are exactly the depth-j
agents of φ(t). Disjointness: a mark- and sector-preserving automorphism fixing the root
is the identity when marks are distinct.)

Tree teams pairwise share at most one agent, so their 3|E| pairs are distinct and their
marks are independent, each with density ξ²/(2q²) on [0,K]. Writing x for the deficits of
all tree pairs, the density of μ at (t,ξ) is

    h(t,ξ) = Σ_φ Π_e ξ_e²/(2q²) · E[ P(W_φ = 0 | x) | marks = ξ ].         (2.1)

The number of φ is a product of 2|E| factors, each at least q − V, so it is at least
q^(2|E|)·(1 − 2|E|V/q)₊.

*Step 2 (classification of the other teams).* Let F be the set of non-tree teams meeting
the interior of φ(t), split by the number of tree pairs they contain.

- F0 (none): |F0| ≤ m·q². Given x each is near-tight with probability s³/6 = λ/q².
- F1 (exactly one, the pair p of tree team e): contained in the (q−1) teams through p.
  Given x it is near-tight with probability (s−x_p)²/2.
- F2 (exactly two): the two pairs share an agent b, so they lie in two tree teams e₁,e₂
  through b, and the team is one of the two "crossed" teams of (e₁,e₂). (Two pairs of the
  same tree team would make it that tree team.) So |F2| ≤ 2·P2, and given x each is
  near-tight with probability ≤ s.
- Three tree pairs: impossible. They would come from three distinct tree teams forming a
  cycle in the incidence graph of a tree.

*Step 3 (Harris).* Conditionally on x, the remaining pair deficits are iid U[0,1], and
each event {α not near-tight}, α ∈ F, is increasing in them. By Harris' inequality

    P(W_φ = 0 | x) ≥ Π_{α∈F} (1 − p_α(x))
       ≥ (1 − λ/q²)^(m q²) · (1 − Σ_{F1} p_α(x) − 2·P2·s)₊
       ≥ e^(−mλ)·(1 − mλ²/q²)·(1 − (q−1)·Σ_{e∈t} Σ_{p⊂e} (s−x_p)²/2 − 2·P2·s)₊ ,

using log(1−u) ≥ −u−u² for u = λ/q² ≤ 1/6.

*Step 4 (density ratio).* Dividing (2.1) by (0.1), and using
(1−a)₊(1−b)₊(1−c)₊ ≥ 1−a−b−c for a,b,c ≥ 0, the ratio f = h/g satisfies

    (1 − f(t,ξ))₊ ≤ 2|E|V/q + mλ²/q² + 2·P2·K·q^(−2/3)
                    + (q−1)·Σ_e Σ_{p⊂e} E[(s−x_p)²/2 | ξ_e].               (2.2)

*Step 5 (total variation needs only the lower bound).* μ restricted to trees is f·ν and
the rest of μ sits on non-trees, where ν has no mass. Hence
d_TV(μ,ν) = ∫ (1−f)₊ dν. Integrate (2.2) against ν:

- by the Mecke formula and (1a)–(1b), the marks of the teams of the ν-tree have density
  3ξ²/K³, so E_ν of the last term is E_ν|E|·(q−1)·3·(3s²/10) ≤ (9/10)K²q^(−1/3)·λS_r;
- E_ν P2 = λ²/2 + (λ + λ²/2)(S_r − 1) ≤ (λ + λ²/2)·S_r (a root with N child teams
  contributes C(N,2), any other interior agent C(N+1,2));
- 2|E|V ≤ V², and with Z_j the number of agents at depth j, E Z_j² = (2λ)^(2j)
  + 2Σ_{i=1..j}(2λ)^(2j−i) ≤ (1+2j)·max(1,2λ)^(2j), which gives the stated bound on Q_r;
- E_ν m = S_r.

This is ε₁. Finally μ(non-trees) = 1 − ∫ f dν ≤ ∫(1−f)₊ dν. ∎

Remarks.
- Forgetting sector labels can only decrease total variation, so the statement holds for
  unlabelled neighbourhoods and for a uniformly random root.
- *Conditioning through pairs.* The only place where "a near-tight team has small pair
  deficits, which helps the teams through its pairs" enters is the F1 and F2 terms. It is
  exactly this effect that produces the K²q^(−1/3) rate; teams sharing only an agent with
  the tree contribute nothing beyond the Poisson term.
- *Why no upper bound on f is needed.* An upper bound would need a Janson-type inequality
  for the clumped family F0. I believe the Boppana–Spencer argument goes through after
  conditioning on the deficits of one team, giving P(W=0) ≤ exp(−EW + Δ), but I have not
  written it out, and total variation does not use it.
- *Sharpness of the K dependence.* d_TV ≥ μ(non-trees) ≥ P(X_o ≥ 1), and E X_o
  = K⁵q^(−1/3)(1−1/q)/20 exactly. For r = 1 the upper bound is 3K⁵q^(−1/3)/20 (and
  K⁵q^(−1/3)/10 if one notes that only the two pairs through o matter at the last level).
  The matching lower bound P(X_o ≥ 1) ≥ E X_o − E C(X_o,2), with
  E C(X_o,2) = O(K^10 q^(−2/3) + K⁷q^(−2/3) + K⁹/q), is a routine count that I have only
  sketched: **gap (minor)**. The numerics in section 6 agree with the first moments.

## 3. Theorem B: two roots, and convergence in probability (**complete**)

**Theorem B.** For distinct agents o ≠ o' with sectors X, X', the joint law μ₂ of
(B_r(o), B_r(o')) satisfies d_TV(μ₂, ν_X ⊗ ν_X') ≤ 4·ε₁(q,K,r).

Proof. Repeat Steps 1–5 for a forest (t,ξ) ⊔ (t',ξ') and pairs (φ,φ') with disjoint
images. Disjoint embeddings give a lower bound for the density of μ₂ (overlapping balls
only add mass). The number of such pairs is at least q^(2|E|+2|E'|)·(1 − 2(|E|+|E'|)(V+V')/q)₊.
Tree pairs of the two images are distinct, so the marks are independent. F0 has at most
(m+m')q² teams; F1 is as before for each tree; a team containing two forest pairs
contains their common agent, so both pairs lie in the same tree and |F2| ≤ 2(P2+P2').
Therefore (1−f₂)₊ is at most the sum of the two one-tree bounds (2.2) with the collision
term replaced by 2(|E|+|E'|)(V+V')/q ≤ (V+V')²/q, whose ν⊗ν-mean is at most 4Q_r/q. ∎

**Corollary B (local weak convergence in probability).** Let h be a measurable function
of r-balls with |h| ≤ 1, and A_q = (3q)^(−1) Σ_v h(B_r(v)). Let ν̄ be ν with a uniformly
random root sector. Then

    |E A_q − ν̄(h)| ≤ 2ε₁,   Var(A_q) ≤ 1/(3q) + 12·ε₁,

so A_q → ν̄(h) in probability, with rate q^(−1/3) at fixed K, r. For the number N of
agents with a given neighbourhood type (h an indicator), Var N ≤ 3q + 54·q²·ε₁
= O(C(K,r)·q^(5/3)) = o(q²).

Proof. For o ≠ o', |E h(B(o))h(B(o')) − ν_X(h)ν_X'(h)| ≤ 2·4ε₁ and
|E h(B(o))·E h(B(o')) − ν_X(h)ν_X'(h)| ≤ 4ε₁ (signed functions cost a factor 2 over total
variation; for indicators the two bounds are 4ε₁ and 2ε₁). So |Cov| ≤ 12ε₁ (6ε₁ for
indicators) for off-diagonal pairs, and the 3q diagonal terms are at most 1 each. ∎

## 4. Uniform integrability (**complete given fact 3**; I did not open FPT)

Input: g(k) := E(3k − T_k) ≤ C₀·k^(1/3) for all k ≥ 1 (fact 3 of the brief with a constant
uniform in k; for large k the brief quotes C₀ ≈ 4.87). I use it as stated there.

**Theorem C (integer optimum).** For K ≥ 64 and all q,

    E[(1/q) Σ_{e∈M} (ξ_e − K)₊] ≤ (2C₀ + 1.26·C₀^(4/3))·K^(−1/6).

Equivalently, for a uniformly random agent o, sup_q E[(ξ_{M(o)} − K)₊] → 0 as K → ∞, so
{ξ_{M(o)}} is uniformly integrable (E[ξ; ξ > 2K] ≤ 2E(ξ−K)₊).

Proof (sprinkling). Fix η ∈ (0,1/2]. Realise each pair deficit as δ = min(δ',δ'') with
δ', δ'' independent, P(δ' > t) = (1−t)^(1−η), P(δ'' > t) = (1−t)^η. Put κ = Kq^(−2/3) and,
for a cost vector c, Φ(c) = min over partial assignments N of
[Σ_{e∈N} d_c(e) + κ·(q − |N|)].

1. Φ(δ) ≤ Σ_{e∈M} min(d(e),κ) = (3q − T) − Σ_{e∈M}(d(e) − κ)₊ (keep M's teams with d ≤ κ).
2. Let N' attain Φ(δ'), S its uncovered agents (k per sector). Completing N' by the
   δ''-optimal assignment on S gives a perfect assignment, and δ ≤ δ', δ ≤ δ'':
   3q − T ≤ Φ(δ') − κk + T''(S).
3. S is a function of δ', independent of δ''. By Bernoulli's inequality
   δ'' = 1 − (1−U)^(1/η) ≤ U/η for a uniform U, so E[T''(S) | δ'] ≤ g(k)/η ≤ C₀k^(1/3)/η.
4. Likewise δ' ≤ U/(1−η), and Φ is monotone and Φ(c/(1−η)) ≤ Φ(c)/(1−η), so
   E Φ(δ') ≤ E Φ(δ)/(1−η).
5. κ·E k ≤ E Φ(δ') and E Φ(δ) ≤ E(3q−T) ≤ C₀q^(1/3). By Jensen,
   E k^(1/3) ≤ (C₀q/((1−η)K))^(1/3).

Combining, E Σ_{e∈M}(d(e)−κ)₊ ≤ E(3q−T) − EΦ(δ)
≤ C₀q^(1/3)·η/(1−η) + (C₀^(4/3)/η)·q^(1/3)·((1−η)K)^(−1/3). Multiply by q^(−1/3) and take
η = K^(−1/6). ∎

**Theorem C′ (LP optimum).** The same bound, with C₀^(4/3) replaced by
C₀·(C₀+1)^(1/3) for K ≤ q, holds for (1/q)Σ_e x*_e(ξ_e − K)₊ where x* is the LP optimum.
Proof: use the fractional Φ (x ≥ 0, load ≤ 1 at every agent, penalty κ(q − Σx_e)). In
step 2 the deficiency vector u of the fractional optimum has equal sector sums k; shrink
x by a factor (1−θ) to make k an integer ⌈k⌉, write each sector's u as a convex
combination of indicators of ⌈k⌉-subsets, and patch with the corresponding mixture of
δ''-optimal assignments. LP cost ≤ IP cost gives E(3q − LP*) ≤ C₀q^(1/3). ∎

**Corollary (repair lemma, the form a transfer from the tree needs).** If x_q is any
(fractional or integral) sub-assignment of the market with expected rescaled cost
(1/q)E Σ ξ_e x_e ≤ a and expected uncovered mass ≤ εq per sector, then for every
η ∈ (0,1): q^(−1/3)·E(3q − LP*) ≤ a/(1−η) + C₀(ε+1/q)^(1/3)/η, and the same for
E(3q − T) when x_q is integral. (Apply x_q's rule to the δ' market, rescaled; steps 2–4.)

What each direction of a cost transfer needs:
- *Finite → tree (lower bounds on a_IP, a_LP through the tree).* Only tightness, which is
  Markov: P(ξ_{M(o)} > K) ≤ C₀/K, and likewise the x*-mass on teams with ξ > K at a
  random agent is ≤ C₀/K. Costs are lower semicontinuous, so limit objects cost at most
  the limit of the finite costs. **complete** (given Theorems A, B and fact 3).
- *Limit of the optimum has exactly the limiting cost.* Theorems C, C′.
- *Tree → finite (upper bound on a_LP from a tree fractional matching).* The repair
  corollary, **plus** a way to realise the tree object on H_q(K) with small deficiency.
  If the tree object is a factor of iid marks this follows from Theorem A/B; if it is only
  unimodular it is the usual soficity-type question: **gap**, outside this lemma.

Caveats: the constant 20·K^(−1/6) (with C₀ = 4.87) is far from the truth; in simulations
the optimal assignment has no team above ξ = 8 at q ≤ 32 (section 6). The rate is not
needed for uniform integrability.

## 5. What is not proved here (**gaps**)

1. Matching lower bound for the total variation (second-moment count for X_o), section 2.
2. Fact 3 with a uniform constant is taken from the brief; FPT not opened by me.
3. Realising a non-factor unimodular tree object on the finite market (section 4).
4. Convergence is for fixed K and r. Nothing here controls K growing with q; by (1c) the
   statement is false once K² is comparable to q^(1/3).

## 6. Numerics (floating point; evidence only)

`pairshare.py`: fresh markets, all q³ teams enumerated. "formula" columns: K³/6; (1c)
mean partners; (1c) probability of at least one partner (Monte Carlo integral over the
simplex); (1d) E X_o.

| q | K | reps | degree sim / K³/6 | partners per team sim / formula | P(≥1 partner) sim / formula | X_o per agent sim / formula | P(1-ball not a tree) sim |
|---|---|---|---|---|---|---|---|
| 32 | 2 | 300 | 1.319 / 1.333 | 1.054 / 1.099 | 0.650 / 0.665 | 0.459 / 0.488 | 0.235 |
| 32 | 4 | 300 | 10.70 / 10.67 | 4.409 / 4.394 | 0.987 / 0.987 | 15.73 / 15.62 | 0.952 |
| 32 | 6 | 300 | 35.73 / 36.00 | 9.820 / 9.886 | 1.000 / 1.000 | 117.0 / 118.6 | 1.000 |
| 64 | 2 | 80 | 1.312 / 1.333 | 0.872 / 0.886 | 0.580 / 0.585 | 0.385 / 0.394 | 0.210 |
| 64 | 4 | 80 | 10.86 / 10.67 | 3.595 / 3.544 | 0.970 / 0.968 | 13.00 / 12.60 | 0.952 |
| 64 | 6 | 80 | 35.94 / 36.00 | 7.954 / 7.973 | 1.000 / 1.000 | 95.08 / 95.68 | 1.000 |
| 128 | 2 | 14 | 1.339 / 1.333 | 0.682 / 0.709 | 0.487 / 0.505 | 0.303 / 0.315 | 0.185 |
| 128 | 4 | 14 | 10.53 / 10.67 | 2.812 / 2.835 | 0.934 / 0.937 | 9.915 / 10.08 | 0.945 |
| 128 | 6 | 14 | 35.87 / 36.00 | 6.361 / 6.379 | 0.998 / 0.998 | 75.99 / 76.55 | 0.999 |

(q = 48 and 96 are in `out_pairshare.txt`; same agreement.) The formulas hold. The
consequence is the important part: at q ≤ 128 the hypertree regime is not visible. With
K = 4 more than 93% of near-tight teams have a near-tight partner through a shared pair,
and 95% of agents have a non-tree 1-ball. Even at K = 2 about 19–23% of 1-balls are not
trees. Per-team pair-sharing below 10% needs q ≳ 3·10⁶ at K = 4 and q ≳ 3·10⁷ at K = 6.

`tail.py` (optimal assignment M and LP optimum x*, instances from `scaling/cg.py`):

| q | reps | mean ξ over M | mean max ξ in M | E[ξ; ξ>K] per team, K = 2, 4, 6, 8 (M) | same for x* |
|---|---|---|---|---|---|
| 16 | 60 | 2.643 | 5.04 | 2.183, 0.668, 0.041, 0.000 | 2.164, 0.617, 0.033, 0.000 |
| 24 | 60 | 2.576 | 5.17 | 2.067, 0.542, 0.027, 0.000 | 2.047, 0.498, 0.029, 0.001 |
| 32 | 40 | 2.586 | 5.43 | 2.090, 0.627, 0.031, 0.000 | 2.057, 0.546, 0.040, 0.000 |

(The q = 48 row was still running when this was written; if it finishes it is in
`out_tail.txt`.) The optimum lives on ξ ≤ 6 (about 1% of its cost lies above 6), so K ≈ 6 is the level a
transfer needs, and at that level pair-sharing per team is 0.9·36·q^(−1/3) ≥ 1 until
q ≈ 3·10⁴. This is consistent with reading (ii) of the sister report: the pair-additive
model can share the local limit of the independent-team model and still differ from it by
order-one amounts at every simulated size.
