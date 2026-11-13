# Theorem 1: the gap D = LP* − T is asymptotically deterministic

Full write-out, 12 November 2026. It replaces `reports/pairdeficit.md` §§1–5 as the proof of
Theorem 1 of `MAIN_THEOREM.md`. The mathematics is that of the first version; what is new is
the bookkeeping the two reviewers asked for (listed in §10). Written by one agent; not read
by a person. Every step below is claimed **complete** unless it is labelled **GAP**; no
step is so labelled. Caveats on the imported statement and on checking are in §11.

Logarithms are natural. "Pair" always means two agents from different sectors.

---

## 0. Setting and notation

Three sectors (1, 2, 3) of q agents each; a, b, c always denote agents of sectors 1, 2, 3.
The 3q² pair values are independent U[0,1]. The **deficit** of a pair is d = 1 − value, so
the 3q² pair deficits are independent U[0,1]. A team e = (a, b, c) has value
v(e) = 3 − d(e), where d(e) = d(a,b) + d(b,c) + d(a,c) is its deficit.

- A (perfect) **assignment** is a set of q disjoint teams. Its **cost** is the sum of the
  deficits of its 3q pairs. T = 3q − (least cost of an assignment). An assignment is
  determined by its set of pairs. M denotes an optimal assignment; it is almost surely
  unique (two distinct assignments have distinct pair sets, so a tie is a nontrivial linear
  equation in the deficits).
- P is the polytope of fractional assignments: x >= 0 on the q³ teams with
  sum over teams containing i of x_e = 1 for each agent i. LP* = max over P of
  sum_e v(e)·x_e. The **support** of x is {e : x_e > 0}. D = LP* − T >= 0.
- The **load** of a pair i under x is l_i(x) = sum of x_e over the teams e containing both
  agents of i. Then 0 <= l_i <= 1, sum_i l_i = 3q, sum_e v(e)x_e = sum_i l_i·(1 − d_i), and
  sum_i l_i·d_i = 3q − sum_e v(e)x_e.
- A **position** of type {i, j} is a pair (x, y) with x in sector i and y in sector j; its
  **third sector** is the remaining one. For 0 < s <= 1, an agent z of the third sector is
  **s-close** to the position (x, y) if d(x,z) + d(y,z) < s. This has probability
  p := s²/2 <= 1/2 (the sum of two independent uniforms is below s).
- Positions are **vertex-disjoint** if no two share an agent. For a set S of positions of
  one type, N(S) is the set of third-sector agents that are s-close to at least one
  position of S. N is monotone: S' ⊆ S implies N(S') ⊆ N(S).

---

## 1. Statement

**Theorem 1.** Fix eps in (0, 1/2) and put

    s = q^(−1/2+eps),   C_eps = 3·(ceil(1/(2 eps)) + 2).

(i) For every K >= 1 there is q_0(eps, K) such that for q >= q_0, with probability at least
1 − 150·2^K·q^(−K): every optimal assignment uses no pair of deficit >= C_eps·s, and every
team in the support of every optimal solution of the LP has deficit < C_eps·s.

(ii) Let C_0 := sup_n E(3n − T_n)/n^(1/3), which is finite by Import F (§9). For
q >= q_0(eps, 2),

    Var T   <=  4·C_eps·C_0·q^(−1/6+eps) + 7000/q,
    Var LP* <=  4·C_eps·C_0·q^(−1/6+eps) + 7000/q,
    Var D   <= 16·C_eps·C_0·q^(−1/6+eps) + 28000/q.

(iii) If E D >= q^(−1/12+eps) for all large q, then
P(core nonempty) = P(D = 0) <= Var D/(E D)² = O(q^(−eps)) → 0. In particular
P(D = 0) → 0 whenever liminf E D > 0.

For eps >= 1/2 statement (i) is trivial (pair deficits are at most 1, team deficits at most
3), and (ii), (iii) for eps >= 1/2 follow from the case eps = 1/4.

The explicit thresholds: q_0(eps, K) is any q_0 such that every q >= q_0 satisfies

    (Q1)  q^(2 eps) >= 4·max(13(K+5), 600)·log q,
    (Q2)  q^(2 eps) >= (20/3)^ceil(1/(2 eps)),
    (Q3)  the exponentially small terms displayed in §8 are together below q^(−K).

These are astronomically large (for eps = 1/4, (Q1) alone needs q > 10^9); the theorem is
asymptotic only.

---

## 2. Elementary tools

**(C1) Chernoff, lower tail.** If Y ~ Bin(n, pi) and 0 <= delta <= 1, then
P(Y <= (1 − delta)·n·pi) <= exp(−delta²·n·pi/2).

**(C2) Chernoff–Hoeffding.** If Y ~ Bin(n, pi) and 0 < a < pi, then
P(Y <= a·n) <= exp(−n·KL(a ‖ pi)), where
KL(a ‖ pi) = a·log(a/pi) + (1 − a)·log((1 − a)/(1 − pi)).

**(C3) Monotonicity.** If pi >= pi_0 then P(Bin(n, pi) <= t) <= P(Bin(n, pi_0) <= t).

**(M) Bipartite matchings.** In a bipartite graph, if y >= 0 on the edges has
sum of y over the edges at each vertex <= 1, then there is a matching with at least
sum_edges y edges, hence (the size being an integer) at least ceil(sum_edges y) edges.
(König's theorem, or integrality of the bipartite matching polytope.)

**(LP) Duality.** The LP of §0 is feasible and bounded, so its dual
min sum_i u_i subject to u(e) := u_a + u_b + u_c >= v(e) for every team e = (a,b,c),
u unrestricted in sign, has an optimal solution u, and for every optimal x and every
optimal u: x_e > 0 implies u(e) = v(e) (complementary slackness).

**(ES) Efron–Stein.** If f is a square-integrable function of independent coordinates
E_1, ..., E_m, then Var f <= sum_i E[Var_i f], where Var_i is the variance in E_i with the
other coordinates held fixed.

**Lemma P (Poincaré inequality for product Exp(1), constant 4).** Let mu be the law of m
independent Exp(1) variables on [0, ∞)^m. Let f : [0, ∞)^m → R be continuous, bounded and
Lipschitz in each coordinate separately. Then each partial derivative ∂_i f exists mu-almost everywhere
and

    Var_mu(f) <= 4·E_mu[ sum_i (∂_i f)² ].

*Proof.* One dimension first. Let g : [0, ∞) → R be bounded and Lipschitz; it is absolutely
continuous, g' exists almost everywhere and is bounded. Put h = g − g(0). Integrating by
parts on [0, R] (h² is absolutely continuous, h(0) = 0),

    ∫_0^R h²·e^(−x) dx = −h(R)²·e^(−R) + 2·∫_0^R h·h'·e^(−x) dx.

Let R → ∞; h is bounded, so the boundary term vanishes, and Cauchy–Schwarz gives
E h² = 2·E[h·h'] <= 2·(E h²)^(1/2)·(E h'²)^(1/2), that is E h² <= 4·E h'². Hence
Var g <= E(g − g(0))² <= 4·E g'².

For m coordinates: for every value of the other coordinates, t ↦ f(.., t, ..) is bounded
and Lipschitz, so Var_i f <= 4·E_i (∂_i f)² by the one-dimensional case, where the
derivative exists for almost every t. The set where ∂_i f fails to exist is measurable and
meets every line parallel to the i-th axis in a null set, so it is mu-null by Fubini. Take
expectations, sum over i and apply (ES). ∎

(The inequality with constant 4 is classical; `MAIN_THEOREM.md` cites Bobkov–Ledoux 1997.
The proof above makes this document independent of that citation.)

---

## 3. Lemma U (uniform expansion)

Consider an instance with sectors of size n and independent U[0,1] pair deficits.

**Lemma U.** Let n >= 2, K >= 1, 0 < s <= 1, p = s²/2, and assume

    (U-i)  n·p >= 13·(K+2)·log n,        (U-ii)  n·p >= 600·log n.

Let U(n, s) be the event that, for each of the three position types and **every** nonempty
set S of vertex-disjoint positions of that type,

    (a) |N(S)| >= 0.3·n·p·|S|   if |S| <= 1/p;
    (b) |N(S)| >= 0.9·n         if |S| >= 3/p.

Then P(U(n, s) fails) <= 6·n^(−K) + 3·exp(−0.009·n).

*Proof.* Fix a type, say {1, 3} with third sector 2, and a set S = {(x_1,y_1), ..., (x_m,y_m)}
of m vertex-disjoint positions. A sector-2 agent z belongs to N(S) iff
d(x_i,z) + d(y_i,z) < s for some i. Because the positions are vertex-disjoint, the 2m pairs
(x_i,z), (y_i,z) are distinct for fixed z, and pairs with different z are distinct. So the
indicators 1{z ∈ N(S)}, over the n agents z, are independent with probability
pi_m := 1 − (1−p)^m, and |N(S)| ~ Bin(n, pi_m).

The number of sets of m vertex-disjoint positions of one type is
binom(n,m)²·m! <= n^(2m).

(a) Let 1 <= m <= 1/p. Since 1 − p <= e^(−p) and t ↦ (1 − e^(−t))/t is decreasing on (0, 1],

    pi_m >= 1 − e^(−pm) >= (1 − e^(−1))·p·m =: pi_0,   1 − e^(−1) = 0.63212...

By (C3) and (C1) with (1 − delta)·(1 − e^(−1)) = 0.3, i.e. delta = 0.52540...,

    P(|N(S)| < 0.3·n·p·m) <= exp(−delta²·(1 − e^(−1))·n·p·m/2)
                           = exp(−0.08724..·n·p·m) <= exp(−0.079·n·p·m).

By (U-i), 0.079·n·p >= 1.027·(K+2)·log n >= (K+2)·log n. The union bound over all S with
1 <= |S| <= 1/p gives at most

    sum_{m>=1} n^(2m)·n^(−(K+2)m) = sum_{m>=1} n^(−Km) <= 2·n^(−K)     (n^K >= 2).

(b) By monotonicity of N it is enough to treat |S| = m_b := ceil(3/p) (a set with
|S| >= 3/p has at least m_b elements; if no such set exists, (b) is vacuous). Then
pi_{m_b} >= 1 − e^(−p·m_b) >= 1 − e^(−3) = 0.95021..., and by (C3), (C2),

    P(|N(S)| < 0.9·n) <= exp(−n·KL(0.9 ‖ 1 − e^(−3))) = exp(−0.020879..·n) <= exp(−0.0208·n).

The number of sets is at most n^(2 m_b), and 2·m_b·log n <= (6/p + 2)·log n. By (U-ii),
(6/p)·log n <= n/100, and (using p <= 1/2) log n <= n·p/600 <= n/1200, so
2·log n <= n/600. The union bound gives at most

    exp((0.01 + 0.001667 − 0.0208)·n) <= exp(−0.009·n).

Add (a) and (b) and multiply by 3 for the three types. ∎

The point of the lemma is that it holds for all S simultaneously, so S may be chosen after
looking at the instance (at M, or at an LP optimum). No conditioning is involved.

---

## 4. Lemma G (growth through the cutoff; deterministic)

Lemma U(a) stops at |S| = 1/p. This lemma carries the growth on to linear size.

**Lemma G.** Assume the event U(n, s) of Lemma U, with its hypotheses, and put
lam := 0.3·n·p. Then lam > 1, and for every nonempty set S of vertex-disjoint positions of
one type, with m = |S|:

    (G1) |N(S)| >= min(lam·m, 0.15·n);
    (G2) if m >= 0.15·n, then |N(S)| >= 0.9·n.

Moreover:

    (G3) if reals x_0, x_1, ... satisfy x_{i+1} >= min(lam·x_i, 0.15·n) for all i, then
         x_i >= min(lam^i·x_0, 0.15·n) for all i.

*Proof.* By (U-ii), n·p >= 600·log 2 > 415, so lam > 1.

(G1) If m <= 1/p this is U(a). If m > 1/p, let m_0 = floor(1/p); then 1 <= m_0 <= m
(p <= 1/2 gives m_0 >= 2), so S has a subset S' with exactly m_0 positions. S' is a set of
vertex-disjoint positions with |S'| <= 1/p, so by monotonicity and U(a),

    |N(S)| >= |N(S')| >= 0.3·n·p·floor(1/p) >= 0.3·n·p·(1/p − 1) = 0.3·n·(1 − p) >= 0.15·n,

the last step because p <= 1/2. (For small p the bound is nearly 0.3·n; 0.15·n is all that
is used.)

(G2) n·p >= 20 gives 0.15·n >= 3/p, so m >= 3/p and U(b) applies.

(G3) Induction on i. The case i = 0 is trivial. If x_i >= min(lam^i·x_0, 0.15·n), then

    x_{i+1} >= min(lam·x_i, 0.15·n) >= min(lam^{i+1}·x_0, 0.15·lam·n, 0.15·n)
             = min(lam^{i+1}·x_0, 0.15·n),

using lam >= 1. ∎

There is no rounding issue: (a), (G1), (G2) are real lower bounds on integer cardinalities,
and the only integer parts taken are floor(1/p) in (G1) and ceil(3/p) in Lemma U(b).

---

## 5. Theorem B (support of the LP)

**Theorem B.** Fix eps in (0, 1/2), K >= 1, s = q^(−1/2+eps), p = s²/2 = q^(−1+2 eps)/2 and

    r_B := ceil(1/(2 eps)) + 1.

Let q satisfy (Q1) and (Q2). Then, on the event U(q, s) of Lemma U (applied to the whole
instance, n = q), for every optimal solution x of the LP and every team f with x_f > 0,

    d(f) < 3·(r_B + 1)·s.

Consequently this holds with probability at least 1 − 6·q^(−K) − 3·exp(−0.009·q).

*Proof.* (Q1) gives n·p = q^(2 eps)/2 >= max(13(K+2), 600)·log q, so Lemma U applies with
n = q and Lemma G holds on U(q, s), with lam = 0.15·q^(2 eps).

Fix an optimal x, an optimal dual u as in (LP), and a support team f = (a, b, c). Recall
u(g) >= v(g) = 3 − d(g) for every team g, with equality on the support.

*Step 1 (start).* For every sector-2 agent b', u(a,b',c) >= v(a,b',c) and u(f) = v(f).
Subtracting,

    u_b' − u_b >= d(a,b) + d(b,c) − [d(a,b') + d(b',c)].

Let B_0 := N({(a, c)}). For b' ∈ B_0 the bracket is below s, so

    u_b' > u_b + d(a,b) + d(b,c) − s          for all b' ∈ B_0.                    (5.1)

*Step 2 (one link).* If g = (a', b', c') is a support team and b'' is s-close to the
position (a', c'), then u(a',b'',c') >= v(a',b'',c') and u(g) = v(g) give

    u_b'' − u_b' >= d(a',b') + d(b',c') − [d(a',b'') + d(b'',c')] > −s,              (5.2)

because deficits are nonnegative and the bracket is below s.

*Step 3 (growth).* Suppose B_i is a nonempty set of sector-2 agents. Let H_i be the
bipartite graph on (sector 1) × (sector 3) whose edges are the positions (a', c') for which
some support team (a', b', c') has b' ∈ B_i, and put
y(a', c') := sum over b' ∈ B_i of x(a', b', c'). Each team has exactly one sector-2 agent
and each agent has total weight 1, so sum over edges of y = |B_i|; and for each a',
sum_{c'} y(a', c') <= sum over teams containing a' of x = 1, and likewise at each c'. By (M)
H_i contains a matching S_i with |S_i| >= |B_i|. S_i is a set of vertex-disjoint positions of
type {1, 3}, each of which is the (sector 1, sector 3) part of a support team whose
sector-2 agent lies in B_i. Define

    B_{i+1} := N(S_i).

By (5.2), every b'' ∈ B_{i+1} satisfies u_b'' > u_b' − s for some b' ∈ B_i. With (5.1) and
induction,

    u_b'' > u_b + d(a,b) + d(b,c) − (i+1)·s      for all b'' ∈ B_i, all i >= 0.     (5.3)

*Step 4 (sizes).* U(a) with |S| = 1 gives |B_0| >= lam. By (G1) and |S_i| >= |B_i|,
|B_{i+1}| >= min(lam·|B_i|, 0.15·q), so by (G3)

    |B_i| >= min(lam^{i+1}, 0.15·q)      for all i >= 0.

Put h := ceil(1/(2 eps)) = r_B − 1. Then 2 eps·(h+1) >= 1 + 2 eps, so

    lam^{h+1} = 0.15^{h+1}·q^(2 eps (h+1)) >= 0.15^{h+1}·q^(1+2 eps) >= 0.15·q,

the last inequality being (Q2): q^(2 eps) >= (20/3)^h. Hence |B_h| >= 0.15·q, so
|S_h| >= 0.15·q, and (G2) gives

    |B_{r_B}| = |N(S_h)| >= 0.9·q.

(Count: h = ceil(1/(2 eps)) applications of (G1) after B_0, then one application of (G2);
total r_B = h + 1 steps after B_0, as in the first version.)

*Step 5 (symmetry).* The same argument with the roles of the sectors permuted (Lemma U
covers all three position types) gives sets A_r, B_r, C_r, one in each sector, each of size
>= 0.9·q, with r = r_B and

    u_a'' > u_a + d(a,b) + d(a,c) − (r+1)·s      for a'' ∈ A_r,
    u_b'' > u_b + d(a,b) + d(b,c) − (r+1)·s      for b'' ∈ B_r,
    u_c'' > u_c + d(a,c) + d(b,c) − (r+1)·s      for c'' ∈ C_r.

*Step 6 (close).* The support teams with a sector-1 agent outside A_r carry total weight
q − |A_r| <= 0.1·q; likewise for B_r and C_r. Since the total weight is q > 0.3·q, some
support team t = (a'', b'', c'') lies in A_r × B_r × C_r. Adding the three inequalities and
using u(t) = v(t) <= 3 and u(f) = v(f) = 3 − d(f),

    3 >= v(t) = u(t) > u(f) + 2·d(f) − 3(r+1)·s = 3 + d(f) − 3(r+1)·s,

so d(f) < 3(r_B + 1)·s. ∎

The argument uses complementary slackness, so it says nothing about M when D > 0; that is
Theorem A.

---

## 6. Theorem A (the optimal assignment)

**Theorem A.** Fix eps in (0, 1/2), K >= 1, s = q^(−1/2+eps), p = s²/2 and

    r_A := ceil(1/(2 eps)) − 1  (>= 1),     kappa := 0.3^{r_A}/(8·4^{r_A+1}).

Let q satisfy (Q1). Then the probability that some optimal assignment uses a pair of deficit
>= (r_A + 2)·s is at most

    eta_A := 3q³·[ 6·(q−1)^(−(K+3)) + 3·exp(−0.009·(q−1))
                   + exp(−q^(2 eps)/32) + exp(−kappa·q^(2 eps)) ].

Since (q−1)^(−(K+3)) <= 2^(K+3)·q^(−(K+3)), the first term is at most 144·2^K·q^(−K); the
other three are smaller than any power of q.

*Proof.* The law of the instance is invariant under permuting the three sectors, so it is
enough to bound the probability for pairs of type {1, 2} and multiply by 3. Fix e = (a, b).
Write delta := (r_A + 2)·s.

*Step 1 (threshold; deterministic).* Let W_0 be the least value of cost − d(e) over
assignments that use e, and W_1 the least cost over assignments that avoid e (both classes
are nonempty for q >= 2). Neither depends on d(e). The least cost of an assignment is
min(W_0 + d(e), W_1). If an optimal assignment uses e, then W_0 + d(e) <= W_1, that is
d(e) <= rho_e := W_1 − W_0. Hence

    {some optimal assignment uses e and d(e) >= delta} ⊆ {rho_e >= delta}.          (6.1)

*Step 2 (decoupling; deterministic).* For a sector-3 agent c let I_c be the instance with
a, b, c removed (sectors of size n := q − 1) and M_c an optimal assignment of I_c, chosen
by a fixed rule (say the lexicographically first optimum), so that M_c is a function of
the pairs of I_c only. An assignment that uses e contains the team (a, b, c) for exactly
one c, and the rest of it is an assignment of I_c. Hence

    W_0 = min over c of [ d(a,c) + d(b,c) + cost(M_c) ].                           (6.2)

*Step 3 (rotations; deterministic).* Write M_c = {(a_j, b_j, c_j) : j = 1, ..., n}. For
k >= 1 and distinct indices j_1, ..., j_k, the **rotation** is the following set of q teams:
the teams of M_c with index outside {j_1, ..., j_k}; (a_{j_1}, b, c_{j_1});
(a_{j_{i+1}}, b_{j_i}, c_{j_{i+1}}) for i = 1, ..., k−1; and (a, b_{j_k}, c). Every agent
appears exactly once (sector 1: all a_j and a; sector 3: all c_j and c; sector 2: b in team
j_1, b_{j_i} in team j_{i+1}, b_{j_k} with a and c, the others unmoved), so this is an
assignment, and it avoids e because b is teamed with a_{j_1} ≠ a. Its pairs are: the pairs
(a_j, c_j) and (a, c); the sector-2 pairs of the unmoved teams; and the 2(k+1) new pairs

    (a_{j_1}, b), (b, c_{j_1});  (a_{j_{i+1}}, b_{j_i}), (b_{j_i}, c_{j_{i+1}}), i < k;
    (a, b_{j_k}), (b_{j_k}, c).

Let N be the sum of the deficits of the new pairs. The cost of the rotation is at most
d(a,c) + cost(M_c) + N <= d(a,c) + d(b,c) + cost(M_c) + N, so
W_1 <= d(a,c) + d(b,c) + cost(M_c) + N for every c and every rotation of M_c. Taking c to
be a minimiser in (6.2): rho_e <= N for every rotation of that M_c. Therefore, with

    R(e, c) := {some rotation of M_c has N < delta},

    {rho_e >= delta} ⊆ union over c of R(e, c)^c.                                   (6.3)

(The minimiser depends on the pairs at a and b; the union over all q fixed values of c
removes that dependence.)

*Step 4 (the three independent pieces, c fixed).* On the index set {1, ..., n} of M_c write
j → l if b_j is s-close to the position (a_l, c_l) (j = l allowed), and define

    L := {l : b_l is s-close to (a, c)}    = {l : d(a, b_l) + d(b_l, c) < s},
    J := {j : b is s-close to (a_j, c_j)}  = {j : d(a_j, b) + d(b, c_j) < s}.

Three disjoint families of pairs are involved: (1) the pairs of I_c, which determine M_c
and the relation →; (2) the pairs (a, b_l), (b_l, c), which with M_c determine L; (3) the
pairs (a_j, b), (b, c_j), which with M_c determine J. Hence, conditionally on I_c, the sets
L and J are independent, and each contains each index independently with probability p.
In particular J is independent of (I_c, L) with that law.

Let Gamma^0 := L and Gamma^{i+1} := Gamma^i ∪ {j : j → l for some l ∈ Gamma^i}. So
Gamma^i is the set of indices from which L is reached by a →-path with at most i steps.

*Claim.* If J ∩ Gamma^{r_A} ≠ ∅ then R(e, c) holds. Indeed, let j be in the intersection
and take a →-path j = j_1 → j_2 → ... → j_k ∈ L with the fewest steps; its vertices are
distinct and k − 1 <= r_A. In the rotation along (j_1, ..., j_k), the two new pairs at b
sum to less than s (j_1 ∈ J), the two new pairs at b_{j_i} sum to less than s (j_i → j_{i+1}),
and the two new pairs at b_{j_k} sum to less than s (j_k ∈ L). So
N < (k + 1)·s <= (r_A + 2)·s = delta.

*Step 5 (expansion on I_c).* I_c is an instance with sectors of size n = q − 1 >= q/2 and
independent uniform deficits. By (Q1),
n·p >= q^(2 eps)/4 >= max(13(K+5), 600)·log q, so Lemma U applies to I_c with K + 3 in
place of K:

    P(U(I_c) fails) <= 6·(q−1)^(−(K+3)) + 3·exp(−0.009·(q−1)),

where U(I_c) is the event U(n, s) for the instance I_c. On U(I_c), for any nonempty index
set Rset the positions {(a_l, c_l) : l ∈ Rset} are vertex-disjoint and of type {1, 3}, and
identifying the sector-2 agent b_j of I_c with its index j,
N({(a_l, c_l) : l ∈ Rset}) = {j : j → l for some l ∈ Rset}. So by (G1), whenever L ≠ ∅,

    |Gamma^{i+1}| >= min(lam·|Gamma^i|, 0.15·n),     lam = 0.3·n·p,

and by (G3), |Gamma^{r_A}| >= min(lam^{r_A}·|L|, 0.15·n). (Only U(a) and the cutoff step
(G1) are used here; U(b) is not needed for Theorem A.)

*Step 6 (probability for fixed e and c).* |L| ~ Bin(n, p), so by (C1) with delta = 1/2,
P(|L| < n·p/2) <= exp(−n·p/8) <= exp(−q^(2 eps)/32). Let

    Good := U(I_c) ∩ {|L| >= n·p/2},     G := Gamma^{r_A}.

Good and G are functions of (I_c, L). On Good,

    |G| >= min( 0.3^{r_A}·(n·p)^{r_A+1}/2 , 0.15·n ),   hence
    p·|G| >= min( 0.3^{r_A}·(n·p)^{r_A+2}/(2n) , 0.15·n·p ).

(Correction, 12 November 2026, after the second referee's read: the first version of this
display omitted the factor p in the first argument, and kappa was too large by a factor 4.
Nothing downstream changes.) Here n·p >= q^(2 eps)/4 and
2 eps·(r_A + 1) = 2 eps·ceil(1/(2 eps)) >= 1, so
(n·p)^{r_A+1} >= q^(2 eps (r_A+1))/4^{r_A+1} >= q/4^{r_A+1}. With n <= q this gives

    0.3^{r_A}·(n·p)^{r_A+2}/(2n) = 0.3^{r_A}·(n·p)·[(n·p)^{r_A+1}/n]/2
                                 >= 0.3^{r_A}·(q^(2 eps)/4)·(1/4^{r_A+1})/2 = kappa·q^(2 eps).

Also 0.15·n·p >= 0.0375·q^(2 eps) >= kappa·q^(2 eps), since kappa <= 0.3/128 for r_A >= 1.
Thus p·|G| >= kappa·q^(2 eps) on Good. Since J is independent of (I_c, L),

    P(J ∩ G = ∅ | I_c, L) = (1 − p)^{|G|} <= exp(−p·|G|).

By the Claim, R(e, c)^c ⊆ {J ∩ G = ∅}, so

    P(R(e, c)^c) <= P(Good^c) + E[ 1_Good · P(J ∩ G = ∅ | I_c, L) ]
                 <= P(U(I_c) fails) + P(|L| < n·p/2) + exp(−kappa·q^(2 eps))
                 <= 6·(q−1)^(−(K+3)) + 3·exp(−0.009(q−1))
                    + exp(−q^(2 eps)/32) + exp(−kappa·q^(2 eps)).                  (6.4)

*Step 7 (union bounds).* By (6.1) and (6.3), for each of the q² pairs e of type {1, 2},
P(some optimal assignment uses e and d(e) >= delta) <= sum over the q agents c of
P(R(e, c)^c), each term bounded by (6.4). Summing over the q² pairs of this type and
multiplying by 3 for the three types gives 3q³ times (6.4), which is eta_A. ∎

Remarks. (1) Uniqueness of M is not used. (2) When 1/(2 eps) is an integer, the set
Gamma^{r_A} is only guaranteed to have size of order kappa'·q with a small constant, not
0.15·n; this is enough, and r_A is as in the first version.

---

## 7. Gradients and the variance bound

Write each pair deficit as d_i = exp(−E_i) with E_i independent Exp(1) (so d_i ~ U[0,1] and
the pair value is 1 − exp(−E_i)), i ranging over the m = 3q² pairs. Regard T, LP*, D as
functions of E ∈ [0, ∞)^m.

**Lemma V.** (1) T, LP* and D are continuous, bounded (by 3q) and Lipschitz in each
coordinate (constant 1 for T and LP*, 2 for D).

(2) Almost surely the optimal assignment M is unique, and ∂T/∂E_i = 1{i ∈ M}·d_i for all i.

(3) Optimal solutions x of the LP need not be unique, but almost surely all of them have
the same load vector l* = (l*_i), and ∂LP*/∂E_i = l*_i·d_i for all i.

(4) Consequently, almost surely,

    |grad T|²   = sum over i ∈ M of d_i²        <= 3q,
    |grad LP*|² = sum_i (l*_i)²·d_i²            <= sum_i l*_i = 3q,
    |grad D|²  <= 2·|grad T|² + 2·|grad LP*|².

*Proof.* For a point l in R^m write phi_l(E) := sum_i l_i·(1 − exp(−E_i)); it is smooth with
∂_i phi_l = l_i·exp(−E_i) = l_i·d_i. T is the maximum of phi_l over the finitely many
indicator vectors l of assignments (viewed as pair sets). LP* is the maximum of
sum_e v(e)x_e = phi_{l(x)}(E) over x ∈ P, that is, the maximum of phi_l over the **load
polytope** Lam := {l(x) : x ∈ P}, a polytope in [0,1]^m (the linear image of P), and hence
the maximum of phi_l over the finitely many vertices of Lam.

(1) For 0 <= l_i <= 1 each phi_l is 1-Lipschitz in each coordinate with values in [0, 3q]
(as sum_i l_i = 3q); a maximum of finitely many such functions is continuous with the same
properties, and D is a difference of two of them with 0 <= D <= 3q.

(2), (3) Let V be a finite set of points of R^m (the assignment vectors, or the vertices of
Lam). In the coordinates X_i = 1 − d_i, phi_l = <X, l>. If the maximiser over V is not
unique at X, then <X, l − l'> = 0 for two distinct l, l' ∈ V: a finite union of hyperplanes,
which is Lebesgue-null. E ↦ X is a diffeomorphism of (0, ∞)^m onto (0, 1)^m, so the
exceptional set is null in the E coordinates too. Where the maximiser l* over V is unique,
it stays the unique maximiser on a neighbourhood (finitely many strict inequalities), so
the maximum equals phi_{l*} there and its partial derivatives are l*_i·d_i. For T, l* is the
indicator of M. For LP*, uniqueness of the maximising vertex of Lam means that the set of
optimal loads is the single point l* (the optimal face of Lam is the convex hull of the
optimal vertices), i.e. every optimal x has l(x) = l*. Non-uniqueness of x itself is
possible because x ↦ l(x) maps R^(q³) to R^(3q²) and is not injective; it is harmless,
since the objective and the gradient depend on x only through l(x).

(4) follows from (2), (3), d_i <= 1, 0 <= l*_i <= 1, sum_i l*_i = 3q, and
(u − w)² <= 2u² + 2w². ∎

**Proposition W (variance bound with the bad events written out).** Let delta_A, delta_B > 0
and let

    G_A := {every pair used by M has deficit < delta_A},
    G_B := {every team in the support of every LP optimum has deficit < delta_B},

with P(G_A fails) <= eta_A and P(G_B fails) <= eta_B. Then

    Var T   <= 4·[ delta_A·E(3q − T) + 3q·eta_A ],
    Var LP* <= 4·[ delta_B·E(3q − T) + 3q·eta_B ],
    Var D   <= 8·[ (delta_A + delta_B)·E(3q − T) + 3q·(eta_A + eta_B) ].

*Proof.* By Lemma V(1), Lemma P applies to T, LP* and D, and by Lemma V(2)–(4)

    Var T <= 4·E sum_{i∈M} d_i²,   Var LP* <= 4·E sum_i (l*_i)² d_i²,
    Var D <= 8·E[ sum_{i∈M} d_i² + sum_i (l*_i)² d_i² ].

On G_A: sum_{i∈M} d_i² <= delta_A·sum_{i∈M} d_i = delta_A·(3q − T). Off G_A the sum is at
most 3q (Lemma V(4)). So E sum_{i∈M} d_i² <= delta_A·E(3q − T) + 3q·eta_A.

On G_B: if l*_i > 0, take any optimal x; l_i(x) = l*_i > 0, so some support team contains
the pair i, and d_i <= (deficit of that team) < delta_B. With l*_i <= 1,

    sum_i (l*_i)² d_i² <= delta_B·sum_i l*_i·d_i = delta_B·(3q − LP*) <= delta_B·(3q − T),

using LP* >= T. Off G_B the sum is at most 3q. So
E sum_i (l*_i)² d_i² <= delta_B·E(3q − T) + 3q·eta_B. ∎

---

## 8. Proof of Theorem 1

Let eps ∈ (0, 1/2), s = q^(−1/2+eps), and q >= q_0(eps, K).

(i) Theorem A (its hypothesis is (Q1)) gives: off an event of probability eta_A, no optimal
assignment uses a pair of deficit >= delta_A := (r_A + 2)·s = (ceil(1/(2 eps)) + 1)·s
<= C_eps·s. Theorem B ((Q1), (Q2)) gives: off an event of probability
eta_B := 6·q^(−K) + 3·exp(−0.009 q), every support team of every LP optimum has deficit
< delta_B := 3(r_B + 1)·s = 3(ceil(1/(2 eps)) + 2)·s = C_eps·s. The total failure
probability is at most

    eta_A + eta_B <= (144·2^K + 6)·q^(−K)
                     + 3q³·[3·exp(−0.009(q−1)) + exp(−q^(2 eps)/32) + exp(−kappa·q^(2 eps))]
                     + 3·exp(−0.009 q).

By (Q3) the exponentially small terms are together below q^(−K), so
eta_A + eta_B <= (144·2^K + 7)·q^(−K) <= 150·2^K·q^(−K).

(ii) Apply Proposition W with these delta_A, delta_B and K = 2, using
E(3q − T) <= C_0·q^(1/3) (Import F):

    delta·E(3q − T) <= C_eps·q^(−1/2+eps)·C_0·q^(1/3) = C_eps·C_0·q^(−1/6+eps),
    3q·eta_A <= 3q·(576 + 1)·q^(−2) <= 1731/q,     3q·eta_B <= 3q·7·q^(−2) = 21/q,

(the exponentially small parts of eta_A and of eta_B each being below q^(−2) by (Q3)).
Hence Var T <= 4·C_eps·C_0·q^(−1/6+eps) + 6924/q, Var LP* <= 4·C_eps·C_0·q^(−1/6+eps) + 84/q
and Var D <= 16·C_eps·C_0·q^(−1/6+eps) + 14016/q, which are within the stated bounds.

(iii) D >= 0, so {D = 0} ⊆ {|D − E D| >= E D} and Chebyshev gives
P(D = 0) <= Var D/(E D)². If E D >= q^(−1/12+eps), then (E D)² >= q^(−1/6+2 eps), and (ii)
with the same eps gives P(D = 0) <= 16·C_eps·C_0·q^(−eps) + 28000·q^(−5/6−2 eps) → 0. If
liminf E D > 0, then E D >= q^(−1/12+1/24) for large q, and the previous sentence applies
with eps = 1/24. The core is nonempty exactly when D = 0 (LP duality: a core payoff is a
dual solution of value T). ∎

---

## 9. What is imported

**Import F (the only non-textbook input).** Frieze, Pegden and Tkocz, "On random
multi-dimensional assignment problems", arXiv:1901.07167v1 (22 January 2019; journal
reference Discrete Applied Mathematics 287 (2020), 1–9). Opened in this session (arXiv v1
PDF, pages 1–2). The model, their display (1): "W_{i,j,k} = a_{i,j} + b_{i,k} + c_{j,k},
1 <= i, j, k <= n, where the a_{i,j}, b_{i,k}, c_{j,k} are independent uniform [0, 1] random
variables"; the objective is to minimise sum_i W_{i,sigma(i),tau(i)} over permutations
sigma, tau; Z_n is the optimum. Verbatim:

> **Theorem 2.** There exist constants c_1, c_2 such that (a) E(Z_n) >= c_1 n^{1/3} and (b)
> GREEDY(n − n^{1/4}) finds a solution of expected value at most c_2 n^{1/3}. In this case
> Step 6 can be completed by choosing an arbitrary completion.

(Step 6 of their Algorithm 1 completes the partial assignment of m = n − n^{1/4} triples to
a full one.) Their Z_n is our 3n − T_n. How it is used: the completed greedy output is a
perfect assignment, so Z_n is at most its cost. Whether or not the bound in (b) is read as
including the completion, the completion adds at most n^{1/4} + 1 triples of cost at most 3
each, so E(3n − T_n) <= c_2·n^(1/3) + 3·n^(1/4) + 3 <= (c_2 + 6)·n^(1/3). If the statement is
read as holding only for large n, then, since 3n − T_n <= 3n always,
C_0 := sup_n E(3n − T_n)/n^(1/3) is still finite. Only finiteness of C_0 is used, and only
in Theorem 1(ii) and (iii). Part (i), Theorems A and B and Lemmas U, G, P, V do not use it.

Robustness: the arXiv listing page was reported to this session as describing the pair
costs as exponential rate 1, unlike the v1 PDF. If a later version is stated for Exp(1)
pair costs, Import F still follows for uniform deficits by the coupling
d = 1 − exp(−E) <= E, under which every assignment is cheaper.

**Textbook facts, quoted without proof:** the Chernoff bounds (C1)–(C3); König's theorem in
the form (M); LP strong duality and complementary slackness (LP); the Efron–Stein inequality
(ES); absolute continuity and almost-everywhere differentiability of Lipschitz functions on
the line; Fubini's theorem; Chebyshev's inequality.

**Not imported:** the Poincaré inequality (proved as Lemma P); Rademacher's theorem (not
needed, since only coordinate-wise derivatives are used); Harris's inequality and the
Fulkerson–Hoffman–Oppenheim theorem (not used in Theorem 1).

---

## 10. What changed relative to `reports/pairdeficit.md`

No mathematical step was found to be wrong. The changes are:

1. **Reviewers' gap (a), rotation estimate.** §6 Steps 4–7: conditioning on (I_c, L), the
   event Good, the Bernoulli miss bound on Good, the two exceptional probabilities added,
   and the union over the q values of c and the 3q² pairs with the final bound eta_A.
   One consequence that the first version left implicit: Lemma U is applied to 3q³
   sub-instances I_c, so it is invoked with K + 3 in place of K.
2. **Reviewers' gap (b), cutoff.** New Lemma G (§4): sub-selection of floor(1/p) positions
   giving 0.3·n·(1 − p) >= 0.15·n, the passage to U(b) via n·p >= 20, and the induction.
   Step counts verified: r_B = ceil(1/(2 eps)) + 1 and r_A = ceil(1/(2 eps)) − 1 both stand
   as first stated. Theorem B needs the extra largeness condition (Q2). Theorem A does not
   use U(b) at all.
3. **(c)** Lemma U(b): rate KL(0.9 ‖ 1 − e^(−3)) = 0.020879, written 0.0208; the union-bound
   exponent is then 0.0091 >= 0.009.
4. **(d)** Variance: Proposition W with the bad-event terms 3q·eta; Lemma P stated with its
   regularity hypotheses and proved; Lemma V on almost-sure uniqueness of loads versus
   non-uniqueness of x.
5. **Lemma U's probability** is stated as 6·n^(−K) + 3·exp(−0.009·n) instead of O(n^(−K)):
   under the stated hypotheses the part (b) term is exponentially small in n but is not
   bounded by n^(−K) with a constant independent of K.
6. **Step 1 of Theorem A** is stated with rho_e >= delta and a strict rotation bound
   N < delta, so almost-sure uniqueness of M is not needed there.
7. The dual variable is called u (the first version called both the dual and the closeness
   probability p).

---

## 11. Remaining gaps and caveats

- **Import F** was read in the arXiv v1 PDF only. The published version (Discrete Applied
  Mathematics 287) was not opened and should be compared before circulation; the coupling
  remark in §9 covers the one discrepancy reported (uniform versus exponential pair costs).
  The proof of FPT's Theorem 2 was not checked here.
- The textbook facts in §9 are quoted from memory, not from an opened source.
- No step of §§2–8 is knowingly incomplete. The document has not been read by a second
  agent or by a person. The first version's numerical checks (`v3/`) concern the same
  Lemma U constants and the deterministic steps of Theorems A and B, and remain applicable;
  Lemma G, the union-bound arithmetic of §6 and §8, and Lemma V are new text and unchecked
  by anyone else.
- The thresholds q_0(eps, K) are astronomically large; nothing here is informative at
  computable sizes.
