# Dual route to vanishing of the core: status (2026-11-11)

Written by a Claude Code lane agent in a 60-minute time box. Not reviewed; not reproduced by a
different model. Floating-point numerics are evidence, not certificates. Scripts and outputs in
this directory (`lex.py`, `width.py`, `spread.py`, `out/`). Notation as in `../../THEOREM_VANISHING.md`:
M = integer optimum (a.s. unique), relabelled so that M = {(x_i, y_i, z_i)}; A_ij = A[x_i, y_j],
B_ik = B[x_i, z_k], C_jk = C[y_j, z_k]; v_ijk = A_ij + B_ik + C_jk; v_i = v_iii; S = sum of the three
bipartite optima (of A, B, C separately).

**Outcome.** No proof of Conjecture V. Four complete (unreviewed) results, of which the main one is
negative for idea 4 of the brief: in the ordinal (infinitely spread) limit the core is *always*
nonempty (Theorem O). Two dual routes are closed or collapse onto Lemma C (D2, D3). Numerics on the
core shape (W1) contradict the "core width shrinks" form of idea 2.

---

## D1. The core in potential coordinates. Status: complete.

Write a_i = p(x_i), b_i = p(y_i), c_i = p(z_i). Then
core(G) = { (a,b,c) >= 0 : a_i + b_i + c_i = v_i for all i, a_i + b_j + c_k >= v_ijk for all i,j,k }.
*Proof.* p(N) = T = sum_i v_i and p(e) >= v(e) on the q M-teams force equality on M. QED.

(i) **Two-equal constraints are three bipartite potential systems.** Using the equalities, the
constraint at (i,i,k) reads c_k - c_i >= W^Z_ik := (B_ik + C_ik) - (B_ii + C_ii); likewise
at (k,i,i): a_k - a_i >= (A_ki + B_ki) - (A_ii + B_ii), and at (i,k,i):
b_k - b_i >= (A_ik + C_ki) - (A_ii + C_ii). Each system alone is feasible: a positive
cycle in W^Z would be a cyclic reassignment of the Z-agents among the fixed (x_i, y_i) pairs that
beats M. So the two-equal constraints never cause emptiness; emptiness needs the q^3 - 3q^2 + 2q
constraints with i, j, k distinct, which couple a, b, c.

(ii) **Dimension.** If M is the unique LP optimum, the core (when nonempty) has dimension exactly 2q.
*Proof.* Goldman-Tucker strict complementarity for the pair (max v.lam, lam >= 0, loads <= 1) /
(min p(N), p >= 0, p(e) >= v(e)): there is a dual optimum p whose constraint e is slack for every
team e outside the support of every primal optimum (all e not in M) and with p_i > 0 for every
player whose load is tight in every primal optimum (all players). So p is a relative-interior point
of the affine space {p : p(M_i) = v_i}, which has dimension 3q - q = 2q. QED.
(A.s. uniqueness of M as LP optimum when the core is nonempty: the LP optimum is a.s. unique for
continuous iid values since two distinct vertices tie on a null set.)

## D2. Pair-price (separable) certificates. Status: complete. Closes idea 1 in its natural form.

Call a core point *separable* if p = p^A + p^B + p^C where p^A is a feasible dual of the XY
assignment problem (p^A(x) + p^A(y) >= A[x,y], supported on X and Y), and likewise p^B (XZ, for B)
and p^C (YZ, for C). This is the "split each M-team's value across its three pairs and price each
pair bipartitely" decomposition.

**Proposition D2.** (a) A separable core point exists iff T = S. (b) For iid continuous entries,
P(T = S) = 1/q!.
*Proof.* (a) If p is separable then p(N) >= (dual objective of p^A) + ... >= S_A + S_B + S_C = S by
bipartite weak duality, while p(N) = T <= LP* <= S; so T = S. Conversely, if T = S take optimal
bipartite duals; their sum p satisfies p(e) >= A + B + C = v(e) on every team and p(N) = S = T.
(b) T = S iff one assignment projects onto optimal assignments of all three bipartite problems,
i.e. sigma_C o sigma_A = sigma_B for the a.s. unique optimal bijections sigma_A: X->Y,
sigma_B: X->Z, sigma_C: Y->Z. Let pi be a uniform permutation of Y independent of (A,B,C), and
A'[x,y] = A[x, pi(y)]. Then (A', B, C) has the law of (A, B, C), and sigma_{A'} = pi^{-1} o sigma_A.
Hence P(T = S) = P(sigma_C pi^{-1} sigma_A = sigma_B) = P(pi^{-1} = sigma_C^{-1} sigma_B sigma_A^{-1})
= 1/q!, since pi is uniform and independent of the
right-hand side. QED.

Consequence: every core point of a large random market is non-separable; the coupling through the
distinct-index constraints of D1(i) is essential. "Three bipartite dual conditions plus one
consistency constraint" is exactly D2 and the consistency fails with probability 1 - 1/q!, but this
only proves that *separable* certificates fail, not that the core is empty (gate 1: the core is
nonempty with substantial probability at q = 10 while P(T = S) = 1/10! ). Any dual obstruction has to
handle non-separable p.

## D3. Dual depth = -mu*. Status: complete. The "random polytope" form of idea 2 is Lemma C.

Define the dual depth delta* = max { t : a_i + b_i + c_i = v_i, a_i + b_j + c_k >= v_ijk + t for
(i,j,k) not all equal } (no sign constraint on t).
Here p is free (the dual of the padded LP with unit loads, PROOF_NOTES Thm 1 Step 1, which has the
same optimum LP* as the formulation with loads <= 1 and p >= 0).
**Claim.** delta* = -mu* (mu* of Lemma 2); in particular core(G) nonempty iff delta* >= 0.
(`width.py` adds p >= 0 to the depth LP; this can only lower the reported depth.)
*Proof.* LP duality for the objective v + mu 1_M over the padded polytope P: M is optimal iff there is
p with p(e) >= v(e) + mu [e in M] for all e and p(N) = T + mu q. Put p' = p - mu/3 on every player:
p'(e) >= v(e) - mu off M, p'(e) >= v(e) on M, p'(N) = T, hence equality on M. So M is optimal for
v + mu 1_M iff delta* >= -mu, and by Lemma 2 iff mu >= mu*. QED.
So "show the core polytope is empty by an in-radius / depth argument" is literally Lemma C; the dual
view does not give an independent scalar. Numerics (W1) below agree: depth ~ the F3 values of -mu*.

## Theorem O (ordinal limit: the core is always nonempty). Status: complete, unreviewed.

**Setting.** Fix q and a strict order e_1 > e_2 > ... > e_{3q^2} on the 3q^2 pairs. For lam in the
padded axial polytope P let w_e(lam) = sum of lam over teams containing pair e. Let lam^lex be a
lexicographic maximiser of (w_{e_1}, w_{e_2}, ...) over P.

**Theorem O.** w_e(lam^lex) is 0 or 1 for every pair e; hence lam^lex is integral (0/1 XY- and
XZ-marginals with unit loads determine an integral lam).
**Corollary O.** For each q there is K_q < infinity such that if pair values satisfy
val(e_k) >= K_q val(e_{k+1}) > 0 for all k, the core is nonempty. If pair values are exp(K E) with E
iid continuous, then for each fixed q, P(core nonempty) -> 1 as K -> infinity.

*Proof of Theorem O.* Run the lexicographic optimisation pair by pair. Maintain the state after k
steps: a partition of the 3q agents into *components* (connected parts of the graph of pairs fixed
at weight 1). Invariant (I_k): every fixed value so far is 0 or 1; each component is a singleton, a
two-agent "pair-component" of two different sectors, or a full team; and the feasible set
P_k = {lam in P : earlier weights fixed} equals the set F_k of lam in P supported on teams that are a
union of components (call such lam *compatible*).

Let q' = number of agents of one sector lying in components that are not full teams,
n_X, n_Y, n_Z = numbers of singleton components per sector, p_XY, p_XZ, p_YZ = numbers of
pair-components. Counting agents of each sector in unfinished components,
q' = n_X + p_XY + p_XZ = n_Y + p_XY + p_YZ = n_Z + p_XZ + p_YZ, so
D := n_Z - p_XY = n_Y - p_XZ = n_X - p_YZ = q' - (p_XY + p_XZ + p_YZ) is a single integer.

*Lemma O1.* F_k is nonempty iff D >= 0, and then it contains an integral point.
In a compatible lam, a pair-component of type XY lies only in teams (XY-component + Z-singleton)
(any other Z-agent is in a component containing an X- or Y-agent already). Summing loads, the
Z-singletons give mass n_Z, of which p_XY is used by XY-components; so n_Z >= p_XY, i.e. D >= 0. If
D >= 0, match each XY-component with its own Z-singleton (and XZ with Y, YZ with X), and group the
remaining D singletons of each sector into D triples: an integral compatible assignment. QED.

*Lemma O2.* Under (I_k), let e = e_{k+1} join agents s in component C1 and t in C2. Then
max {w_e(lam) : lam in F_k} is 0 or 1, as follows:
 (a) C1 = C2: value 1 (already fixed).
 (b) C1, C2 contain agents of a common sector: value 0 (a team containing e would contain both).
 (c) C1 cup C2 is a full team (pair-component + complementary singleton): value 1 (merge; D unchanged,
     Lemma O1 completes it integrally).
 (d) C1, C2 singletons of sectors, say, X and Y: value 1 if D >= 1 (merge into an XY-component, D
     drops by 1, Lemma O1 completes); value 0 if D = 0, since then the n_Z = p_XY Z-singleton mass is
     used up by XY-components, so no compatible team contains three singletons, and two singletons
     x, y meet only in such teams.
These cases are exhaustive: two components of size >= 2 always share a sector, and a component of
size 2 plus a singleton is either case (b) or (c).
*Maintaining (I_{k+1}).* If the value is 1, merge C1 and C2; P_{k+1} = {lam in P_k : w_e = 1} =
compatible lams for the new partition (a compatible lam with w_e = 1 puts all mass of s on teams
containing t, hence on teams containing C1 cup C2). If the value is 0, every lam in
P_k = F_k already has w_e = 0, so P_{k+1} = P_k = F_k = F_{k+1} (partition unchanged). Either way
(I_{k+1}) holds. At k = 0 all components are singletons and P_0 = P = F_0. QED (Theorem O).

*Proof of Corollary O.* P has finitely many vertices; lam^lex is the unique lex-maximal vertex in
marginal space (Theorem O makes all its pair weights 0/1, so it is the integral assignment they
define). For any other vertex w, let k be the first index with w_{e_k} != lam^lex_{e_k}; then
lam^lex_{e_k} - w_{e_k} >= g_q > 0, where g_q = min |w_e - w'_e| over pairs e and vertices w, w' of P with w_e != w'_e. Since every
pair weight lies in [0,1], the objective difference is >= val(e_k) (g_q - sum_{j>=1} K_q^{-j})
= val(e_k) (g_q - 1/(K_q - 1)) > 0 once K_q > 1 + 1/g_q. So the LP optimum is the integral lam^lex: LP* = T. For exp(K E) values
with E continuous, the ratio of consecutive ordered values exceeds K_q with probability -> 1 as
K -> infinity (fixed q, finitely many gaps, each positive a.s.). QED.

Numerical check (`lex.py`, sequential HiGHS LPs on random orders): 0 fractional sequential values in
2,000 instances at q = 3 and q = 4, 100 at q = 2, 5, 6, 8, and 60 at q = 10, 12.

**What Theorem O says about the target.**
- Idea 4 ("a heavier-tailed law may be easier") points the wrong way: at fixed q, spreading the law
  drives P(empty) to 0, not 1 (see W2 for the finite-K sweep). The gate-1 Exp contrast (share falling
  more slowly than U[0,1]) is the first step of this.
- Barrier (informal): Conjecture V is false for the ordinal model, so no argument that uses only the
  relative order of the 3q^2 pair values can prove it. Emptiness is created by *cardinal* near-ties:
  sums of different pairs competing at the scale q^{-2/3} (THEOREM_VANISHING F1/F2). Any proof must
  quantify that competition, which is what Lemma C does through mu*.
- The limits do not commute: lim_q P_U(empty) (conjecturally 1) versus lim_{K} P_K(empty) = 0 at every
  fixed q. A law-free monotonicity argument in q (idea 3) would have to fail for the ordinal model,
  where P(empty) = 0 for all q; so any monotonicity must use the law.

---

## Numerics

### W1. Core shape on nonempty U[0,1] instances (`width.py`, `out/width_q*.jsonl`)

Width = range over the core of p_i - (mean of p over i's sector), which removes the sector-shift
directions (Theorem 2 shift lemma); wmax = max over the 3q agents, wmed = median. depth = delta* of
D3. gapMmean = mean over M-teams of 3 - v_i. n(eps) = number of non-M teams with slack <= eps at the
depth-LP optimum. Gate-1 seeds; nonempty instances only.

| q | nonempty / n | wmax | wmed | depth delta* | gapMmean | n(0.01) | n(0.03) | n(0.1) | 2q |
|---|---|---|---|---|---|---|---|---|---|
| 6 | 102/120 | 0.665 | 0.400 | 0.0376 | 0.725 | 2.5 | 6.3 | 14.1 | 12 |
| 10 | 85/120 | 0.617 | 0.322 | 0.0180 | 0.539 | 8.2 | 16.9 | 33.2 | 20 |
| 16 | 38/120 | 0.493 | 0.225 | 0.0074 | 0.397 | 23.0 | 38.1 | 69.6 | 32 |
| 24 | 20/160 | 0.431 | 0.179 | 0.0033 | 0.305 | 47.2 | 63.5 | 128.2 | 48 |

(Means over nonempty instances; conditioning on nonempty biases everything; 20 instances at q = 24.)

Reading: the core is not uniformly thin. The median agent range shrinks slowly (about q^{-0.6}
from q = 6 to 24) while the depth (= -mu*, D3) falls about q^{-1.7}, and the number of non-M teams
within slack 0.01 of the deepest point grows from 0.2 to about 1.0 times the dimension 2q. The
nonempty core is a flat 2q-dimensional polytope, thin in the few directions of the mesoscopic
odd-cycle structures of F2 and long in the others. So the "width shrinks like q^-alpha faster than
the fluctuations" form of idea 2 has to be stated with the depth, not the agent ranges, and then it
is Lemma C (D3). The near-tight count reaching about 2q is suggestive of a random-LP feasibility
threshold (constraints near tight ~ dimension), but no proof was attempted.

### W2. Spread sweep (`spread.py`, `out/spread.json`)

Pair values exp(K (U - 1)); share of instances with fractional LP optimum (empty core).

| q \ K | 0 (U[0,1]) | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|---|
| 4 | 0.040 | 0.018 | 0.010 | 0.003 | 0.003 | 0.000 |
| 8 | 0.225 | 0.195 | 0.165 | 0.110 | 0.030 | 0.007 |
| 12 | 0.405 | 0.385 | 0.352 | 0.253 | 0.130 | 0.025 |

(400 instances per cell, gate-1 seeds reused across K, so cells within a row are coupled.)
P(empty) decreases monotonically in K at every q, as Corollary O predicts in the limit; at each fixed
K it still increases with q (up to 12), so the conjecture may well hold for each fixed K, with the
onset pushed to larger q as the spread grows.

---

## Most promising route, and what remains

1. Through the dual nothing independent of Lemma C was found: D3 identifies the dual depth with -mu*,
   D2 kills separable certificates, W1 kills uniform-width arguments. The dual's genuine content is
   D1(i): emptiness lives entirely in the distinct-index constraints.
2. Theorem O is the main new structural fact: emptiness is a cardinal near-tie phenomenon. Its proof
   suggests a route for U[0,1] that is different from Lemma C: run the same component/D bookkeeping
   in the regime where pair values are within O(q^{-2/3}) of 1 and show that the lexicographic
   (greedy) bookkeeping breaks down at Theta(q) places where two *sums* of near-1 pairs compete
   within o(q^{-2/3}), each giving an odd-cycle (half-integral) improving direction. This is a
   heuristic plan only; nothing of it is proved.
3. Gaps: Conjecture V open; Lemma C open; no uniform lower bound on P(empty); Theorem O and D1-D3
   unreviewed; all tables are floating point.

Confidence: high in D1, D2, D3 and Theorem O (short, checked against simulation); low that the dual
view alone yields Conjecture V.
