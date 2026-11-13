# Round 1 report: pairdeficit line

Saved by the orchestrator from the agent's message (the agent's own write of a report
file was refused by the harness). Scripts and raw outputs: `../pairdeficit/` (pd.py,
pd_team.py, test_chain.py, test_rot.py, pd_q*.jsonl, logs). Nothing here is verified.

## Result

The lemma is proved with exponent 1/2 instead of 2/3, for M (pairs) and for the LP support
(whole teams). 1/2 > 1/3, which is all that is needed: Var T, Var LP*, Var D =
O(q^(-1/6+eps)) -> 0, so P(core nonempty) -> 0 follows from E D >= q^(-1/12+eps), e.g. from
liminf E D > 0. Exponent 2/3 is open for the maximum; for the LP a version "all but
q^(1-3eps) of the mass" at scale q^(-2/3+eps) is sketched (Prop C).

## Notation

d = deficit. A "position" is a pair of agents from two different sectors. A third-sector
agent z is s-close to a position if its two deficits to the position's agents sum to less
than s. p := s^2/2 = P(U1 + U2 < s).

## 1. Sufficient exponent (complete given facts 1 and 3 of the brief)

In E-coordinates dT/dE_i = 1{i in M}·d_i and dLP*/dE_i = load_i·d_i. So
Var T <= 4 E sum_M d_i^2, Var LP* <= 4 E sum load_i^2 d_i^2, and
Var D <= 8 E[sum_M d_i^2 + sum load_i^2 d_i^2].

If used and loaded deficits are at most delta off an event of probability eta, then
Var D <= 16·delta·E(3q - T) + 48·q·eta. With delta = q^(-alpha): Var D = O(q^(1/3 - alpha)),
and P(D = 0) -> 0 once E D >> q^(1/6 - alpha/2).

- alpha = 0: q^(1/6) (Proposition R).
- alpha = 1/3: E D -> infinity at any rate.
- alpha = 1/2 - eps (proved below): E D >> q^(-1/12 + eps/2).
- alpha = 2/3 - eps: E D >> q^(-1/6 + eps/2).

So the threshold for Var T = o(1) is any exponent above 1/3.

First-moment check: forcing one pair of deficit >= delta into an assignment of total
deficit t multiplies the union bound of fact 4 by 3q(1 - delta/t)^(3q-1), about
3q·exp(-3·delta·q^(2/3)/C). Against exp(Theta(q)) entropy this needs delta ~ q^(1/3) > 1.
First moments cannot work.

## 2. Lemma U (uniform expansion; claimed complete)

Sector size n, p = s^2/2, with np >= 13(K+2) log n and n >= (600/p) log n. With probability
1 - O(n^(-K)), for all three position types and EVERY set S of vertex-disjoint positions,
writing N(S) for the third-sector agents s-close to some position of S:

(a) |N(S)| >= 0.3·np·|S| if 1 <= |S| <= 1/p;
(b) |N(S)| >= 0.9·n if |S| >= 3/p.

Proof. For fixed S of size m the events use distinct pairs, so
|N(S)| ~ Bin(n, 1 - (1-p)^m). (a) The success probability is at least 0.632·pm; Chernoff
gives exp(-0.079·npm); there are at most n^(2m) sets. (b) For |S| = ceil(3/p) the success
probability is at least 0.95; P(Bin < 0.9n) <= exp(-0.0209·n); there are at most
n^(2·ceil(3/p)) sets. QED.

Because it holds for all S at once, S may depend on M or on the LP optimum; this replaces
conditioning. It needs np of order at least log n, i.e. s of order at least
(log q / q)^(1/2): the barrier of the method.

## 3. Theorem B (LP; claimed complete)

Fix eps, s = q^(-1/2+eps), r = ceil(1/(2 eps)) + 1. With probability 1 - O(q^(-K)) for
every K, every support team f of every LP optimum has d(f) < 3(r+1)s.

Proof, on the event of Lemma U (n = q). Let x be optimal, p an optimal dual (p(g) >= v(g),
with equality on the support), and f = (a,b,c) in the support.

Start: p(a,b',c) >= v(a,b',c) and p(f) = v(f) give
p_b' - p_b >= d(a,b) + d(b,c) - [d(a,b') + d(b',c)], so
p_b' > p_b + d(a,b) + d(b,c) - s for b' in B_0 := N({(a,c)}).

Chain: if g = (a_g, b', c_g) is a support team and b'' is s-close to (a_g, c_g), then
p_b'' >= p_b' + d(a_g,b') + d(b',c_g) - [d(a_g,b'') + d(b'',c_g)] > p_b' - s.

Growth: the support teams whose sector-2 agent is in B_i carry mass |B_i|; their
(a,c)-positions carry it as a fractional bipartite matching, so they contain a matching S_i
with |S_i| >= |B_i|. Put B_(i+1) = N(S_i') for a sub-matching of size
min(|B_i|, floor(1/p)), or ceil(3/p) once |B_i| >= 3/p. By Lemma U,
|B_i| >= min((0.3np)^(i+1), 0.15n), then >= 0.9n one step later. So |B_r| >= 0.9q, and
p_b'' > p_b + d(a,b) + d(b,c) - (r+1)s on B_r.

Symmetry: there are sets A_r, C_r of size >= 0.9q with
p_a'' > p_a + d(a,b) + d(a,c) - (r+1)s and p_c'' > p_c + d(a,c) + d(b,c) - (r+1)s.

Close: the support teams with an agent outside A_r, B_r, C_r carry mass at most 0.3q < q,
so some support team t lies in A_r x B_r x C_r. It is tight:
3 >= v(t) = p(t) > p(f) + 2d(f) - 3(r+1)s = 3 + d(f) - 3(r+1)s. QED.

This uses complementary slackness, so it says nothing about M when D > 0.

## 4. Theorem A (M; claimed complete)

Fix eps, s = q^(-1/2+eps), r = ceil(1/(2 eps)) - 1. With probability 1 - O(q^(-K)), M uses
no pair of deficit >= (r+2)s.

Proof for e = (a,b).

Threshold: let W_0 be the least cost of an assignment using e with d(e) counted as 0, and
W_1 the least cost of an assignment avoiding e. Both are independent of d(e), and M uses e
iff d(e) < rho_e := W_1 - W_0. So P(M uses e, d(e) >= delta) <= P(rho_e > delta).

Decoupling (this is what avoids conditioning on M): let I_c be the instance minus a, b, c,
and M_c = {(a_j, b_j, c_j)} its optimal assignment. Then
W_0 = min over c of [d(a,c) + d(b,c) + cost(M_c)], and M_c is independent of every pair
touching a, b or c.

Rotation: for distinct j_1..j_k move b into team j_1, b_(j_i) into team j_(i+1), and
b_(j_k) into (a,c). The result avoids e and costs at most
d(a,c) + d(b,c) + cost(M_c) + N, where N is the sum of the 2(k+1) new pair deficits. At the
minimising c: rho_e <= N. So
P(rho_e > delta) <= sum over c of P(no rotation of M_c has N <= delta).

Fixed c: write j -> l if b_j is s-close to (a_l, c_l); L = {l : b_l is s-close to (a,c)};
J = {j : b is s-close to (a_j, c_j)}. Given I_c, L and J are independent Bernoulli(p)
subsets of M_c's teams (L uses pairs touching a or c, J uses pairs touching b). Let
Gamma^r(L) be the teams with a directed path of length <= r into L. Lemma U(a) on I_c
gives |Gamma^i(L)| >= min((0.3np)^i·|L|, 0.15n), and |L| >= np/2 off probability
exp(-np/8). J is independent of (I_c, L), so
P(J misses Gamma^r(L)) <= exp(-p·min((0.3np)^r·np/2, 0.15n)) = exp(-Omega(q^(2 eps))).
If j is in J and in Gamma^r(L), the rotation along a shortest path has every move < s, so
N < (r+2)s. Union over the q values of c and the 3q^2 pairs. QED.

With k moves the same proof gives the threshold (C log q / q^(k-1))^(1/(2k)): k = 2
(2-swaps) gives (C log q / q)^(1/4), k = 3 gives (C log q / q^2)^(1/6).

## 5. Corollary

Var T, Var LP*, Var D = O(q^(-1/6+eps)); P(D = 0) <= Var D/(E D)^2 -> 0 if
E D >= q^(-1/12+eps). Rests on fact 1 (re-derived by the agent) and fact 3 (FPT; the agent
did not open the paper).

## 6. Proposition C (LP at the true scale; SKETCH, a gap until written out)

Let u := 1 - p, so u(g) <= d(g) with equality on the support.

(C1) Janson on triangles of the random tripartite graph {d < kappa/3}: whp any three sets
of at least beta·q agents, one per sector, contain a team with d < kappa =
(3.9/beta)·q^(-2/3). Top-(1-alpha) sets give Q_alpha := sum over sectors of the
alpha-quantile of u < kappa; for alpha > 2/3 the bottom-alpha sets contain a tight team, so
Q_alpha >= 0. Hence quantiles 0.67 to 0.99 of u in each sector lie in a window of width
about 390·q^(-2/3): 32% of each sector has flat prices.

(C2) Uniformly, all but O(q/K^3) agents have a team of deficit < K·q^(-2/3) into the
windows, for K^3 >= 530 log q, which bounds their u from above; tight teams then bound u
from below for all but O(q/K^3) more.

Conclusion: all but O(q/K^3) of the LP mass sits on teams of deficit at most
(390 + 3K)·q^(-2/3). With Theorem B for the rest, Var LP* = O(q^(-1/4+eps)).

## 7. Gaps

(i) Maximum used deficit at q^(-2/3+eps): open. One-sector moves cost two new pairs per
factor q, so they stop at q^(-1/2); single exceptional agents are the obstruction.
(ii) Prop C is not written to referee standard.
(iii) M has no analogue of Prop C beyond Markov.
(iv) Optimised logs (s ~ (C log q / q)^(1/2), r ~ log q / log log q) not written.

## 8. Numerics (floating point, evidence only; seeds tag 31)

n = 100, 100, 100, 50, 40 at q = 16, 24, 32, 40, 48.

- Largest used pair deficit × q^(2/3), mean: M 3.17, 3.37, 3.55, 3.67, 3.65; LP support
  3.35, 3.67, 3.99, 3.97, 4.14.
- E sum_M d^2 × q^(1/3): 3.89, 3.87, 3.94, 3.97, 3.86 (so the Poincaré bound is
  Var T <= about 15.6·q^(-1/3); measured Var T is about 0.23).
- Share of M's pairs above K·q^(-2/3): K = 2: 8.4-9.0%; K = 3: 1.4-1.7%; K = 4: 0.14-0.24%;
  K = 6: about 0; stable in q. LP load shares are the same.
- E sum load^2 d^2: 1.15, 0.84, 0.72, 0.62, 0.56. E|grad D|^2: 0.57, 0.65, 0.65, 0.64, 0.61
  (measured Var D is 0.003).
- Theorem B logic test (test_chain.py; q = 16, 24, 32; 80 instances; four thresholds):
  507,874 price inequalities, 0 violations; 5,476 closing cases, 0 violations.
- Theorem A: worst 2-swap repair cost × q^(1/4) = 2.87, 2.86, 2.86, 2.95, 2.95 (constant,
  as predicted for k = 2); the cheapest rotation is always >= d(a,b). At q <= 48 the
  cheapest rotation (about 1.4) is far above the largest used deficit (about 0.4): the
  theorems are asymptotic and their constants are not visible at these sizes.
