# Vanishing of the core on the random triangle market: status (2026-11-11)

Written by a Claude Code lane agent in a 60-minute time box. Not reviewed; not reproduced by a
different model. Floating-point numerics are evidence, not certificates. Files: `vanishing/`.

## Target statement (OPEN, not proved here)

**Conjecture V.** Three sectors X, Y, Z of q agents; v(x,y,z) = A[x,y] + B[x,z] + C[y,z] with all
3q^2 entries iid U[0,1]. Then P(core nonempty) -> 0 as q -> infinity.

Equivalently (draft Lemma 2 / PROOF_NOTES Theorem 1, LP formulation): P(LP* = T) -> 0, where LP* is
the optimum of the fractional axial 3-index assignment LP and T the integer optimum.

**Outcome of this session:** neither (1) nor a weaker complete limit theorem. Obtained: three
complete elementary results (Lemmas 1-2, Proposition 3), a precise reduction of Conjecture V to a
concentration statement about one scalar (Reduction 4), and numerical evidence that changes the
proof strategy: both the "global" route and the bounded-size "local" route are contradicted by the
data. Gaps are listed at the end.

---

## Complete results

Throughout, M denotes an optimal integer assignment (a.s. unique, a perfect set of q disjoint
teams; with v >= 0 a partial optimum extends to a perfect one at no loss). P = the axial polytope
{lam >= 0 : every agent has load exactly 1} (padding, PROOF_NOTES Thm 1 Step 1, shows LP* is the
max over P). For S a set of teams of M, G_S is the sub-market on the 3|S| agents of S.

**Lemma 1 (localization). Status: complete.**
(i) M restricted to S is an integer optimum of G_S. (ii) LP*(G) - T(G) >= LP*(G_S) - T(G_S).
(iii) The core of G is empty iff some S has empty sub-core (S = M gives "only if").
*Proof.* (i) If a perfect assignment M'_S of G_S beat M|S, then M'_S together with M|(M\S) would beat
M. (ii) Let lam_S be an LP optimum of G_S; lam_S + (indicator of M\S) is feasible for G with value
LP*(G_S) + T(G) - T(G_S). (iii) "If" by (ii) and the LP criterion. QED.

**Lemma 2 (stability margin). Status: complete.** For a vertex w of P put
d(w) = sum over teams e of M of (1 - w_e) (> 0 for w != M). Define
mu* = max over vertices w != M of v.(w - M) / d(w).
Then M is an optimal vertex of the LP with objective v + mu 1_M iff mu >= mu*, and
core(G) is empty iff mu* > 0.
*Proof.* (v + mu 1_M).w = v.w + mu (q - d(w)) and (v + mu 1_M).M = v.M + mu q, so M is optimal iff
v.(w - M) <= mu d(w) for every vertex w, i.e. mu >= mu*. Core empty iff LP* > T iff some vertex w
has v.w > v.M iff mu* > 0. QED.
(The maximizing w can be taken adjacent to M on P: standard simplex fact, applied to v + mu* 1_M.)

**Proposition 3 (P(empty) > 0 for every q >= 3). Status: complete (gives no uniformity in q).**
For every q >= 3, P(core empty) >= (1/(200 q))^{3 q^2} > 0.
*Proof.* Let (A3, B3, C3) be the exact q = 3 certificate idx = 54 in
`gate1_ext/out/certificates_q3.txt` (entries in [0,1]; re-checked in exact rationals by
`vanishing/verify_cert54.py`: a half-integral lam on six teams has value 113/20, the best of the 36
perfect assignments is 559/100, gap 3/50). Build a q-market: agents 0,1,2 of each sector carry
A3, B3, C3; for i, j >= 3 put A[i,j] = B[i,j] = C[i,j] = 1 if i = j, else 0; every pair with one
agent in {0,1,2} and the other in {3,...,q-1} gets 0.
Upper bound on T: split the teams of any assignment into pure-block, pure-rest and mixed. A mixed
team has exactly two of its three pairs crossing block/rest (any 2-colouring of three members that
is not constant has exactly two bichromatic pairs), so its value is <= 1. Pure-rest teams have value
<= 3; pure-block teams form an assignment of the q = 3 market, total <= T3 = 559/100. If m teams are
mixed, some sector has >= ceil(m/3) rest agents in mixed teams, so at most q - 3 - ceil(m/3)
pure-rest teams exist. Total <= T3 + 3(q - 3) - 3 ceil(m/3) + m <= T3 + 3(q - 3), attained by
M3 plus the rest diagonal. Lower bound on LP*: lam3 plus the rest diagonal gives LP3 + 3(q - 3).
So LP* - T >= 3/50. Both LP* and T are maxima over sets of total team weight q of sums of three
pair entries, so each moves by at most 3q delta when every pair entry moves by at most delta;
LP* - T moves by at most 6 q delta, so for delta <= 1/(200 q) the gap stays >= 3/50 - 3/100 = 3/100 > 0. The
delta-box around the construction intersected with [0,1]^{3q^2} has volume >= delta^{3q^2}. QED.

---

## Numerical findings that redirect the strategy

Gate 1 seeds throughout (`SeedSequence([20261111, 0, 0, q, idx])`), scipy HiGHS, 6 workers.

### F1. The global route is refuted (checked, `vanishing/global_check.py`, `out/global_check.json`)

S = sum of the three bipartite optima (an upper bound on LP*). 300 instances per q (T only where
stored, i.e. all at q = 24, nonempty-only at q >= 28).

| q | 3q - S | 3q - LP* | 3q - T | S - LP* | LP* - T (q=24, all) |
|---|---|---|---|---|---|
| 24 | 4.54 | 7.39 | 7.46 | 2.84 | 0.071 |
| 28 | 4.64 | 7.83 | - | 3.19 | - |
| 32 | 4.59 | 8.10 | - | 3.51 | - |
| 40 | 4.67 | 8.75 | - | 4.07 | - |
| 48 | 4.72 | 9.22 | - | 4.50 | - |

3q - S tends to 3 zeta(2) ~ 4.93 (bipartite Aldous limit, U[0,1] near 1); 3q - LP* grows like
~2.55 q^{1/3} (ratio to q^{1/3}: 2.56, 2.58, 2.55, 2.56, 2.54). The LP pays the full q^{1/3}
integrality-type loss; the LP/IP gap is O(1). The coordinator's reading is confirmed: discard.

Interpretation (heuristic, not a proof): with eps ~ q^{-2/3}, the pairs with value >= 1 - eps form a
random tripartite graph of degree ~ q^{1/3}; a near-optimal assignment is a triangle factor of it.
The threshold for a triangle factor and for a *fractional* triangle factor are both governed by
"every vertex lies in a triangle" (Johansson-Kahn-Vu scale n^{-2/3} (log n)^{1/3}), which is why LP*
and T lose the same leading q^{1/3} term. The core question is a second-order one.

### F2. The bounded-size local route is contradicted (`vanishing/local_count.py`, `out/local_q*.jsonl`)

N3 (N4) = number of 3-subsets (4-subsets) S of M whose sub-market has empty core (Lemma 1). N3 is
exact against the full list of extreme points of the 3x3x3 marginal image (`verts3.py`: 1,386
vertices of the axial polytope, 1,116 extreme marginal images, of which 1,080 fractional: 432 with
weight 1/2 on six teams, 648 with weights 1/3 x5 and 2/3 x2; cdd in floating point). N4 by one LP
per subset, only for q <= 12. k* = number of M-teams disrupted by the LP-optimal vertex.

| q | n | empty | N3 > 0 | N4 > 0 (when N3 = 0) | k* median / mean / max | LP* - T median |
|---|---|---|---|---|---|---|
| 6 | 150 | 22 | 1 | 10 | 5 / 4.95 / 6 | 0.047 |
| 8 | 60 | 16 | 1 | 3 | 6 / 6.12 / 8 | 0.046 |
| 10 | 150 | 43 | 1 | 2 | 8 / 7.79 / 10 | 0.045 |
| 12 | 150 | 62 | 1 | 4 | 9 / 9.02 / 12 | 0.046 |
| 16 | 150 | 97 | 1 | not run | 12 / 11.48 / 16 | 0.056 |
| 20 | 150 | 117 | 0 | not run | 16 / 15.27 / 20 | 0.062 |
| 24 | 150 | 133 | 0 | not run | 20 / 18.89 / 24 | 0.071 |

Sanity: no nonempty instance had N3 > 0; no 3-subset of M was beaten by an integer reassignment
(max violation 3e-15). The most *efficient* improving edge (the maximizer in Lemma 2, found by
bisection on mu, `efficient_move.py`, `out/eff_q*.jsonl`) is also mesoscopic: its number of
disrupted M-teams has median about 5 (q=6, 14 empty), 8 (q=10, 32), 12 (q=16, 68), 19 (q=24, 90). Weights are mostly 1/2
(half-integral "odd-cycle" structures with 2k teams on k M-teams), with some 1/3 and 1/4.

Reading: improving structures supported on a bounded number of M-teams essentially never occur, and
their frequency falls with q; the obstruction lives on a positive fraction (~0.75-0.8) of the
market. A second-moment count of bounded local half-integral configurations relative to M cannot
work as stated, because the expected number of such configurations appears to tend to 0, not
infinity.

Heuristic explanation (not a proof): conditional on M, the M-pairs sit within ~q^{-2/3} of 1. The
q = 3 half-integral structure replaces 4 units of M pair-weight by 8 new pairs at weight 1/2, so to
improve it needs all 8 new pairs within O(q^{-2/3}) of 1: probability ~ q^{-16/3} per triple of
M-teams against ~q^3 triples, expected count -> 0. A structure on k M-teams using r new pairs has
expected count ~ q^{k - 2r/3}; local structures have r/k too large. Mesoscopic structures win by
exploiting the q^{-2/3} deficits spread over many M-teams, so the relevant object is global-in-M.

### F3. The stability margin concentrates (`vanishing/margin.py`, `out/margin_q*.jsonl`)

Signed mu* of Lemma 2, by bisection, 120 instances per q:

| q | P(mu* > 0) | mean mu* | sd mu* | mean/sd | 10% | 50% | 90% |
|---|---|---|---|---|---|---|---|
| 4 | 0.042 | -0.0644 | 0.0560 | -1.15 | -0.146 | -0.051 | -0.005 |
| 6 | 0.150 | -0.0294 | 0.0402 | -0.73 | -0.087 | -0.016 | 0.010 |
| 10 | 0.292 | -0.0090 | 0.0214 | -0.42 | -0.046 | -0.005 | 0.016 |
| 16 | 0.683 | 0.0052 | 0.0120 | 0.43 | -0.008 | 0.005 | 0.021 |
| 24 | 0.892 | 0.0075 | 0.0060 | 1.24 | -0.0001 | 0.0075 | 0.015 |

(P(mu* > 0) agrees with the gate-1 empty shares 0.030, 0.108, 0.331, 0.634, 0.878.) The mean turns
positive between q = 10 and 16 and stays near 0.005-0.01, while the sd falls roughly like
q^{-1.4} (10 -> 24). Small samples; the exponent is a rough reading.

---

## Reduction 4 (precise; the lemma is OPEN)

**Lemma C (open).** There are deterministic m_q > 0 with Var(mu*) = o(m_q^2) and
E mu* >= m_q for all large q.

**Claim.** Lemma C implies Conjecture V. *Proof (complete given Lemma C):* by Lemma 2,
P(core nonempty) = P(mu* <= 0) <= P(|mu* - E mu*| >= m_q) <= Var(mu*)/m_q^2 -> 0 (Chebyshev).

Support: F3 (mean/sd rises from -1.15 to +1.24 over q = 4..24, sd shrinking faster than the mean).
What a proof of Lemma C needs, and why it is hard:
- *Positivity of the mean.* A lower bound mu* >= v.(w - M)/d(w) for one explicit w built from M;
  F2 says w must be mesoscopic (k ~ 0.8 q M-teams), so this is a structural statement about the
  random optimum, not a local count.
- *Variance.* Efron-Stein over the 3q^2 entries needs each entry's resampling influence on mu* to be
  o(m_q / q) on average. The maximizer is a ratio over vertices with d(w) >= 1/2 (weights are
  multiples of 1/k), so a single entry's naive influence is O(1); the needed bound is that a random
  entry lies on the support of M or of the maximizing w with probability O(1/q) and that, when it
  does, the change is small. Not attempted.

Alternative reductions considered and set aside: the LP* - T gap has an atom at 0 and a continuous
part of order 0.05-0.07 roughly constant in q, so a concentration statement for the gap itself
cannot separate it from 0; mu* is the natural scalar because it is signed.

---

## Gaps (strict list)

1. Conjecture V itself: open. No limit theorem was proved (not even limsup P(nonempty) < 1).
2. Lemma C: open; supported only by 600 instances at q <= 24 (F3).
3. The vertex list in F2 came from cdd in floating point (weights 1/3 printed as floats); counts
   1,386 / 1,116 / 1,080 are not certified. N3 conclusions rely on it; N4 used direct LPs.
4. The "triangle-factor threshold" explanation of F1 and the q^{k - 2r/3} count in F2 are
   heuristics, not proofs.
5. Uniform lower bound P(empty) >= c > 0 for all q >= 3: not proved (Proposition 3 decays in q).
   Planting a bounded bad block does not give it: a bounded set of entries can only be pushed
   within delta of a target at probability ~delta^{O(1)}, and F2 shows bounded blocks stop being
   used by M once q is large.

## Suggested next step

Study the maximizer of Lemma 2 directly: its support graph on M (which M-teams, which new pairs,
odd-cycle structure), whether it decomposes into many nearly disjoint pieces of size growing with q
(then a second moment over *mesoscopic* pieces, e.g. of size q^a, could work), and the per-entry
influence on mu* (for Efron-Stein). Both are cheap numerically with `efficient_move.py`.
