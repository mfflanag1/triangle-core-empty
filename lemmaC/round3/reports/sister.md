# Round 1 report: sister line

Saved by the orchestrator from the agent's message (the agent's own write of a report
file was refused by the harness). Scripts and data: `../sister/` (iid_axial.py,
iid_sum3.py, out/*.jsonl). Nothing here is verified.

## 1. Proved (claimed complete by the agent)

**Theorem S.** Model R(n,r): parts A, B, C of n vertices, r points per vertex, two
independent uniform bijections (A-points to B-points, A-points to C-points); one edge per
A-point. This is an r-regular 3-partite 3-uniform multi-hypergraph. nu = largest matching,
nu* = largest fractional matching.

Statement: nu* = n always. For integer m <= n, mu = m/n:

  P(nu >= m) <= (rn+1)^2 exp(n g_r(mu)),
  g_r(mu) = 3H(mu) + 3 mu ln r - 2r H(mu/r)   (H = natural-log binary entropy).

g_r(1) = ln r - 2(r-1) ln(r/(r-1)) = -0.693, -0.523, -0.340, -0.176, -0.032 for r = 2..6,
positive from r = 7 (increasing in r). So for r <= 6: nu* - nu >= (1 - mu_r) n - 1 with
probability 1 - exp(-Omega(n)), where mu_r is the largest root of g_r:
1 - mu_r = 0.0819, 0.0528, 0.0292, 0.0125, 0.0016 (roots located numerically; hand check
g_2(0.95) = -0.197).

Proof: weight 1/r on every edge loads each vertex exactly 1, so nu* = n. X_m = number of
m-sets of A-points forming a matching. There are C(n,m) r^m m-sets on distinct A-vertices;
the image of a fixed m-set under a uniform bijection is a uniform m-subset of the rn
B-points, on distinct vertices with probability C(n,m) r^m / C(rn,m); same independently
for C. So E X_m = (C(n,m) r^m)^3 / C(rn,m)^2. Use C(n,m) <= exp(nH(mu)),
C(rn,m) >= exp(rn H(mu/r))/(rn+1), Markov, and {nu >= m} = {X_m >= 1}. QED.

**Corollary.** Edge values in [1-delta, 1], any joint law: LP >= (1-delta)n, IP <= nu, so
LP - IP >= (1 - delta - mu)n.

Sanity (exact MILP): r=2, n=60, 10 instances: nu/n mean 0.89, max 0.90 (bound 0.918).
r=3, n=40, 6 instances: mean 0.929, max 0.95 (bound 0.947 is asymptotic).

**Proposition A**, model (a) = n^3 iid Exp(1) or U[0,1] team costs:
P(min assignment cost <= c/n) <= (n!)^2 (c/n)^n / n! <= e sqrt(n) (c/e)^n.
So n·IP >= e - o(1) with probability 1 - exp(-Omega(n)).

GAP: a fractional solution of cost (e - eps)/n in model (a). That alone would prove the
constant gap there. Numerically LP is about 2.30/n, 15% of room.

## 2. Numerics (floating point)

Model (a), costs scaled by n^2, per agent; 300 instances per size (104 at n=16):

| n | 3 | 4 | 5 | 6 | 8 | 10 | 12 | 16 |
|---|---|---|---|---|---|---|---|---|
| P(LP integral) | .743 | .623 | .573 | .387 | .220 | .117 | .067 | .010 |
| LP/n | 1.95 | 2.05 | 2.06 | 2.19 | 2.34 | 2.30 | 2.28 | 2.30 |
| IP/n | 2.07 | 2.20 | 2.20 | 2.40 | 2.60 | 2.59 | 2.58 | 2.65 |
| gap/n | .116 | .146 | .138 | .210 | .263 | .282 | .298 | .350 |

Model (a3): q^3 iid team deficits, each a sum of three fresh uniforms (same one-team law,
same scale and same union bound 2.149 as the pair-additive model; pair-sharing removed).
Deficits divided by q^(1/3); thr = finite-q first-moment threshold:

| q | 8 | 12 | 16 | 32 | 40 |
|---|---|---|---|---|---|
| instances | 200 | 200 | 200 | 100 | 100 |
| P(LP integral) | .260 | .075 | .015 | 0 | 0 |
| LP deficit | 2.116 | 2.125 | 2.154 | 2.192 | 2.167 |
| IP deficit | 2.175 | 2.202 | 2.244 | - | - |
| E gap | .118 | .178 | .227 | - | - |
| thr | 2.026 | 2.055 | 2.072 | 2.103 | 2.110 |

(se of LP deficit about 0.01; sd of gap about 0.11; thr limit 2.149; q = 32, 40 are LP only.)
Pair-additive values for comparison: P(LP integral) = .355 at q=16, .125 at q=24;
E D = .039 at q=16; LP deficit about 2.53, IP about 2.57.

## 3. Sources opened (arXiv abstract pages only)

- Frieze, Sorkin, arXiv:1008.0390: for the axial problem "even the scaling is not known";
  cost Omega(1/n); an algorithm reaches O(n^(-1+o(1))).
- Martin, Mezard, Rivoire, arXiv:cond-mat/0507180 and cond-mat/0407623: cavity method; for
  d > 2 a completely frozen glass phase with vanishing entropy.

## 4. Transfer

The proof of Theorem S uses two properties: (R) an explicit fractional point of known
value; (W) an annealed window: the first-moment threshold for integer solutions lies
strictly beyond that fractional value.

(W) holds in R(n,r), r <= 6 (proved); holds numerically in model (a) (LP 2.30 against
e = 2.718); is absent or at most 1% wide in (a3) (LP about 2.17 against 2.149); is absent
in the pair-additive model (2.149 against LP about 2.53). There the union bound sits below
the LP optimum itself, so no fractional point can go under it. A transfer needs an integer
lower bound within about 1.5% of the truth, i.e. a quenched (conditional or truncated)
count instead of the annealed one.

Local structure does transfer: in the pair-additive near-tight hypergraph (deficit
<= t q^(-2/3)) the mean degree is t^3/6, the expected number of pairs of tight teams
through a fixed agent that share a second agent is q^3 s^5/40 = O(q^(-1/3))
(s = t q^(-2/3)), and of Berge triangles through an agent O(1/q). These expectation counts
are complete; local weak convergence to the Poisson hypertree of (a3) is a gap.

Global constants do not transfer at q <= 64: same one-team law and same local limit, yet
(a3) has LP deficit about 2.17 against 2.53 and a gap about six times larger at q=16. Two
readings, not separated by this data: (i) pair-sharing changes the constants through long
cycles; (ii) the pair-additive model converges to the (a3) constants with O(q^(-1/3))
corrections, which are 25-40% at q <= 64. Under (ii) the true pair-additive gap constant
would be the larger (a3) one, and the fit with a negative intercept is pre-asymptotic.

## 5. Confidence and next step

Theorem S, corollary, Proposition A: high. Absence of the window in the pair-additive
model rests on the numerical 2.53. Next step: run the pair-additive model and (a3) on one
interpolating family (a fraction theta of pairs shared) to decide between (i) and (ii).

## Correction from the agent (after reading _lit/LIT_CHECK_1.md)

The iid-model figure "LP about 2.30/n, 15% of room below e" was pre-asymptotic. Extended
LP-only runs (Exp(1) costs, per agent, scaled by n^2):

| n | 16 | 24 | 32 | 48 | 64 | 96 |
|---|---|---|---|---|---|---|
| instances | 141 | 100 | 100 | 60 | 40 | 20 |
| LP/n | 2.308 | 2.419 | 2.452 | 2.478 | 2.499 | 2.463 |
| se | .033 | .030 | .033 | .033 | .030 | .039 |

- The LP constant appears to level off near 2.47-2.50, still below e = 2.718 (about 7
  standard errors at n = 64). The annealed window in the iid model has about 8% of room.
- With 3.126 (Martin-Mezard-Rivoire) as the integer limit, the predicted relative gap in the
  iid model is about 20%, against about 2% in the pair-additive model at q <= 64.
- The lookup's "RS value roughly in [2.72, 2.9]" was read off a figure and does not match
  the LP value 2.5; do not build on it.
- The same-marginal model's LP deficit (2.17-2.19 at q = 32, 40) may still drift.

Nothing changes in Theorem S, its corollary, Proposition A or the transfer verdict.
