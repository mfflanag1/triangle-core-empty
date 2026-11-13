# Red team report (round 3): is E D = o(q^(1/3))?

Status: every line of attack on the conjecture failed, and two of them turned into evidence for it. Nothing here is a proof.
All numerics are floating point (evidence only). Scripts and outputs are in this folder.

## 1. Rounding scheme: dead (evidence, plus one complete first-moment computation)

Idea tested: if the LP optimum were integral except on a vanishing fraction f of agents, re-matching the
fractional part would cost about a·(fq)^(1/3) = o(q^(1/3)).

**Observation (lpfrac.py, own seeds, basic optimal LP solutions, LP only).** Fraction of agents covered
fractionally:

| q | 16 | 24 | 32 | 48 | 64 | 96 | 128 |
|---|---|---|---|---|---|---|---|
| instances | 40 | 40 | 40 | 30 | 30 | 20 | 12 |
| fraction fractional | 0.45 | 0.70 | 0.82 | 0.88 | 0.90 | 0.88 | 0.91 |
| (3q − LP*)/q^(1/3) | 2.63 | 2.59 | 2.55 | 2.56 | 2.50 | 2.51 | 2.49 |
| median smallest positive weight | 0.5 | 0.2 | 0.035 | 0.009 | 0.003 | 0.002 | 0.0014 |

Standard errors: about 0.01 on the fraction (q ≥ 48) and 0.015 to 0.03 on the constant. The fraction
rises to about 0.9 and stays there. Half-integral solutions disappear by q = 48. So the number of
"repairs" is of order q, not o(q): there is no sparse set of odd cycles to break, in contrast with matching
on the complete graph.

**No small exchanges (complete as a first-moment bound).** Call a team c-cheap if its deficit is at most
c·q^(−2/3); P = c³/(6q²) for c ≤ q^(2/3). For a fixed perfect assignment, the expected number of ways to
replace s of its teams by s other c-cheap teams on the same agents is at most
C(q,s)·(s!)²·(c³/6q²)^s ≤ (s·c³/(6q))^s. This tends to 0 uniformly for 2 ≤ s ≤ εq when ε < 6/c³. Any two
perfect assignments made of c-cheap teams therefore differ in at least about 6q/c³ teams (in expectation
count; a w.h.p. statement over all assignments needs a union bound I did not do: **gap**). Rounding moves
of intermediate size do not exist; this extends the brief's "bounded size" obstacle to all sizes below εq.

**Quantitative rounding cost.** By the same count, re-matching a fraction f of the agents together with a
chosen fraction σ of integral teams needs teams of deficit c·q^(−2/3) with
c³ ≥ (6e²/μ²)(σ/e)^(σ/μ), μ = f + σ, so the cost is at least about (6e²)^(1/3)·μ^(1/3)·(σ/e)^(σ/3μ)·q^(1/3).
For f → 0 this can be pushed to f^(2/3−o(1))·q^(1/3); for f ≈ 0.9 it is a constant times q^(1/3)
(first-moment heuristic, **gap**: the LP also pays for those agents, so this is not a lower bound on D).

## 2. The "local" argument for a_LP = a_IP, and why I no longer trust it

Local structure (re-derived): two c-cheap teams at an agent share a second agent with probability
O(q^(−1/3)), so the cheap-team hypergraph converges locally to a Poisson hypertree in which the teams at an
agent have scaled costs forming a Poisson process of intensity t²/2. On a hypertree LP = IP. If the optimum
were a local function (as for the 2-D assignment, and for matching on K_n; see the addendum), a_LP = a_IP
would follow.

**Tree value (rs_rde.py, deterministic grid).** The fixed point of X = min_k(ξ_k − X'_k − X''_k) gives
a_tree = 2.2496 (normalisation check 1.0000; stable to 1e−4 between grid steps 0.02 and 0.01; the same code
returns 1.64495 = π²/6 for the 2-D assignment).

- The first-moment bound (fact 4, re-derived: 3·e^(−1/3) = 2.1496) is BELOW a_tree. So in this model the
  first moment does not refute the local value. For the 3-index model with independent uniform team costs
  the same code gives 2.706 (extrapolated in the grid step; 2.702 and 2.704 at steps 0.02 and 0.01), 0.45%
  below the first-moment bound e = 2.718. See the addendum: this margin is thin and disagrees with a relayed
  footnote of Martin, Mézard and Rivoire, so I do not rely on it.
- But the data are not at the tree value: (3q − LP*)/q^(1/3) is flat near 2.50 for q = 64 to 128, and
  3q − T is near 2.57·q^(1/3). A linear extrapolation in q^(−1/3) of the LP column gives roughly 2.38.
  Either finite-size corrections are large or neither optimum is local.
- **Instability (complete, elementary).** Shifting every input of the tree recursion by δ shifts the
  output by −(d−1)δ exactly (d = 3: factor −2; d = 2: factor −1). The sector shifts with zero sum are a
  symmetry of the dual; the total-price direction is amplified by 2 per level. The fixed point is not an
  attractor of the plain recursion (my population dynamics diverged until damped), so boundary values at
  depth L influence the root at order 2^L. This is a structural difference from the graph case and removes
  the reason to expect the fractional and integral limits to coincide.
- **Entropy crisis in the tree model (rs_finiteT.py; Monte Carlo, population 12,000, 20 averaged sweeps).**
  Bethe free energy f, energy e and entropy s = β(e − f) per team:

  | β | ∞ (grid) | 30 | 12 | 6 | 4 | 3 | 2 |
  |---|---|---|---|---|---|---|---|
  | f | 2.2496 | 2.242 | 2.246 | 2.263 | 2.274 | 2.289 | 2.323 |
  | s (± se) | | −0.40 ± 0.58 | −0.15 ± 0.26 | −0.15 ± 0.14 | −0.21 ± 0.07 | −0.23 ± 0.05 | −0.02 ± 0.03 |

  The free energy RISES with temperature from 2.250 to about 2.32 at β ≈ 2 (noise on f about 0.01), which
  requires negative entropy; the direct estimate is negative at 3 to 5 standard errors for β = 3, 4 and
  returns to zero near β = 2. A set of assignments cannot have negative entropy, so the local (tree) solution
  is inconsistent below T* ≈ 0.5. This is the entropy-crisis signature that, in the statistical physics of
  multi-index matching, gives a frozen phase with integer constant max_β f(β) ≈ 2.32, above the
  zero-temperature tree value 2.250. If the LP constant is the tree value (heuristic: the LP is the
  zero-temperature Bethe minimum), the predicted gap is a_IP − a_LP ≈ 0.07, the same order as the brief's
  fitted 0.11 and as the observed D/q^(1/3) = 0.053 at q = 64 and rising. **Gap:** Monte Carlo only (the
  normalisation check is off by up to 0.9%; correcting for it leaves s negative at β = 3, 4); the
  identification of both constants is physics heuristic; the finite-size constants (2.50, 2.57) are both far
  above the predicted limits (2.25, 2.32).

## 3. Slow growth

Weighted fits of the brief's table to D = a·g(q) − b (fits.py, se = 0.05/√n, 5 dof):

| g(q) | q^(1/2) | q^(1/3) | q^(1/4) | q^(1/6) | log q | q^(1/3)/log q |
|---|---|---|---|---|---|---|
| chi² | 8.0 | 12.4 | 15.3 | 18.7 | 26.9 | 99 |
| predicted D(128) | 0.34 | 0.31 | 0.30 | 0.29 | 0.27 | n/a |
| predicted D(200) | 0.46 | 0.40 | 0.38 | 0.36 | 0.32 | n/a |

q^(1/3)/log q is flat on 16..64 (0.91 to 0.96) and cannot produce a fivefold rise. The data prefer exponents
at or above 1/3; log q is disfavoured (p ≈ 1e−4) and q^(1/6) is poor (p ≈ 0.002), but a two-parameter family
with a free offset is a weak test. Every surviving form has D/q^(1/3) still rising at q = 64, so the
conjecture fails only if that ratio turns over beyond the data.

## 4. What would settle it

- Against the conjecture: integer assignments within 0.31 of LP* at q = 200 (the q = 64 ratio 0.053
  carried over); the q^(1/3) fit predicts 0.40. With sd 0.05, about 10 instances suffice. A heuristic
  integer solution gives an upper bound on D, so this test can be run without exact T.
- For the conjecture: the same experiment showing D(200) ≈ 0.40, or a deterministic (grid) version of the
  finite-temperature tree computation confirming that f(β) has an interior maximum above 2.2496.

## 5. Not done

- Restricted-support integer experiment (LP support plus k teams per agent): not run. With 90% of agents
  fractional and at most 3q − 2 support teams, the support alone rarely contains a perfect assignment.
- Literature: see the addendum.

## Addendum: validation and literature (after the first report to the orchestrator)

**Validation of rs_finiteT.py against a published number.** For independent uniform team costs (d = 3,
Poisson intensity 1) my code gives entropy +0.22 ± 0.01 at β = 0.3, −0.03 ± 0.01 at β = 0.412,
−0.23 ± 0.03 at β = 0.6, −0.33 ± 0.04 at β = 1, and f = 3.20 at β = 0.412 (noise about 0.05), f = 2.95 at
β = 1. Martin, Mézard and Rivoire report the entropy vanishing at β_s = 0.412 ± 0.001, frozen value
3.126 ± 0.002, and f_RS ≈ 2.9 at β = 1 (the last read off their Fig. 6 by the search agent). My zero
crossing is at about 0.40 and my f is 2% high, within Monte Carlo error. So the method reproduces the
published entropy crisis, and the same code gives the crisis in our model (intensity t²/2) at β ≈ 2 with
max f ≈ 2.32 against the zero-temperature value 2.2496.

**Literature (opened by the search agent, not by me; its note is lemmaC/round3/_lit/LIT_CHECK_1.md).**
- Martin, Mézard, Rivoire (arXiv cond-mat/0407623, cond-mat/0507180): independent team costs only; replica
  symmetric equations "are inconsistent when d > 2"; f(β) = f_rs(β_s) for β ≥ β_s; they say densities ~ c^r
  "can be treated along the same lines, but gives different quantitative results"; nothing on LP or
  fractional relaxations; branch and bound at M ≤ 22 gives 3.06 to 3.09 ± 0.03 against the predicted 3.126.
- Unresolved: the agent relays a footnote saying the replica symmetric value violates the annealed bound
  only from d = 4. My grid value for d = 3 uniform costs is 2.706 < e. The margin is 0.45%; one of the two
  is off or they compare different quantities. Nothing in this report depends on it.
- Frieze, Pegden, Tkocz (arXiv 1901.07167) is our model; no LP comparison. No source was found that states
  whether LP and integer constants agree for random multidimensional assignment.
- Thresholds (Kahn, arXiv 1909.06834 and 2008.01605; Devlin and Kahn, arXiv 1703.03513, restating
  Krivelevich): in the random r-graph the fractional and integral perfect matching both appear when isolated
  vertices disappear. In k-out models they can separate: the integer problem needs k exponential in r, the
  fractional one is not known to need more than 2. Thresholds agreeing says nothing about optimisation
  constants at bounded degree, which is our regime (O(1) cheap teams per agent).
- Matching on K_n (Hessler and Wästlund, "friendly model"): an LP relaxation has expected cost within
  O(log n/n) of the integer optimum. This is the d = 2 case, where the recursion's shift factor is −1.

**LP-only run, later sizes.** q = 160 (8 instances): fraction fractional 0.905 ± 0.010,
(3q − LP*)/q^(1/3) = 2.47 ± 0.03. The constant is flat or drifting down slowly (2.50, 2.51, 2.49, 2.47 at
q = 64, 96, 128, 160); the fraction stays at 0.9.

## Assessment

Probability the conjecture is true: about 80% (from 60%). The attempt to refute it failed on every line:
the LP optimum is fractional on 90% of agents, exchanges of fewer than εq teams do not exist, slow-growth
forms fit worse than q^(1/3), and the tree model that would justify a_LP = a_IP shows an entropy crisis with
a predicted gap of about 0.07·q^(1/3) (method validated against the published uniform-cost value). What
remains for the other side: all of this is first-moment, floating-point or physics heuristic, and D/q^(1/3)
has not been seen to level off.

Note: lpfrac.py runs q = 160 and 200 last; it appends to lpfrac.jsonl and `summ.py` prints the table.
