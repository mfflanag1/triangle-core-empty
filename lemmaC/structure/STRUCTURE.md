# Structure of the improving object on the random triangle market (2026-11-11)

Lane agent (Claude Code), 60-minute time box. Not reviewed, not reproduced by another model.
**All numbers are floating point (scipy HiGHS dual simplex + scipy milp). They are evidence, not certificates.**
Seeds: gate-1 `instance('triangle','unif',q,idx)`, idx = 0..n-1, so instances match `vanishing/` and `gate1_ext/`.

Files (all in this folder):
- `mustruct.py`: per-instance structure run -> `out/struct_q{08..48}.jsonl` (n = 200 at q = 8..40, n = 120 at q = 48).
- `blocks.py`: block-superadditivity probe -> `out/blocks_q{16,24,32}.jsonl` (n = 40 each).
- `influence.py`: Efron-Stein resampling probe -> `out/influence_q*.jsonl` (n = 16, K = 12 per stratum) and
  `out/influence2_q*.jsonl` (n = 12, K = 10, adds T and LP*).
- `analyse.py`: prints every table below (`out/analysis_partial.txt` is an earlier snapshot).

## 0. Two method changes that are themselves useful

**(a) mu\* is one LP, not a bisection.** By Lemma 2 and the tangent cone of P at M:

    mu* = max  v.y - v_M.z   s.t.  for each agent i:  sum_{f containing i, f not in M} y_f = z_{M(i)},
                                    sum_{e in M} z_e = 1,   y, z >= 0.

Dual: **mu\* = min mu such that some free price vector p has p(f) >= v_f for every team f not in M and
p(e) <= v_e + mu for every e in M.** So mu\* is the smallest uniform subsidy to M's teams that makes M
LP-optimal. An optimal basic solution is an extreme ray, which gives the adjacent vertex w = M + (y - z)/max z.
Checks: it agrees with the bisection values in `vanishing/out/margin_q24.jsonl` to 4e-8 (120 instances) and with
primal/dual/ratio to 1e-14. Max dual-feasibility violation is 2.5e-14.

**(b) Exact integer optimum by reduced-cost fixing.** Take the full-LP dual p_lp and r_f = v_f - p_lp(f) <= 0.
Every perfect x has v.x = LP\* + sum r_f x_f, so an optimal team has r_f >= T_inc - LP\*. Run the MILP on
{r_f >= -delta} and double delta until T_inc >= LP\* - delta. The candidate set is 1.4-8% of the q^3 teams.
It agrees with the full-MILP T of `gate1_ext` at q = 24 (200 instances, max difference 3e-9). This is what made
q = 32-48 cheap enough to run (13 s per instance at q = 48).

## 1. Scalars and scaling (item 4)

| q | n | P(mu\*>0) | mean mu\* | sd mu\* | mean/sd | mean gap LP\*-T | sd gap | mean/sd gap | mean 3q-T | sd 3q-T | mean 3q-LP\* | M pair deficit x q^(2/3) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 8 | 200 | 0.255 | -0.0119 | 0.0226 | -0.53 | 0.011 | 0.026 | 0.42 | 4.89 | 0.60 | 4.88 | 0.82 |
| 12 | 200 | 0.385 | -0.0035 | 0.0160 | -0.22 | 0.022 | 0.036 | 0.60 | 5.89 | 0.58 | 5.87 | 0.86 |
| 16 | 200 | 0.645 | 0.0043 | 0.0111 | 0.38 | 0.039 | 0.046 | 0.84 | 6.44 | 0.52 | 6.40 | 0.85 |
| 24 | 200 | 0.875 | 0.0076 | 0.0065 | 1.17 | 0.069 | 0.053 | 1.30 | 7.48 | 0.48 | 7.41 | 0.86 |
| 32 | 200 | 0.965 | 0.0095 | 0.0050 | 1.90 | 0.108 | 0.054 | 1.99 | 8.21 | 0.46 | 8.10 | 0.86 |
| 40 | 200 | 0.985 | 0.0090 | 0.0038 | 2.39 | 0.127 | 0.053 | 2.41 | 8.90 | 0.46 | 8.77 | 0.87 |
| 48 | 120 | 1.000 | 0.0088 | 0.0030 | 2.92 | 0.156 | 0.054 | 2.91 | 9.35 | 0.43 | 9.19 | 0.86 |

Fits (log-log over q = 16..48, 90% bootstrap bands):
- sd mu\* ~ q^-1.17 [-1.30, -1.04]. **mean mu\* plateaus near 0.009 for q >= 32.** mean/sd of mu\* ~ q^1.82 [1.55, 2.19].
  Quantiles at q = 48: 5% 0.0040, 50% 0.0088, 90% 0.0125.
- 3q - T ~ q^0.341 [0.331, 0.349], ratio to q^(1/3) 2.56-2.60. 3q - LP\* ~ q^0.332, ratio 2.53-2.57. sd(3q - T) is flat
  at about 0.45. The M pair deficit (mean of 1 - pair value over M's 3q pairs) x q^(2/3) is 0.86, flat.
- **The gap LP\* - T is not O(1).** Its mean grows (slope 1.27 over 16..48, still above 1 at the top end), while its sd
  is flat at 0.053-0.054 from q = 24 on. mean/sd of the gap ~ q^1.15 [1.03, 1.27]. gap/q^(1/3) = 0.015, 0.024,
  0.034, 0.037, 0.043 at q = 16, 24, 32, 40, 48, still rising.
  This corrects THEOREM_VANISHING.md ("gap O(1), cannot separate from 0"), which rested on the q = 24 value alone.
  The gap itself is a usable concentrating scalar.

## 2. Decomposition of the improving structure (item 1)

**mu\*-maximiser w** (adjacent vertex, all instances):

| q | one component | kM (M-teams touched) mean / median | kM/q | off-M teams / kM | share of off teams keeping an M-pair | new pairs / kM | aux graph non-bipartite | w half-integral |
|---|---|---|---|---|---|---|---|---|
| 8 | 1.00 | 5.0 / 5 | 0.63 | 1.16 | 0.67 | 2.62 | 0.84 | 0.35 |
| 12 | 1.00 | 7.5 / 8 | 0.63 | 1.26 | 0.57 | 2.94 | 0.91 | 0.35 |
| 16 | 1.00 | 10.2 / 11 | 0.64 | 1.30 | 0.51 | 3.11 | 0.94 | 0.40 |
| 24 | 1.00 | 17.3 / 18 | 0.72 | 1.48 | 0.44 | 3.64 | 0.99 | 0.14 |
| 32 | 1.00 | 24.2 / 26 | 0.76 | 1.61 | 0.39 | 4.00 | 1.00 | 0.07 |
| 40 | 1.00 | 31.6 / 32 | 0.79 | 1.72 | 0.37 | 4.29 | 1.00 | 0.02 |
| 48 | 1.00 | 39.4 / 40 | 0.82 | 1.79 | 0.33 | 4.51 | 1.00 | 0.01 |

In empty-core instances, d(w)/q = 0.48, 0.41, 0.37, 0.33, 0.29, 0.28, 0.26 and the gain v.(w - M) is 0.04 -> 0.11.

**LP-optimal vertex lam\* - M** (empty-core instances):

| q | n_empty | components | kM/q | d(lam\*)/q (M-weight moved) | off-M teams / q | keep-M-pair share | lam\* half-integral |
|---|---|---|---|---|---|---|---|
| 8 | 51 | 1 in 51/51 | 0.80 | 0.49 | 1.07 | 0.62 | 0.65 |
| 12 | 77 | 1 in 77/77 | 0.76 | 0.44 | 1.10 | 0.54 | 0.44 |
| 16 | 129 | 1 in 129/129 | 0.72 | 0.42 | 1.04 | 0.50 | 0.39 |
| 24 | 175 | 1 in 175/175 | 0.79 | 0.42 | 1.25 | 0.44 | 0.15 |
| 32 | 193 | 1 in 193/193 | 0.84 | 0.43 | 1.49 | 0.39 | 0.02 |
| 40 | 197 | 1 in 197/197 | 0.85 | 0.44 | 1.59 | 0.36 | 0.02 |
| 48 | 120 | 1 in 120/120 | 0.89 | 0.46 | 1.73 | 0.33 | 0.00 |

- **The improving structure is always one giant component**, not many pieces. The support of lam\* - M is
  agent-connected in all 942 empty instances. It touches 72-89% of M's teams, and kM/q rises with q.
  It moves a stable positive density of M-weight (d/q of about 0.43). For the maximiser, "the best single component"
  is the whole structure.
- **Block superadditivity fails empirically.** Lemma 1 gives gap(G) >= sum_j gap(G_Sj) for any partition of M into
  blocks. For random partitions of M into blocks of size s in {4, 6, 8, 12, q/4, q/2} (3 partitions per instance,
  40 instances at each of q = 16, 24, 32), **not one block had a positive sub-gap** (0 of every block count), although
  31, 31 and 39 of those instances have empty cores. A random half of M never carries an improving structure.

## 3. Combinatorial pattern (item 2)

- In the auxiliary graph (nodes = M-teams; each off-M support team links the M-teams of its 3 agents), the 2-section
  is non-bipartite in 100% of LP-optimal structures and in 94-100% of maximisers at q >= 16. So odd cycles are
  present, as any fractional vertex needs.
- **The vertices are not half-integral.** The half-integral share falls from 0.65 (q = 8) to 0.00 (q = 48) for lam\*,
  and from 0.40 (q = 16) to 0.01 (q = 48) for w. Denominators of 3-24 are common, and at q >= 32 most vertices have
  some weight with denominator > 24 (e.g. 1/83 seen at q = 32). The "half-integral odd cycle of length 2k on k
  M-teams" picture from F2 is a small-q phenomenon only.
- Off-M support teams are of two kinds:
  - "keep an M-pair" teams: two agents from the same M-team, one swapped in. Their share falls from 0.6 to 0.33.
  - "fully new" teams: 3 distinct M-teams.

  The number of off-M teams per touched M-team rises (1.1 -> 1.8), and new (non-M) pairs per touched M-team rise
  (2.6 -> 4.5). The structure becomes denser and more "fully new" as q grows.

## 4. Reduced costs and the near-tight hypergraph (item 3)

- On supp(w) \ M the reduced cost w.r.t. the mu\*-LP dual is 0 (complementary slackness; max |r| 5e-15). Under the
  mu\*-dual, **94-96% of M-teams sit exactly at p(e) = v_e + mu\***: the subsidy is binding almost uniformly over M.
- Under the full-LP vertex dual p_lp, only 9-14% of M-teams carry negative reduced cost (the slack summing to the gap).
- Near-tight off-M teams per agent, N(eps) = #{f not in M containing the agent : v_f - p_lp(f) > -eps}:

| q | 0.005 | 0.01 | 0.02 | 0.05 | 0.1 | 0.2 | eps = 0.5 q^(-2/3) | eps = q^(-2/3) | min over agents at q^(-2/3) |
|---|---|---|---|---|---|---|---|---|---|
| 8 | 1.54 | 1.60 | 1.73 | 2.12 | 2.87 | 4.87 | 3.35 | 6.06 | 1.30 |
| 16 | 1.90 | 1.99 | 2.19 | 2.86 | 4.28 | 8.19 | 3.62 | 6.34 | 0.89 |
| 24 | 2.14 | 2.27 | 2.55 | 3.52 | 5.63 | 11.78 | 3.90 | 6.64 | 0.74 |
| 32 | 2.19 | 2.34 | 2.67 | 3.94 | 6.64 | 15.06 | 3.92 | 6.59 | 0.67 |
| 40 | 2.20 | 2.38 | 2.79 | 4.27 | 7.67 | 18.66 | 3.86 | 6.55 | 0.54 |
| 48 | 2.25 | 2.47 | 2.93 | 4.64 | 8.70 | 22.34 | 3.88 | 6.53 | 0.57 |

  At eps = c q^(-2/3) the count is **q-independent** (slopes 0.02-0.05): about 3.9 at c = 1/2 and 6.6 at c = 1.
  The near-tight hypergraph is sparse with bounded degree, so a local weak limit is plausible. At fixed eps the excess
  over the LP support (about 2 per agent) fits N(eps) - N(0.005) of about 4.5 q eps^(3/2). The ratio is 4.3-5.0 for
  eps in [0.05, 0.2] and q >= 16, which is consistent with a single scaling function of eps q^(2/3).
  The expected "q^(1/3) per agent" holds only for eps of order q^(-4/9) (eps fixed gives about q).
  The minimum over agents falls below 1 at q >= 16: some agents have no near-tight off-M team at all.

## 5. Efron-Stein probe (relevant to the variance half of Lemma C)

Resample single pair entries, recompute, and estimate 0.5 sum_i E(Z - Z^(i))^2. The pairs are stratified as
S1 = pairs used by M or supp(lam\*) and S2 = all other pairs. The Var columns come from the struct runs.
Small samples and heavy tails make the ES estimates noisy to within a factor of about 2-3.

| q | n | ES_T | Var T | ES_LP | Var LP\* | ES_gap | Var gap | (mean gap)^2 | ES_mu | Var mu\* | (mean mu\*)^2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 12 | 12 | 0.54 | 0.34 | 0.46 | 0.33 | 0.082 | 0.0013 | 0.00046 | 0.0081 | 2.6e-4 | 1.2e-5 |
| 16 | 12 | 1.88 | 0.27 | 1.66 | 0.26 | 0.123 | 0.0021 | 0.0015 | 0.0064 | 1.2e-4 | 1.8e-5 |
| 24 | 12 | 0.82 | 0.23 | 0.39 | 0.23 | 0.218 | 0.0029 | 0.0048 | 0.0059 | 4.2e-5 | 5.8e-5 |
| 32 | 12 | 0.87 | 0.21 | 0.43 | 0.20 | 0.221 | 0.0030 | 0.0117 | 0.0011 | 2.5e-5 | 9.1e-5 |

ES/mean^2 (the Chebyshev-via-Efron-Stein ratio, which must go below 1): gap 176, 81, 45, 19; mu\* 650, 350, 103, 12
at q = 12, 16, 24, 32. A separate run (`influence_q*`, n = 16, K = 12) gives ES_gap of 0.090, 0.128 and 0.230 at
q = 12, 16, 24, consistent with these values.

- **First-order Efron-Stein currently exceeds mean^2 for both the gap and mu\*, but the ratio is falling fast.**
  ES_gap is 40-80 x Var(gap): the entry influence on the gap is O(0.03-0.05) whenever the entry touches the support
  (about 6q pairs), while the true variance is much smaller. Yet ES_gap levels off (0.22 at q = 24 and 32) while
  mean^2 grows, and ES/mean^2 falls by about q^-3 over 12-32 for the gap and faster for mu\*. Taken at face value,
  the extrapolation crosses 1 near q of about 40-50 for mu\* and about 60-80 for the gap. With 12 instances per q this
  is a weak reading. Plain Efron-Stein for Lemma C is therefore **not ruled out**, but it is not supported either.
- For T and LP\* separately, ES is O(1) (0.4-1.9, no growth) against Var of about 0.2-0.34. Measured against
  the q^(1/3) scale of the means this is negligible.

## 6. Proposed statement to prove (item 5)

Facts 2-5 rule out every route that assembles the gap from many small or independent pieces:
- bounded local structures, which die (F2);
- half-integral odd cycles, which disappear by q of about 32;
- unions of disjoint mesoscopic pieces: there is always one component;
- block superadditivity over partitions of M: random halves of M never have empty sub-cores;
- plain Efron-Stein on mu\* or on the gap is not ruled out, but it is still above mean^2 by a factor of 12-19 at q = 32 (Section 5).

What is stable:
- the LP optimum differs from M on a **positive density** of agents (d/q of about 0.43, kM/q of about 0.85);
- the near-tight hypergraph at scale q^(-2/3) has a **q-independent bounded degree**;
- T and LP\* each have O(1) fluctuations around means of order 2.55 q^(1/3);
- gap/q^(1/3) rises steadily from 0.015 to 0.043.

**Proposed statement (Lemma G, two constants).** Write 3q - T = q^(1/3) a_IP(q) and 3q - LP\* = q^(1/3) a_LP(q).
Then
(i) Var(T) = o(q^(2/3)) and Var(LP\*) = o(q^(2/3)), and
(ii) liminf_q (E a_IP(q) - E a_LP(q)) > 0, i.e. E(LP\* - T) >= c q^(1/3).

Lemma G implies Conjecture V. P(LP\* = T) <= P(|T - ET| >= c q^(1/3)/2) + P(|LP\* - E LP\*| >= c q^(1/3)/2) -> 0 by
Chebyshev on T and LP\* **separately**. This avoids concentrating the gap or mu\*, whose Efron-Stein bounds are still far above mean^2 at the q measured.

Why this is the most promising formulation:
1. Part (i) is the standard kind of statement: data give sd of about 0.45 and ES of order 1, so Efron-Stein with an
   O(q^(-2/3))-to-O(q^(-1/3)) rerouting cost per used pair should give Var = O(q^(1/3)) or better.
2. Part (ii) is a statement about the **local weak limit**. Rescale values near 1 by q^(2/3). The near-tight
   hypergraph converges, with bounded degrees (Section 4), to a Poisson-weighted 3-partite hypertree (the analogue
   of Aldous's PWIT). Then a_IP and a_LP become the per-agent costs of the optimal integral and fractional perfect
   triangle factors on that limit object. A positive-density disagreement between lam\* and M (Section 2) with
   per-agent gain of order q^(-2/3) is exactly what a_LP < a_IP predicts. A lower-order gap would rather predict a
   vanishing disagreement density. The question becomes "the fractional and integral recursive distributional
   fixed points differ", which the cavity/objective method can address (Aldous; Wastlund-style LP arguments), and
   which can be checked numerically by population dynamics before any proof is attempted.
3. Fallback if (ii) fails (a_LP = a_IP, gap of lower order). Then the gap or mu\* has to be concentrated directly.
   Two options:
   - Lemma C on mu\*: mean mu\* plateaus near 0.009 while sd falls like q^(-1.2), and ES_mu/mean^2 falls 650 -> 12
     over q = 12..32.
   - Lemma C for the gap: E gap grows and Var(gap) is flat at about 0.0029.

   In both, the work is a bound on the entry influence restricted to the about 6q support pairs. The influence would
   need to average o(mean / sqrt(q)) over those pairs, i.e. O(0.03-0.05) per pair would have to shrink with q.

**Caveats.** q <= 48, 120-200 instances per q, all floating point. gap/q^(1/3) is still rising at q = 48, so
"converges to a positive constant" and "grows" cannot be told apart, and neither can "is o(1) eventually".
The ES numbers rest on 7-16 instances and 10-12 resamples per stratum. The local-weak-limit picture is a heuristic.
The 3q - T ~ a q^(1/3) law is empirical: it is not proved here, nor (to my knowledge) for pair-additive values.
