# Typical-case protection in random multi-sided markets: ledger

Append dated entries. Never edit an earlier entry; correct it with a new one.

## 2026-11-09: lane proposed

Proposed by the coordinating Claude Code session as target 1 ("average-case protection")
of the 9 November next-targets slate. Max approved opening it.

## 2026-11-09: prior-art sweep (gate 0), verdict "open only if narrowed"

Full report: `GATE0_REPORT.md`. Produced by a sweep agent; its deductions are unverified.

- No paper studies protection of individual core maxima or a floor-constrained subsidy in
  random multi-sided markets.
- The two headline questions as first posed are close to settled: Hassidim and Romm
  (2015) with P1's cut lemma give a protection ratio tending to 1/K on trees; Frieze,
  Pegden and Tkocz (2020) give a subsidy share vanishing at rate q^(-2/3).
- The tree/cycle divide probably does not survive in the ratio; it moves to core
  existence, which is the open problem (only iid-triple arrays at tiny sizes are known,
  Lucas 1995).
- Correction: the nonempty-core probabilities are from Lucas 1995, not Marchi-Oviedo.
- Local pilot on record: 190 of 190 nonzero instances with two agents per sector had a
  core; 77 of 78 with three (one empty case, float-flagged, not exactly certified).

## 2026-11-09: coordinator decisions recorded in the plan

- Question restated around core existence and a multi-sided law of one price.
- Gate 1 frozen with a three-way decision rule. Not run.

## 2026-11-09: research-handoff receipt enqueued

Enqueued by the coordinator and synced; recorded as pending review. The registry
does not change until a reviewed receipt supersedes it. Evidence bound: `PLAN.md` and
`GATE0_REPORT.md` as they stood at enqueue time.

- `typical-case-protection`: receipt `b7b2bfa270d54748a931bb04fc8f9ebf` (new lane).

## 2026-11-11: Lemma C push (four sub-agents) and synthesis

Max: "keep going on the open lemma - deploy subagents and solve it". Not solved. Details in
`lemmaC/SYNTHESIS_2026-11-11.md` and the four sub-folders; nothing reviewed by another model
family.

- Prior art: open (about 80% confidence); scoop risk low.
- Structure (q up to 48): the improving structure is one giant component touching 63-89% of
  M, not half-integral at large q, with no improving block of up to q/2 teams. The gap
  LP* - T grows (mean 0.039 at q=16 to 0.156 at q=48, sd flat near 0.054), which corrects
  the "O(1)" reading in THEOREM_VANISHING.md. P(mu* > 0) = 1.000 at q=48 (n=120).
- Concentration: Var T and Var LP* <= 4 E(3q - T) (Poincare for exponential coordinates),
  plus a conditional theorem through surrogates G_kappa.
- Dual: ordinal-limit Theorem O (the core is always nonempty when only the order of values
  matters); P(T = S) = 1/q!; dual depth = -mu*.
- Coordinator's synthesis, Proposition R: if E(LP* - T)/q^(1/4) -> infinity (q^(1/6) given
  E(3q - T) = O(q^(1/3))), then P(core nonempty) -> 0. The theorem now rests on a
  first-moment "two constants" conjecture, a_LP < a_IP.

## 2026-11-11: second push on the vanishing theorem (four sub-agents)

Max: "finish what's left to prove". Not finished. See the update section of
`lemmaC/SYNTHESIS_2026-11-11.md`. Proposition R was checked (FPT Theorem 2 covers this
model; one non-uniqueness wording fixed) and now reads: P(core nonempty) -> 0 if
E(LP* - T)/q^(1/6) -> infinity. What remains is the two-constants conjecture a_LP < a_IP.
Evidence for it: cavity a_LP = 2.197 and a_IP >= 2.2495 (non-rigorous); data to q = 64 with
E D = 0.212 and rising, offset fit c of about 0.11. Rigorous partial results in
`lemmaC/proof/`, all unreviewed.

## 2026-11-12: round 3 (multi-agent search); the vanishing theorem is claimed as proved

Max pasted the orchestration prompt saved in `lemmaC/PROMPT_TWO_CONSTANTS.md` (8 agents at
once, 40 runs). The orchestrating session ran it: 38 agent runs and four Codex reads. Records:
`lemmaC/round3/MAIN_THEOREM.md` (statement and proof map), `lemmaC/round3/CLAIMS.md` (every
claim with its status and who checked it), `lemmaC/round3/final/` (full write-outs, hardened
certificate, errata), `lemmaC/round3/codex/` (the Codex reads).

Claimed theorem: in the random triangle market P(core nonempty) -> 0. More precisely
E(LP* - T) >= (0.050 - o(1))·q^(1/3) and (LP* - T)/E(LP* - T) -> 1 in probability.

Proof in four parts, each written by one agent, attacked by two further fresh agents of the
same model family, and read at least once by Codex. No read found an error affecting a
theorem; fixes are in `lemmaC/round3/final/ERRATA.md`.

1. Variance: whp the optimal assignment and every LP optimum use only pairs and teams of
   deficit below C·q^(-1/2+eps), so Var(LP* - T) = O(q^(-1/6+eps)). The theorem then needs
   only liminf E(LP* - T) > 0.
2. Integer lower bound: whp 3q - T >= 2.302·q^(1/3). A first moment conditioned on agents'
   local stars in the model with independent team deficits (interval-arithmetic
   certificate, 2.30224586), carried to the pair-additive market by a comparison theorem.
3. LP upper bound: for any invariant matching m on the Poisson hypertree,
   limsup E(3q - LP*)/q^(1/3) <= c(m). Local convergence in total variation plus a patching
   argument with a split of each pair deficit.
4. A cheap invariant matching exists: the recursion X = min_k(xi_k - X' - X'') has a
   solution (the map is monotone in the dispersive order, which removes the shift
   instability) and every solution has cost 9·E X <= 2.2519 (computer-assisted certificate;
   numerically 2.2495456).

What this replaces: the two-constants conjecture of `lemmaC/SYNTHESIS_2026-11-11.md` is
proved in the form needed (integer constant at least 2.302, LP constant at most 2.2519);
Proposition R of that file is superseded by part 1. The finite-size fit "0.11·q^(1/3)" was
pre-asymptotic: near-tight teams at cost level K share a pair with about 0.9·K²·q^(-1/3)
others, so sizes below about 10^4 do not show the limit. The pair-additive LP constant
measured at q = 2,000 is 2.31, still falling towards about 2.2.

Other findings: counting price certificates fails (overshoots by exp(1.26·q)); two-round
exposure alone cannot give a positive gap; the LP optimum is fractional on about 90% of
agents; the cheap tree matching is not a function of the tree (it follows from parts 2 and
3). mpmath 1.3.0's interval exp and log are not strict enclosures; both certificates use
verified replacements.

Status and limits: no person has read any of it. Every reader was an AI model (one family,
plus Codex). The certificates trust the interval library's four basic operations and IEEE
arithmetic. All bounds are asymptotic (some thresholds exceed 10^50,000). The proof imports
Frieze-Pegden-Tkocz Theorem 2(b) (arXiv:1901.07167), read by several agents and by Codex as
an expectation bound for uniform pair costs. `draft/paper.tex` predates this round and does
not contain the theorem.

## 2026-11-12: consolidated manuscript assembled

Max: "Yes into one makes sense right?" Six agents wrote the sections in parallel from shared
notation (`manuscript/NOTATION.md`, `manuscript/macros.tex`); the coordinating session
assembled them. Result: `manuscript/main.tex`, 121 pages, compiles with no undefined
references or citations (`manuscript/build/main.pdf`). The earlier `draft/` is untouched and
predates the theorem.

- Sections: introduction, model, variance, independent-team bound, transfer, tree,
  existence, local convergence, LP bound, proof of the main theorem, certificates,
  discussion. The corrections in `lemmaC/round3/final/ERRATA.md` are applied.
- New text written during consolidation, each read once by Codex and found sound: a reworked
  step in the super-solution lemma (after a from-memory citation was removed); positive parts
  in the local-convergence density bound; a complete proof of the corollary that no invariant
  matching on the tree of cost below 2.302 is a function of the marked tree (the one-line
  argument in `MAIN_THEOREM.md` was incomplete).
- A rate follows from the pieces: P(core nonempty) = O(q^(-5/6+eps)).
- The paper states its results at beta_0 = 2.302, the constant covered by the hardened
  certificate; the row for 2.30224 rests on a referee's exact-arithmetic script.
- Two Codex reads of the full manuscript (`lemmaC/round3/codex/CODEX_READ_manuscript_*.md`):
  no fatal mathematical error; sixteen findings in the introduction, model, certificates and
  discussion text (overstatements, a wrong optimisation direction, an imprecise local-limit
  statement, rounding wording, one broken sentence). All were corrected by the coordinator.
- Bibliography: details taken from Crossref records; Prekopa 1973 has no record and is
  marked for checking. Eight marked items remain (published wording of the
  Frieze-Pegden-Tkocz theorem, three changes owed to cert.py, one unsourced reference,
  textbook theorem numbers).
- Known inconsistency left in place: the independent-team section names the sectors
  V_1, V_2, V_3 (B is its number of bins); the other sections use A, B, C.
- Status unchanged: written and checked by AI models only; no person has read it.

## Next entry

Gate 1 result and the decision rule's outcome, recorded before any further analysis.

## 2026-11-11: gate 1 run (Claude Code lane agent), numbers and decision rule

Run as frozen in `PLAN.md`. Floating-point evidence only (scipy 1.15.3 HiGHS LPs and
`milp`); nothing certified exactly; not yet reproduced by a different model.
Code `gate1/gate1.py`, `gate1/summarize.py`; drivers `gate1/run_controls.sh` (run first),
`gate1/run_main.sh`; per-instance records `gate1/out/*.jsonl`, aggregates
`gate1/out/summary.json`, `gate1/out/summary_table.txt`. Seeds:
`numpy SeedSequence([20261111, graph, dist, q, idx])` (graph 0 triangle, 1 path; dist 0
U[0,1], 1 Exp(1)); idx 0..n-1. Wall time about 21 minutes (controls 8, main 13) on 6
workers while the machine carried load average 100-175 from other work; 6,796
per-instance CPU-seconds.

**Controls (run first): all pass.**
- Path (C = 0), U[0,1], 500 instances at each q in the grid: core nonempty 5000/5000;
  rho = 1/3 in every instance (min 0.33333, never below); s = 0 in every instance
  (|s| < 2e-13).
- Triangle q = 2, U[0,1], 2000 instances (the main-run seeds): core nonempty 2000/2000.
- Internal checks on every group: LP duality (core min sum p = LP value, rel err <= 3e-14);
  MILP T equals integral-vertex T where both exist (diff <= 1.1e-14); no MILP left
  unsolved; no fractional-vertex-with-zero-gap instance; the P7 link s = 0 iff
  rho >= 1/3 holds in every nonempty instance; rho <= T / sum H everywhere.

**Triangle, U[0,1].** ne = nonempty-core count; CI = Clopper-Pearson 95%; gap =
(LP value - T)/T; deficit = 1/3 - rho; spread = largest within-sector payoff range
(absolute units, pair values in [0,1]), at the min-spread rho-optimal point (and at
the HiGHS rho-LP vertex).

| q | n | ne | share | CI | gap mean (all) | gap median (empty) | rho mean | rho min | P(rho<1/3 \| ne) | P(s=0 \| ne) | s/T max | spread median (vertex) | T/sumH |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 2000 | 2000 | 1.000 | [0.998,1.000] | 0 | - | 0.33333 | 0.33333 | 0.000 | 1.000 | 0 | 0.303 (0.366) | 0.371 |
| 3 | 2000 | 1987 | 0.994 | [0.989,0.997] | 0.00004 | 0.00480 | 0.33317 | 0.3064 | 0.019 | 0.981 | 0.052 | 0.441 (0.524) | 0.384 |
| 4 | 2000 | 1939 | 0.970 | [0.961,0.977] | 0.00016 | 0.00384 | 0.33295 | 0.3056 | 0.058 | 0.942 | 0.038 | 0.486 (0.561) | 0.388 |
| 5 | 2000 | 1873 | 0.936 | [0.925,0.947] | 0.00025 | 0.00307 | 0.33289 | 0.3134 | 0.087 | 0.913 | 0.020 | 0.516 (0.591) | 0.391 |
| 6 | 2000 | 1784 | 0.892 | [0.878,0.905] | 0.00033 | 0.00225 | 0.33276 | 0.3163 | 0.121 | 0.879 | 0.014 | 0.532 (0.593) | 0.392 |
| 8 | 2000 | 1568 | 0.784 | [0.765,0.802] | 0.00055 | 0.00201 | 0.33258 | 0.3167 | 0.183 | 0.817 | 0.012 | 0.525 (0.584) | 0.392 |
| 10 | 2000 | 1338 | 0.669 | [0.648,0.690] | 0.00075 | 0.00180 | 0.33249 | 0.3147 | 0.241 | 0.759 | 0.013 | 0.512 (0.563) | 0.390 |
| 12 | 1000 | 582 | 0.582 | [0.551,0.613] | 0.00077 | 0.00147 | 0.33217 | 0.3190 | 0.357 | 0.643 | 0.007 | 0.492 (0.538) | 0.386 |
| 16 | 1000 | 366 | 0.366 | [0.336,0.397] | 0.00103 | 0.00146 | 0.33234 | 0.3205 | 0.391 | 0.609 | 0.003 | 0.463 (0.495) | 0.383 |
| 20 | 1000 | 211 | 0.211 | [0.186,0.238] | 0.00112 | 0.00131 | 0.33213 | 0.3229 | 0.559 | 0.441 | 0.004 | 0.421 (0.444) | 0.376 |

Largest deficit 1/3 - rho by q: 0, 0.0269, 0.0277, 0.0199, 0.0170, 0.0166, 0.0186,
0.0143, 0.0129, 0.0104. Mean deficit over nonempty instances never exceeds 0.0012.
rho never exceeded 1/3 (max 1/3 + 2e-9) in any instance of any group.

**Contrast, triangle, Exp(1), 500 per q.** Share nonempty: 1.000, 0.996, 0.994, 0.984,
0.972, 0.932, 0.908, 0.824, 0.782, 0.658 (q = 2..20). rho mean 0.3323-0.3333, min
0.3058-0.3333, P(rho < 1/3 | ne) 0 to 0.17. Median min-spread grows: 1.06, 1.63, 2.18,
2.41, 2.69, 3.05, 3.37, 3.55, 3.81, 4.10 (one price fails, as expected).

**Decision rule, applied mechanically.**
- Vacuous (share < 0.05 by q = 10-12): **not triggered.** Share 0.669 at q = 10 and
  0.582 at q = 12 (lower CI 0.648, 0.551).
- Original proposal revived (rho bounded below 1/3): **not triggered.** rho approaches
  1/3: the largest deficit falls from 0.028 (q = 4) to 0.010 (q = 20) and the mean
  deficit stays at or below 0.0012.
- Paper shape: **the outcome, by elimination, provisional.** Of its three clauses, only
  "ratio tends to 1/3" is supported. "Share bounded away from zero" is not shown: the
  share falls with no sign of a plateau (0.784, 0.669, 0.582, 0.366, 0.211 at
  q = 8..20); a straight-line extrapolation, not a fit, would cross 0.05 near q = 24.
  "Spread tends to zero" is weak: median min-spread falls only from 0.532 (q = 6) to
  0.421 (q = 20), and the tree control behaves alike (0.396 at q = 6, 0.262 at q = 20).
  The rule does not cover this case (vacuity threshold missed at q = 10-12 but the share
  heading toward it soon after), so the label should be read as "not vacuous yet".

**Frozen-definition ambiguities and how they were resolved (fixed before the main run).**
1. The rho-optimal point is not unique. Reported the spread at the vertex HiGHS returns
   for the rho LP and the minimum over the rho-optimal face (t >= rho - 1e-9) of the
   largest within-sector range; the rule is read on the minimum, the reading most
   favourable to "spread tends to zero". Both give the same picture.
2. "Spread" = max over the three sectors of (max - min payoff in the sector), absolute.
3. Core nonempty iff (LP value - T)/T <= 1e-7, T exact by `milp` for every instance at
   q <= 6 and for every fractional-vertex instance at q > 6; integral dual-simplex
   vertex otherwise. Every MILP solved, so the gap is available for all instances.
4. Counts above q = 10 set to 1000 (frozen: at least 500). Contrast size not frozen: set
   to 500 per q on the same grid. Path control size not frozen: 500 per q.
5. The q = 2 triangle control and the q = 2 main row are the same 2000 instances.
6. Core includes p >= 0 (P1, P7). H_i > 0 threshold 1e-9; smallest H_i seen 0.034.

**Observations recorded for later, not analysed.** P(s = 0 | core nonempty) falls from
1.000 to 0.441 over q = 2..20, against the plan's question 3 guess that it tends to one;
on Exp(1) values rho stays near 1/3 although the spread grows, so rho near 1/3 is not
being driven by small dispersion.

Not yet recorded in the research handoff; the coordinating session files receipts.

## 2026-11-11: POST-GATE EXTENSION (not part of the frozen gate 1): longer core-existence curve, exact empty-core certificates, two proofs

Approved by Max via the coordinator after gate 1. Claude Code lane agent. Nothing below changes
the gate 1 entry or its decision-rule outcome. Floating-point results are evidence, not
certificates; nothing here reproduced by a different model yet. Files under `gate1_ext/`
(`ext.py`, `run_ext.sh`, `fit.py`, `cert.py`, `rho_probe.py`; outputs in `gate1_ext/out/`) and
`PROOF_NOTES.md`.

**1. Core-existence share, triangle, U[0,1], same seeds as gate 1** (SeedSequence([20261111, 0,
0, q, idx]), idx 0..n-1). Share only; no ratio or subsidy work. Core test as gate 1 (dual-simplex
LP vertex integral => nonempty). q = 24: fractional instances also solved by MILP (gap known,
0 fractional-vertex/zero-gap cases). q >= 28: **LP-only mode** (fractional vertex => empty, no
MILP, no gap), because under machine load (load average 60-120) not one q = 40 MILP finished in
5 minutes; this relies on the a.s.-unique LP optimum, and gate 1 saw 0 fractional-vertex/zero-gap
cases in 21,000 triangle instances. CP = Clopper-Pearson 95%.

| q | n | nonempty | share | CP 95% | mode |
|---|---|---|---|---|---|
| 24 | 500 | 61 | 0.122 | [0.095, 0.154] | LP + MILP (median empty gap 0.00117) |
| 28 | 500 | 40 | 0.080 | [0.058, 0.107] | LP only |
| 32 | 500 | 19 | 0.038 | [0.023, 0.059] | LP only |
| 40 | 2000 | 16 | 0.008 | [0.005, 0.013] | LP only (first 500: 6) |
| 48 | 1000 | 2 | 0.002 | [0.000, 0.007] | LP only (q = 48 not requested; added for the fit) |

**Fits** (binomial maximum likelihood over all 15 points q = 2..48, gate 1 rows included;
`gate1_ext/out/ext_table_and_fits.txt`, `ext_summary.json`). Deviance against the saturated model:
- exp(-c (q-2)^a): c = 0.0115, a = 1.68; deviance 20.3 on 13 dof. Fits. Predicts 0.006 at
  q = 40, 0.0008 at q = 48, 0.0005 at q = 50.
- same plus a floor f: f = 0.003, deviance 15.5 on 12 dof (likelihood-ratio 4.8 on 1 dof, a floor
  at the 0.3% level is not excluded but not needed; it is driven by q = 40-48 having 18 nonempty
  of 3000 against about 8 predicted).
- logistic in log q (power tail, ~q^-3.1): deviance 111 on 13. Rejected.
- exp(-c q^a) unshifted: a = 2.08, deviance 190. Logistic in q: deviance 588. exp(-c (q-2)): 1505.
  Logistic plus floor: 494. All rejected.
Reading: the share falls faster than exponentially in q (a stretched/compressed exponential with
exponent about 1.7), passes 0.05 near q = 30 and is below 1% by q = 40. On this evidence it goes
to 0; a floor above about 1% is excluded by q = 40-48, and nothing here rules out a tiny residual
floor. Not a proof.

**2. Exact empty-core certificates, q = 3** (`gate1_ext/cert.py`, output
`gate1_ext/out/certificates_q3.txt`). Gate 1 seeds with pair values rounded to hundredths (so each
certified instance is its own rational instance, not the float one). For each: a fractional
primal lam (all six masses 1/2), a dual y >= 0 with y(e) >= v(e) for all 27 teams and
sum y = LP* (so LP* is exact), the integer optimum T by enumeration of all 36 perfect matchings
(values >= 0), all checked in Fractions; the core is empty because any core point would be dual
feasible with total T < LP*.
- idx 54: LP* = 113/20, T = 559/100, gap 3/50 (relative 0.0107).
- idx 1507: LP* = 291/50, T = 579/100, gap 3/100 (0.0052).
- idx 2225 (beyond gate 1's idx range; same seed scheme): LP* = 309/50, T = 609/100, gap 9/100
  (0.0148).
Twelve other float-empty candidates in idx 0..2225 failed basis recovery after rounding (the
rounding can change the instance) and were skipped, not certified.

**3. Proofs** (`PROOF_NOTES.md`; complete as written, unreviewed).
- Theorem 1: every nonnegative pair-additive triangle market with q = 2 has a nonempty core.
  The LP objective factors through three 2 x 2 doubly stochastic marginals (s, t, u); the image of
  the fractional polytope is exactly the tetrahedron with vertices (1,1,1), (1,0,0), (0,1,0),
  (0,0,1), whose vertices are the four perfect matchings, so LP* = T.
- "rho <= 1/3 on every triangle market" is FALSE as stated: exact q = 2 path-market
  counterexample with a one-point core and rho = 1 (players with H_i = 0 drop out). Theorem 2
  (corrected): with equal sector sizes, any team values and any sector graph, if every H_i > 0
  then rho <= 1/K (= 1/3 on the triangle). Proof by the sector-shift lemma: shifting payoff
  between sectors keeps a core point in the core, so H at each sector's minimiser is at least the
  sum of the sector minima. With P1's tree bound this gives rho = 1/K exactly on trees when all
  H_i > 0 (matches gate 1's path control).
- Side observation (float, `gate1_ext/rho_probe.py`): 1,795 nonempty small markets under
  Bernoulli, small-integer, Pareto and U^4 laws; the largest rho was 1, at the degenerate
  instance above (how many others exceed 1/3 was not tabulated); 90 small perturbations of it
  all give rho = 1/3.

Not recorded in the research handoff; the coordinating session files receipts.

## 2026-11-11: first draft working paper

First draft (13 pp., author TBD) in `draft/` (`paper.tex`, `paper.pdf`, figures from stored outputs via `draft/make_figures.py`); all theorems marked "proof not yet independently checked", simulations described as floating-point evidence; not reviewed, not circulated.

## 2026-11-11: vanishing-core attempt (Claude Code lane agent, 60-minute box): no limit theorem; strategy redirected

Goal: P(core nonempty) -> 0 on the U[0,1] triangle market. Not proved. Write-up
`THEOREM_VANISHING.md`; scripts and outputs `vanishing/` (floating point unless stated; not
reproduced by a different model).
- Complete (unreviewed): Lemma 1 (localization to sub-collections of the integer optimum M);
  Lemma 2 (signed stability margin mu* = max_w v.(w-M)/d(w); core empty iff mu* > 0);
  Proposition 3 (P(empty) >= (1/(200q))^{3q^2} > 0 for every q >= 3, by embedding the exact q = 3
  certificate idx 54, re-checked in rationals by `vanishing/verify_cert54.py`). No uniformity in q.
- Global route refuted (`global_check.py`): S - LP* grows 2.84 -> 4.50 (q = 24 -> 48), 3q - LP* ~
  2.55 q^{1/3}, LP* - T ~ 0.07 at q = 24.
- Bounded local route contradicted (`local_count.py`, exact 3x3x3 vertex list from cdd `verts3.py`,
  float): empty-core sub-collections of 3 M-teams found in <= 1 empty instance per q (0 at q = 20,
  24); of 4 M-teams 10/22 (q = 6) down to 4/62 (q = 12). The LP-optimal and the most efficient
  improving structures disrupt ~0.75-0.8 q M-teams (median 19-20 at q = 24); weights mostly 1/2.
- Reduction: Conjecture follows by Chebyshev from Lemma C (E mu* >= m_q > 0, Var mu* = o(m_q^2)),
  open. Evidence (`margin.py`, 120 per q): mean/sd of mu* = -1.15, -0.73, -0.42, 0.43, 1.24 at
  q = 4, 6, 10, 16, 24; sd falls about q^{-1.4}.
- Gaps: Lemma C; cdd vertex counts not certified; heuristics (triangle-factor threshold, q^{k-2r/3}
  count) unproved; no uniform lower bound on P(empty).
