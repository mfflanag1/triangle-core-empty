# Typical-case protection in random multi-sided markets: plan

Frozen 2026-11-09 by the coordinating session, before any lane computation.
Proposed lane ID: `typical-case-protection` (not yet registered in the research handoff).
Status: prior-art sweep passed, narrowed. The lane does not open as first proposed.
Gate 1 not run.

Read `GATE0_REPORT.md` in this folder first. It is the evidence for everything below.
Its deductions are the sweep agent's own and are unverified.

## Why the original question changed

The first proposal asked whether P1's tree/cycle divide in the protection ratio and P7's
subsidy survive for random markets. The sweep found that both headline answers are close
to settled by prior art plus short arguments:

- **Trees.** Hassidim and Romm (arXiv 1404.6103, EC 2015) prove an approximate law of one
  price in random bilateral assignment games. With P1's cut lemma this implies the
  protection ratio tends to 1/K: the worst-case tree constant is also the typical value.
- **Cycles.** With equal sector sizes, the same shift argument gives a ratio near 1/K on
  any graph whenever the core is nonempty and within-sector payoffs are nearly equal. The
  divide moves from the ratio to core existence, which is certain on trees and unknown on
  cycles.
- **Subsidy.** For bounded iid values, s(M) <= 3q - T, and Frieze, Pegden and Tkocz
  (DAM 2020) give E[3q - T] of order q^(1/3). So s/T vanishes at rate q^(-2/3) for free.
- Both P1's and P7's quantities are undefined when the core is empty, which is the
  generic risk for random cyclic markets.

A correction to carry into any write-up: the published probabilities of nonempty cores
for three-sided games are in Lucas (Duke Math. J. 1995, reporting Crowley and Graham),
for iid triples, not in Marchi and Oviedo.

## Question now

For a connected sector graph G on K sectors with q agents per sector and iid pair values
from a bounded law with continuous density (benchmark U[0,1]):

1. **Core existence.** What is the probability that the core is nonempty, and how does it
   behave as q grows on cyclic G? Known: 1 on trees; 1 on the triangle at q = 2 (sweep
   agent's argument); apparently below 1 from q = 3. This is the real open problem.
2. **Multi-sided law of one price.** How fast does within-sector payoff dispersion vanish,
   and does the rate depend on the sector's degree in G? The protection ratio tending to
   1/K on every G would follow as a corollary.
3. **Second order on the triangle.** Given a nonempty core, does the probability that the
   P7 subsidy is zero tend to one?
4. **Optional contrast.** Unbounded (exponential) values, where the one-price law fails.

What this would give P1 and P7 even in the smallest outcome: a robustness remark that the
sharp tree constant is typical, not knife-edge, and that the cyclic worst case is driven
by instances that are not typical.

## Gate 1 (frozen; run exactly this before anything else)

**Setting.** Triangle, pair values iid U[0,1], q in {2, 3, 4, 5, 6, 8, 10, 12, 16, 20},
2,000 instances per q for q <= 10 and at least 500 above. Fixed, recorded seeds.

**Compute.**

1. Solve the primal LP (maximise total team value over fractional matchings) and record
   whether the optimum is integral. With continuous values the optimum is almost surely
   unique, so the core is nonempty exactly when it is integral. Confirm against an exact
   integer solve for q <= 6. Report the nonempty-core share by q with a binomial interval,
   and the relative gap (LP value minus T) / T where T is available.
2. On nonempty-core instances: every H_i, the protection ratio (one extra LP: maximise t
   subject to p in the core and p_i >= t * H_i), the P7 subsidy s, and the largest
   within-sector payoff spread at the ratio-optimal point.
3. Controls that must pass before any result is read: the path (third matrix zero) gives
   nonempty core always, ratio at least 1/3 and s = 0; q = 2 on the triangle gives a
   nonempty core always. Contrast run: exponential values.

**Decision rule, applied mechanically.**

- *Vacuous:* the nonempty-core share is below 0.05 by q = 10 to 12. P1's and P7's
  quantities are then undefined for typical cyclic markets. Stop; report; the lane
  survives only if Max wants it restated on the LP-dual optimal face.
- *Paper shape:* the share stays bounded away from zero, the ratio tends to 1/3 and the
  spread tends to zero. The lane is "a multi-sided law of one price and the probability
  of a nonempty core", with protection as a corollary.
- *Original proposal revived:* the ratio stays bounded below 1/3. The sweep agent's
  deduction is then wrong and the lane opens as first proposed.

Record the numbers and the rule's outcome in `LEDGER.md` before any further analysis.

## After gate 1, if it does not stop the lane

- Check the two unverified deductions as written proofs: the tree limit via Hassidim and
  Romm (confirm their balanced-case proof uses no sign constraints), and nonemptiness at
  q = 2 for every nonnegative pair-additive triangle market.
- Then the dispersion rate, and core existence on other cyclic graphs (C4, K4).

## Code to reuse (all read-only; copy into this folder)

- `structural_classification_2026-09-06/sharp_subsidy/discovery/probe.py`, function `lp`:
  row construction and the H and s LPs. Its total is brute force over permutations and
  unusable past q of about 4 to 5; replace it with the primal LP. Its import of
  `certificates` points at a `/private/tmp` path that no longer exists; the module is at
  `structural_classification_2026-09-06/simultaneous_protection/certificates.py`.
- `pair_additive_campaign_2026-09-05/experiment/pilot.py`, function `float_check`.
  **Naming hazard:** `rho` there is the Lmin-Lmax distance, not P1's protection ratio.
- No general solver for P1's protection ratio exists in the repository. Write it.
- Interpreter with scipy, numpy, z3 and ortools: `/Users/maxf/miniforge3/bin/python3`.

## Ground rules

- Write only in this folder. Everything else in the repository is read-only, and the P1
  and P7 submission packages under `/Users/maxf/Documents/ECON/artifacts/` must not be
  touched: both are awaiting submission approval.
- The repository has 23 uncommitted entries from other work. Do not commit, and do not
  run `git add -A`.
- Floating-point LP results are evidence, not certificates. Say so, and certify exactly
  anything that goes into a theorem.
- No outreach, posting or submission.
- A result is not verified until a different model has independently reproduced it.
  Record milestones with the research-handoff skill
  (`/Users/maxf/.codex/skills/research-handoff/SKILL.md`), as pending until reviewed.

## Sweep limits

Conference proceedings were covered by title keywords through OpenAlex and Crossref, not
DBLP. Google Scholar, SSRN full text and the Crowley and Graham dissertations were not
searched or read. Mordant (arXiv 2608.05123) and Yokote (arXiv 2604.02559) were read as
abstracts only and may supply tools.
