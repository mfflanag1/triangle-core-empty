# Lemma C: synthesis of the four-agent push, 11 November 2026

Written by the coordinating session from `prior_art/PRIOR_ART.md`, `structure/STRUCTURE.md`,
`concentration/CONCENTRATION.md` and `dual/DUAL.md`. Nothing here has been reviewed by
another model family. Setting as in `THEOREM_VANISHING.md`.

## The reduction (coordinator's combination of A2 with the structure data)

Write each pair value as X_i = 1 - exp(-E_i), with E_i iid Exp(1). T and LP* are maxima of
linear functions of X, so almost everywhere

- dT/dE_i = 1{i is an M-pair} (1 - X_i), giving |grad T|^2 <= sum over M-pairs of (1 - X_i)
  = 3q - T;
- dLP*/dE_i = lambda_i (1 - X_i), where lambda_i in [0,1] is the pair's weight in the LP
  optimum, giving |grad LP*|^2 <= sum_i lambda_i (1 - X_i) = 3q - LP* <= 3q - T.

The Poincare inequality for product Exp(1) measure, with constant 4, then gives

  Var T <= 4 E(3q - T)  and  Var LP* <= 4 E(3q - T).

This is the concentration agent's A2 in its sharpest form. With E(3q - T) = O(sqrt q),
proved there by an explicit construction, both variances are O(sqrt q). If
E(3q - T) = O(q^(1/3)), as the gate-0 sweep reports from Frieze-Pegden-Tkocz (this still
needs checking against the paper) and as the data show (3q - T is about 2.55 q^(1/3)),
both variances are O(q^(1/3)).

Chebyshev on the gap D = LP* - T, with Var D <= 2 Var T + 2 Var LP*, gives:

**Proposition R (complete given A2; unreviewed).** If E(LP* - T) / q^(1/4) -> infinity,
then P(core nonempty) -> 0. With the q^(1/3) bound on E(3q - T), q^(1/6) suffices.

So the theorem now reduces to a **first-moment statement** about the gap, with no further
concentration needed. Since D <= S - T = O(q^(1/3)), where S is the sum of the three
bipartite optima, the natural target is the structure agent's Lemma G(ii):

**Conjecture (two constants).** E(LP* - T) >= c q^(1/3). Equivalently, the per-agent
deficit constants of the fractional and the integral optimum differ: a_LP < a_IP.

## Evidence

- E D = 0.039, 0.069, 0.108, 0.127, 0.156 at q = 16, 24, 32, 40, 48 (200/200/200/200/120
  instances). D / q^(1/3) is still rising. The empirical sd of D is flat at about 0.054,
  far below the rigorous O(q^(1/6)) bound.
- P(mu* > 0) = 1.000 at q = 48 (120 instances).
- The improving structure is a single component that touches 63-89% of M, has
  non-half-integral vertices at large q, and does not decompose: no block of up to q/2
  M-teams has a positive sub-gap. Local and second-moment-over-pieces routes look dead.
- Near-tight off-M teams per agent at slack c q^(-2/3) are bounded in q, so a local limit (a
  weighted three-partite hypertree, like Aldous's PWIT) is plausible.

## What is known not to work

- Ordinal arguments (Theorem O in `dual/`): if only the order of pair values matters, the
  core is always nonempty. Any proof must use how close the values are, at scale q^(-2/3).
- Split pair-price certificates: P(T = S) = 1/q! exactly.
- Bounded local structures: their expected count appears to go to 0.

## Prior art

Open, with confidence about 80% (`prior_art/PRIOR_ART.md`). Nearest analogue: random
matching on K_n, where the LP and integer costs share the limit pi^2/12 and differ at
order 1/n. That points to a different regime (O(1) odd cycles), so it probably doesn't
transfer. Tools: Hessler-Wastlund concentration; the replica and cavity cycle expansion
(Lucibello-Malatesta-Parisi-Sicuro 2018; Parisi-Perrupato-Sicuro 2020); Haber-Krivelevich
and Johansson-Kahn-Vu for first-order factor thresholds.

## Next steps

1. Check the Frieze-Pegden-Tkocz statement for E(3q - T) in this model, and have a
   different model family check A2 and Proposition R.
2. Run the cavity / population-dynamics computation of a_LP and a_IP in the local limit,
   at the q^(-2/3) scale. If the fixed points differ, the two-constants conjecture is
   supported and the proof route is an objective-method argument (Aldous; Wastlund).
3. Extend the E D series beyond q = 48, which is cheap with the structure agent's
   reduced-cost MILP.

## Update after the second push (same day, four more sub-agents)

Files: `check/CHECK.md`, `proof/PROOF_ATTEMPT.md`, `cavity/CAVITY.md`, `scaling/SCALING.md`.
Nothing has been reviewed by another model family.

**Proposition R, checked (same model family).**
- Frieze-Pegden-Tkocz (arXiv:1901.07167; DAM 287, 2020), Theorem 2: under the map U -> 1 - U
  their model is ours exactly, so E(3q - T) = Theta(q^(1/3)). The q^(1/6) threshold
  therefore holds, and the earlier hedge is withdrawn.
- One fix: the LP optimum is not a.s. unique in team weights. Two vertices can carry
  identical pair loads, as in an exact q = 2 example. The optimal pair-load vector is a.s.
  unique, which is all the argument needs.
- Poincare applied to D directly gives the sharper bound
  P(core nonempty) <= 8 E(3q - T) / (E D)^2.

**Statement.** P(core nonempty) -> 0 if E(LP* - T)/q^(1/6) -> infinity.

**What remains: one conjecture, a_LP < a_IP**, i.e. E(LP* - T) >= c q^(1/3).
- Rigorous partial results (`proof/`, labelled complete, unreviewed):
  - 3q - T >= 2.15 q^(1/3) whp, by a union bound;
  - E(3q - LP*) >= 1.36 q^(1/3);
  - conditioning on M lowers LP* (Harris);
  - D > 0 forces an odd Berge cycle in the LP support.
- Routes that are dead: separating the constants with current bounds (needs 1.7% precision,
  best is 18% off); bounded constructions; local blocks.
- Cavity, non-rigorous: a_LP = 2.197 +- 0.002 and a_IP >= a_RS = 2.2495, so
  a_IP - a_LP >= 0.052. The integral fixed point is unstable to fractional perturbations,
  and a bipartite control behaves correctly.
- Data to q = 64: E D = 0.212 at q = 64, sd about 0.04, P(D > 0) = 1 from q = 48.
  E D/q^(1/3) is still rising (0.053). The offset fit c q^(1/3) - d gives c of about 0.11,
  consistent with the cavity bound and with replica-symmetry breaking raising a_IP. The
  finite-q exponent (about 1.1) is pre-asymptotic, so it cannot certify the hypothesis.

**Proof program for the conjecture.**
1. Local convergence of (near-tight hypergraph, M) at scale q^(-2/3) to the Poisson
   3-uniform hypertree (intensity xi^2/2).
2. On that tree, show the integral cavity solution is not the cheapest invariant fractional
   perfect matching. This is a single-object statement and the cleanest next target.
3. Transfer the limit constants back to E(3q - T) and E(3q - LP*), objective-method style.

A useful first lemma: whp M uses no pair with deficit above q^(-2/3+eps). That gives
Var T, Var LP* = o(1), so Proposition R would need only E D -> infinity, and it is also
the main input to step 1.
