# Corrections owed to the round-3 files

Collected by the orchestrating session on 12 November 2026 from the referees' and Codex's
reads. None changes a theorem. Items marked "applied" are already in the file named; the
rest are owed in a consolidated write-up. Sources of each item are in `../CLAIMS.md`.

## final/THEOREM_1.md (Theorem 1)
- Applied: Theorem A, Step 6. The display for p·|G| had dropped a factor p; corrected, with
  kappa = 0.3^{r_A}/(8·4^{r_A+1}). Only one exponentially small term changes.
- Applied: "(Q3) the three exponentially small terms" now reads "the exponentially small
  terms".
- Owed: state the definition of the core (so that "core nonempty iff D = 0" is
  self-contained); note that for eps >= 1/2 the hypothesis of part (iii) is vacuous.

## firstmoment/THEOREM_I.md (Theorem I)
- Owed: section 5 must cite `final/certI/certify_rigorous.py`, not `certify.py`, whose
  library exponential and logarithm are not strict enclosures. The certified value is
  Psi_0 <= -3.84751328322 and -Psi_0/theta >= 2.30224586.
- Owed: in (2.2) keep the factor prod_b (n_b + 1)^2 inside the maximum over degree arrays
  (it is at most (1 + (lambda_max + delta)·q)^(2B) on typical arrays); state q^(2/3) >= K.
- Owed: define "simple" so that it excludes repeated triples across bins as well.
- Owed: two displayed numbers in §5.3 overstate by one unit in the last digit (unused).

## interface/THEOREM_I_PA.md (Theorem I-pa)
- Owed: the prefactor constant is C = C_0(K, B)·(1 + lambda_max + delta_1 + 2·tau)^(2B),
  not a function of K and B alone (`final/certI/NOTES.md` §4).
- Owed: write the condition on A with the true ell_max = 8.2718 (so the coefficient in
  Lemma Q is <= 0, not = 0).
- Owed: add a table row for beta = 2.30224 (one that works: c' = 21, delta_1 = 5e-9,
  delta_2 = 5e-7, tau = 5e-10, exponent <= -7.5e-6).
- Owed: the "size of q" remark is far too small. The error term needs (tau/3)·ln q > 0.203,
  i.e. q above about 10^52,900 at beta = 2.30. The theorem is asymptotic only.

## sc/TRANSFER.md (Theorem T)
- Owed: call the comparison model "independent Bernoulli teams", not "Poissonised".
- Owed: Lemma L2 is the Boppana-Spencer proof of Janson's inequality; drop the phrase
  "using Harris only".

## localconv/LOCALCONV.md
- Owed: the remark that the first term dominates whenever the bound is below 1 is false for
  small K; the fluctuation rate in §2 is O_P(q^(-1/6)), only the bias is O(q^(-1/3)).
- Owed: the "repair corollary" of §4 is not justified for an arbitrary rule and is not used
  anywhere; replace it by Lemma 6 of LP_LOCALITY.md (`final/TREE_AND_LOCALITY.md` (e)).

## lplocal/LP_LOCALITY.md (Theorem 3)
- Applied in `final/TREE_AND_LOCALITY.md` (d): section 3 is replaced by the argument from
  total-variation convergence of marked balls; the 1/(3r) error terms disappear.

## reports/tree.md and rde/RDE.md (Theorem 4)
- Applied in `final/TREE_AND_LOCALITY.md` (a)-(c): the convention at x_f = 1; the hypotheses
  of Lemma D; the full construction of the message process on the infinite-degree tree.
- Owed: Theorem K(i) of reports/tree.md needs the no-tie hypothesis g^{k*}_{n-1} < g^{k*}_n
  for its "only if" (not used in the main line).

## rdeB/ROUTES.md (existence)
- Owed: add the bound E X_-^k <= E Y_+^(k+3)/((k+1)(k+2)(k+3)) so that a fixed point has all
  absolute moments.
- Owed: in Lemma 4 replace the from-memory citation "a sum of IFR variables is IFR" by the
  direct argument: E(t + Y)_+^2/2 is the integral of P(X_1 > -u_1)·P(X_2 > -u_2) over
  u_1 + u_2 <= t, which is log-concave by Prekopa.
- Owed: "all exponential moments" should read "exponential moments for every positive
  parameter"; the super-solution's multiplier is existential, not numerical.
- Owed: the printed c_n are low by about 1e-4; use 3.41528213, 2.66021540, 2.39808802,
  2.30169894, 2.26786005, 2.25596056, 2.25179199, 2.25033206, 2.24982093, 2.24964199.

## cert/CERTIFICATE.md and cert.py (cost certificate)
- Owed: quote the bound as 9·E X <= 2.2519. The printed 2.251887 is a round-to-nearest
  format of an upward-rounded internal bound; assert the internal bound against an exact
  rational instead.
- Owed: make cert.py refuse lattice parameters other than the two reported (h = 2^-10 and
  2^-12), for which the exact-argument justification is written; save the stage files of
  the fine run (regenerated copies are in `../v10/`).

## One numerical loose end (not on the proof path)
- One referee's test of the cost-comparison lemma on NON-IFR pairs with slopes up to 400
  returned +4e-5; a bounded-slope rerun was negative. The lemma is applied only to IFR laws.
  Worth re-evaluating that pair.

## final/TREE_AND_LOCALITY.md (after the fourth Codex read)
- Applied: §D.3 (the optional per-shape alternative, not used on the main line): a false
  inequality repaired by inserting a positive part; the zero-dimensional shape (isolated
  root) handled separately.
- Applied: §E (the corrected side corollary, not used on the main line): the extra
  randomness Z must make the pair (u, Z) independent of w; "Z independent of w" was not
  enough.
- Owed: ship the weights and export files with `final/certI/` so that the hand-off can be
  reproduced from the folder alone (they are in `firstmoment/` and `final/certI/`; the Codex
  packet did not include them).
