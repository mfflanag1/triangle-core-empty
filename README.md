# Asymptotic core emptiness in a random three-sided assignment market

## The result

Three sectors of *q* agents. Each of the 3*q*² pairs of agents from different sectors has an
independent U[0,1] value; a team (one agent per sector) is worth the sum of its three pair
values. *T* is the best total value of a perfect assignment into teams, LP\* the optimum of
the fractional relaxation, and *D* = LP\* − *T*. *D* is the minimum outside subsidy that makes the grand coalition stable. The core is nonempty exactly when *D* = 0.

**Theorem.** As *q* → ∞,

- E *D* ≥ (0.0501 − o(1))·*q*^(1/3), so the expected aggregate subsidy has order *q*^(1/3) and the expected subsidy per agent has order *q*^(−2/3).
- *D* / E *D* → 1 in probability.
- exp(−*Cq*) ≤ P(core nonempty) ≤ exp(−*cq*^(1/3)) for some positive constants *c*, *C* and all sufficiently large *q*.

The normalized fractional deficit (3*q* − LP\*) / *q*^(1/3) converges in mean square to the **minimum expected root cost of an invariant perfect matching** on the limiting Poisson hypertree. The minimum is attained and lies between 6^(1/3)Γ(4/3) and 2.2519. Its numerical value is open.

With high probability, every optimal fractional assignment splits a positive fraction of each sector’s agents across teams, with less than 99% of their weight on any one team. Corollary 10.3 also covers near-optimal assignments of total deficit at most 2.26·*q*^(1/3).

## Proof and supporting materials

| Part | Statement | Source file | Manuscript |
|---|---|---|---|
| 1. Variance | Whp no used pair or LP-support team has deficit above C·q^(−1/2+ε); Var *D* = O(q^(−1/6+ε)) | `lemmaC/round3/final/THEOREM_1.md` | §3, Appendix A |
| 2. Integer bound | Whp 3q − T ≥ 2.302·q^(1/3) (independent-team model, then transfer) | `lemmaC/round3/firstmoment/THEOREM_I.md`, `sc/TRANSFER.md`, `interface/THEOREM_I_PA.md` | §4–5 |
| 3. Fractional limit | Tree-to-market construction and marked local limits identify the fractional constant | `manuscript/sections/09_lp.tex` | §8–9 |
| 4. Existence and cost | The tree recursion has a solution with cost 9·E X ≤ 2.2519 | `lemmaC/round3/final/TREE_AND_LOCALITY.md`, `rdeB/ROUTES.md`, `cert/` | §6–7 |
| 5. Probability bounds | Gaussian concentration of the square-root deficits and an event that forces one integral optimum | `manuscript/sections/10_proofmain.tex` | §10 |

## Where to start

1. [Read the 117-page working paper](manuscript/build/main.pdf), dated 17 November 2026. Sources in `manuscript/`.
2. `lemmaC/round3/MAIN_THEOREM.md` — proof map for the original theorem.
3. `lemmaC/round3/CLAIMS.md` — register of every lemma and the checks it received.
4. `lemmaC/round3/final/ERRATA.md` — historical corrections to the component proofs.
5. `lemmaC/round3/codex/` — the Codex reads, verbatim.

## Computer-assisted steps

Two inequalities are certified with interval arithmetic:

- `lemmaC/round3/final/certI/certify_rigorous.py` — Ψ₀ ≤ −3.84751328322 in the independent-team model. Including the correction for shared pair values certifies 2.302 in the market.
- `lemmaC/round3/cert/cert.py` — 9·E X ≤ 2.2519 for every fixed point of the tree recursion.

mpmath 1.3.0's interval `exp` and `log` were found not to be strict enclosures; both
certificates use verified replacements. The trusted base is mpmath's four basic interval
operations and IEEE arithmetic in the summations. Run with Python 3, numpy, scipy and
mpmath 1.3.0; commands and expected output are in `lemmaC/round3/cert/CERTIFICATE.md`,
`lemmaC/round3/final/certI/NOTES.md` and Appendix B of the manuscript.

The independent Arb implementations verify the corrected integer exponent at 2.302 and a tree cost below 2.250535. Together they give the stronger lower bound of 0.051465 on the limiting inferior of the expected subsidy divided by *q*^(1/3), as recorded in §11. The [proof audit and reproduction commands](lemmaC/independent/README.md) refer to the preceding 116-page draft. The current working paper also gives Corollary 10.3 on fractional assignments and clarifies the transfer to shared pair values in §§4–5.

## Supporting evidence (floating point, not part of the proof)

- `gate1/`, `gate1_ext/` — simulations of the nonempty-core share (0.67 at q = 10, 0.008 at q = 40) and exact empty-core certificates at q = 3.
- `lemmaC/round3/pscout/`, `lemmaC/round3/redteam2/` — large-q constants in the independent-team model and the market.
- `lemmaC/round3/v1/` … `v16/` — scripts written by the adversarial readers.
- `lemmaC/` (other folders), `THEOREM_VANISHING.md`, `vanishing/` — earlier rounds, failed routes and partial results.
- `LEDGER.md`, `PLAN.md`, `GATE0_REPORT.md` — the dated record, the frozen plan and the prior-art sweep.

## Caveats

- The probability lower bound has an explicit finite-size formula in §10. The upper bound uses an asymptotic separation, with some thresholds in the proof exceeding 10^50,000.
- Imports Frieze, Pegden and Tkocz, *On random multi-dimensional assignment problems*, Discrete Appl. Math. 287 (2020), Theorem 2(b).
- Files reference absolute paths on the author's machine (`/Users/maxf/...`) and relative paths of the form `lemmaC/round3/...`; the latter resolve from this repository's root.
- The canonical tree certificate checks its final bound against an exact rational. The coarse and fine stage files are in `lemmaC/round3/cert/stages/`.
