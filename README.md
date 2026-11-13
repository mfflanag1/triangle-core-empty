# The core of a random three-sided assignment market is empty

## The result

Three sectors of *q* agents. Each of the 3*q*² pairs of agents from different sectors has an
independent U[0,1] value; a team (one agent per sector) is worth the sum of its three pair
values. *T* is the best total value of a perfect assignment into teams, LP\* the optimum of
the fractional relaxation, and *D* = LP\* − *T*. The core of the market is nonempty exactly
when *D* = 0.

**Theorem.** As *q* → ∞,

- E *D* ≥ (0.050 − o(1))·*q*^(1/3),
- *D* / E *D* → 1 in probability, and
- P(core nonempty) → 0, in fact O(*q*^(−5/6+ε)).

## Proof in four parts

| Part | Statement | Source file | Manuscript |
|---|---|---|---|
| 1. Variance | Whp no used pair or LP-support team has deficit above C·q^(−1/2+ε); Var *D* = O(q^(−1/6+ε)) | `lemmaC/round3/final/THEOREM_1.md` | §3 |
| 2. Integer bound | Whp 3q − T ≥ 2.302·q^(1/3) (independent-team model, then transfer) | `lemmaC/round3/firstmoment/THEOREM_I.md`, `sc/TRANSFER.md`, `interface/THEOREM_I_PA.md` | §4–5 |
| 3. LP bound | limsup E(3q − LP\*)/q^(1/3) ≤ c(m) for any invariant matching m on the Poisson hypertree | `lemmaC/round3/localconv/LOCALCONV.md`, `lplocal/LP_LOCALITY.md` | §8–9 |
| 4. Existence and cost | The tree recursion has a solution with cost 9·E X ≤ 2.2519 | `lemmaC/round3/final/TREE_AND_LOCALITY.md`, `rdeB/ROUTES.md`, `cert/` | §6–7 |

## Where to start

1. `manuscript/build/main.pdf` — the consolidated paper (121 pages). Sources in `manuscript/`.
2. `lemmaC/round3/MAIN_THEOREM.md` — statement and proof map.
3. `lemmaC/round3/CLAIMS.md` — register of every lemma and the checks it received.
4. `lemmaC/round3/final/ERRATA.md` — corrections applied and still owed.
5. `lemmaC/round3/codex/` — the Codex reads, verbatim.

## Computer-assisted steps

Two inequalities are certified with interval arithmetic:

- `lemmaC/round3/final/certI/certify_rigorous.py` — Ψ₀ ≤ −3.84751328322 (integer constant ≥ 2.30224586).
- `lemmaC/round3/cert/cert.py` — 9·E X ≤ 2.2519 for every fixed point of the tree recursion.

mpmath 1.3.0's interval `exp` and `log` were found not to be strict enclosures; both
certificates use verified replacements. The trusted base is mpmath's four basic interval
operations and IEEE arithmetic in the summations. Run with Python 3, numpy, scipy and
mpmath 1.3.0; commands and expected output are in `lemmaC/round3/cert/CERTIFICATE.md`,
`lemmaC/round3/final/certI/NOTES.md` and §11 of the manuscript.

## Supporting evidence (floating point, not part of the proof)

- `gate1/`, `gate1_ext/` — simulations of the nonempty-core share (0.67 at q = 10, 0.008 at q = 40) and exact empty-core certificates at q = 3.
- `lemmaC/round3/pscout/`, `lemmaC/round3/redteam2/` — large-q constants in the independent-team model and the market.
- `lemmaC/round3/v1/` … `v16/` — scripts written by the adversarial readers.
- `lemmaC/` (other folders), `THEOREM_VANISHING.md`, `vanishing/` — earlier rounds, failed routes and partial results.
- `LEDGER.md`, `PLAN.md`, `GATE0_REPORT.md` — the dated record, the frozen plan and the prior-art sweep.

## Caveats

- All bounds are asymptotic; some thresholds in the proof exceed 10^50,000.
- Imports Frieze, Pegden and Tkocz, *On random multi-dimensional assignment problems*, Discrete Appl. Math. 287 (2020), Theorem 2(b).
- Files reference absolute paths on the author's machine (`/Users/maxf/...`) and relative paths of the form `lemmaC/round3/...`; the latter resolve from this repository's root.
- Eight items remain marked `\todo` in the manuscript (bibliographic details, three small changes to `cert.py`).
