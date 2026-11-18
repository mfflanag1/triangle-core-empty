# Conventions for the consolidated manuscript

Every section file follows this. Macros are in `macros.tex`; use them and do not redefine.

## What the manuscript is

The original round-3 proof files under `../lemmaC/round3/` document the November 12 draft. The revised manuscript includes two bounds on the probability of exact core stability and an attained variational limit for the normalized fractional deficit. The manuscript is the source for these results.

## Setting and symbols

- Three sectors A, B, C of q agents; agents a_i, b_j, c_k. A *pair* is two agents from
  different sectors; there are 3q² pairs with independent U[0,1] values; the *deficit* of a
  pair p is `\defi(p)` = 1 − value. A *team* e = (a, b, c) has value the sum of its three
  pair values and deficit `\Defi(e)` = sum of its three pair deficits.
- A *perfect assignment* is q disjoint teams covering all agents. `\Tq` is the best total
  value, `\Mopt` the optimal assignment (a.s. unique), `\Phiq` = 3q − `\Tq` its deficit.
- The fractional relaxation: weights x ≥ 0 on teams with every agent's load equal to 1;
  optimum `\LPq`; `\PhiLP` = 3q − `\LPq`. Gap `\Dq` = `\LPq` − `\Tq` = `\Phiq` − `\PhiLP` ≥ 0.
- The core: payoffs p ≥ 0 on agents with p(e) ≥ value(e) for every team and total `\Tq`.
  Nonempty iff `\Dq` = 0 (LP duality). State this once, in the model section.
- Rescaled cost `\cost_e` = q^(2/3)·`\Defi(e)`. Near-tight hypergraph H_q(K): teams with
  `\cost_e` ≤ K, marked by their costs.
- Independent-team model: the q³ teams have iid deficits, each the sum of three independent
  U[0,1] variables. Call it "the independent-team model" (never "Poissonised").
- Poisson hypertree `\tree`: each agent has teams whose costs form a Poisson process of
  intensity `\Lam`(dξ) = ξ²/2 dξ; each team has three agents; the incidence graph is a tree.
  "Invariant" means the team-level mass-transport identity (MT) of
  `final/TREE_AND_LOCALITY.md` Definition 1.1.
- The recursion map `\RDE`: `\RDE`(μ) is the law of min_k(ξ_k − X_k′ − X_k″), ξ a Poisson
  process of intensity `\Lam`, X′, X″ iid μ. (The sources call it T or Θ; use `\RDE`.)
  For G = −ln P(X > ·): G(x) = E[(x + X′ + X″)_+³]/6.
- Dispersive order `\disp` on laws written as Q(E), E ~ Exp(1): Q_0 `\disp` Q_1 iff
  Q_1 − Q_0 is nondecreasing. IFR = increasing failure rate.
- Constants: `\betaint` = 2.302 (integer bound; the certificate gives 2.30224586);
  `\alphalp` = 2.2519 (tree cost bound; numerically the cost is 2.2495456). c_2 is the
  Frieze–Pegden–Tkocz constant and C_0 = sup_n E Φ_n / n^(1/3).

## Results and labels (use exactly these labels)

| Label | Statement |
|---|---|
| `thm:main` | E `\Dq` ≥ (`\betaint` − `\alphalp` − o(1))·q^(1/3); `\Dq`/E`\Dq` → 1 in probability; exp(−Cq) ≤ P(core nonempty) ≤ exp(−cq^(1/3)) |
| `thm:variance` | whp no used pair above C_ε q^(−1/2+ε), no LP-support team above C_ε q^(−1/2+ε); Var `\Dq` = O(q^(−1/6+ε)) |
| `thm:iid` | independent-team model: whp every perfect assignment has deficit ≥ `\betaint`·q^(1/3) |
| `thm:transfer` | comparison inequality between the market and the independent-team model |
| `thm:integer` | the market: whp `\Phiq` ≥ `\betaint`·q^(1/3) |
| `thm:localconv` | total-variation convergence of marked balls of H_q(K) to those of `\tree` |
| `thm:lp` | limsup E `\PhiLP`/q^(1/3) ≤ c(m) for every invariant perfect matching m on `\tree` |
| `thm:lp-limit` | The normalized fractional deficit converges in mean square to the attained minimum invariant matching cost on the tree. |
| `thm:decomp` | invariant fractional perfect matchings on `\tree` are marginals of invariant random integral ones |
| `thm:process` | a fixed point of `\RDE` gives an invariant perfect matching of cost 9·E X |
| `thm:mono` | `\RDE` is monotone in `\disp` on IFR laws |
| `thm:existence` | `\RDE` has a fixed point, and every fixed point has 9·E X ≤ `\alphalp` |

Section labels: `sec:intro`, `sec:model`, `sec:variance`, `sec:iid`, `sec:transfer`,
`sec:tree`, `sec:existence`, `sec:localconv`, `sec:lp`, `sec:proofmain`, `sec:certificates`,
`sec:discussion`. Label lemmas `lem:<section>-<name>`, equations `eq:<section>-<name>`.

## Style

- LaTeX, `amsart`-style environments: theorem, lemma, proposition, corollary, definition,
  remark (all defined in `main.tex`; do not define them). One file per section in
  `sections/`, starting with `\section{...}\label{sec:...}`; no preamble.
- Complete proofs, not sketches. Every constant and every union bound explicit, as in the
  sources. Numerical evidence goes in remarks, clearly labelled as floating point.
- Computer-assisted steps: state exactly what the computation certifies and which script
  does it (path relative to the lane folder); summaries go to `sec:certificates` and reproduction details to Appendix B (`sec:certificates-detail`).
- British spelling. No first-person singular. Plain sentences.
- Use the bibliography keys defined in `main.tex` and cite the version in which the stated result was checked.
- Describe the actual scope of each verification.
