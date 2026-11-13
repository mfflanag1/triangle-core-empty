# Gate 0 — Average-case protection (random multi-sided pair-additive markets)

Sweep date: 2026-11-09. No experiments were run; nothing under `/Users/maxf/Documents` or `/Users/maxf/projects` was modified; no research-handoff receipt was enqueued (per instruction).

## 1. Verdict

**OPEN, but only narrowed. Do not open it as proposed.**

No paper studies protection of individual core maxima, or a floor-constrained subsidy, in random multi-sided markets. But the two headline questions are close to answered by prior art plus short arguments, and both P1 and P7 quantities are undefined on empty-core instances, which is the generic risk for random cyclic markets.

| Headline question | Status after the sweep |
|---|---|
| Typical ρ on trees | Not in the literature, but it follows from Hassidim–Romm (2015) plus P1's own cut lemma: ρ(M) → 1/K with high probability. The worst-case tree constant is also the typical value. (My deduction, unverified.) |
| Does the tree/cycle divide in ρ survive? | Probably not in ρ. With equal sector sizes, sector shifts make every H_i close to its team's value, so ρ ≈ 1/K on any graph whenever the core is nonempty and payoffs within a sector are nearly equal. The divide moves to core existence: always for trees, unknown for cycles. (My deduction, unverified.) |
| Does s(M)/T vanish, and how fast? | Yes, trivially, for bounded iid values: s(M) ≤ 3q − T, and Frieze–Pegden–Tkocz (2020) give E[3q − T] = Θ(q^{1/3}), so s/T = O(q^{-2/3}). The exact order of s is open; it may be 0 with high probability. |

**Narrow to:**
- **X1.** Core-existence probability π_G(q) and LP gap for random pair-additive cyclic markets. Only iid-triple arrays at tiny sizes are known (Lucas 1995). This is the real open problem.
- **X2.** A multi-sided law of one price (within-sector payoff dispersion, plausibly depending on the sector's degree in G), with ρ → 1/K as a corollary for every G.
- **X3.** Second order on the triangle: P(s(M) = 0) = P(ρ(M) ≥ 1/3), given a nonempty core.
- **X4.** Optionally, a model where the quantities are not degenerate: unbounded values, where Hassidim–Romm show one price fails.

**Kill condition.** If the first experiment (section 6) shows π_triangle(q) collapsing to 0 by q ≈ 10, then ρ and s as defined in the manuscripts are vacuous for typical cyclic markets. The lane then survives only if restated on the LP-dual optimal face (the extended/aspiration core).

## 2. Closest prior results

### (a) Hassidim and Romm — random bilateral assignment games
"An Approximate 'Law of One Price' in Random Assignment Games." [arXiv:1404.6103](https://arxiv.org/abs/1404.6103) (v1, 2014); two-page abstract in EC'15, pp. 279–280, [doi:10.1145/2764468.2764531](https://doi.org/10.1145/2764468.2764531). **Full text read** (arXiv v1).

- **Model.** n firms, n + k(n) workers, pair values iid on [0,1] with continuous, strictly positive density.
- **Theorem 1.** With probability ≥ 1 − O(1/n), in every core allocation all firms' profits lie within c·log²n/n of each other, and likewise all workers' salaries.
- **Corollary 7 (balanced).** Under the firm-optimal core allocation, the workers' share of surplus is O(log²n/n).
- **Theorem 8 (U[0,1], balanced).** The expected sum of workers' salaries at the firm-optimal point is ≤ log n and Ω(log n).
- **Proposition 10.** For exponential values the law fails: two workers differ by Ω(log n).
- **Consequence for us (my deduction).** For K = 2, every H_i = 1 − O(log²n/n), so ρ ≤ T/ΣH → 1/2, while ρ ≥ 1/2 always.
- **Leaves open.** K ≥ 3, full-team production, cyclic sector graphs, any protection ratio or subsidy, the conjectured tight log n/n rate, unbounded distributions.

### (b) Lucas, reporting Crowley and Graham — random three-sided cores
W. F. Lucas, "Core theory for multiple-sided assignment games," Duke Math. J. 81 (1995) 55–65, [doi:10.1215/S0012-7094-95-08106-X](https://doi.org/10.1215/S0012-7094-95-08106-X). **Full text read** for pp. 63–64, from the local copy `/Users/maxf/Documents/ECON/artifacts/source_comparison_2026-10-31/text/lucas1995.txt`. The primary sources (Crowley–Graham Stony Brook talk 1994; Crowley's 1995 and Graham's 1996 Claremont dissertations) were **not read**.

Payoffs a_ijk are iid uniform on [0,1] over triples, not pair-additive.

| Size | P(core nonempty) |
|---|---|
| 2×2×2 | 209/240 = 0.8708 (exact) |
| 3×3×3 | 73.1% |
| 2×2×r | > 95% for r ≥ 7, > 99% for r ≥ 40 (simulation) |
| 2×5×5 | 96.4% |
| 2×20×20 | 98.9% |
| p fixed, q and/or r → ∞ | → 1, with gains going to the short side |
| 2^m (two per sector, m sectors) | → 0 as m → ∞ |

- **Leaves open.** Pair-additive values; the balanced q×q×q limit (only q = 2, 3 are stated, and the probability falls); anything about payoffs.
- **Correction to the brief.** The "probabilities of nonempty cores" result is Lucas 1995, not Marchi–Oviedo. Marchi and Oviedo, "The core of a further m-sided assignment game," EJOR 98 (1997) 617–625, is a deterministic generalization (**abstract-level only**, via the IAOR/IDEAS listing).

### (c) Frieze, Pegden and Tkocz — random pair-additive 3-D assignment
"On random multi-dimensional assignment problems," Discrete Applied Mathematics 287 (2020) 1–9, [doi:10.1016/j.dam.2020.07.013](https://doi.org/10.1016/j.dam.2020.07.013), [arXiv:1901.07167](https://arxiv.org/abs/1901.07167). **Full text read** (arXiv v1).

- **Model.** W_ijk = a_ij + b_ik + c_jk with a, b, c iid U[0,1]; Z_n is the minimum of Σ_i W_{i,σ(i),τ(i)} over permutations. This is exactly the random triangle market in minimization form.
- **Theorem 2.** E(Z_n) ≥ c₁n^{1/3}, and a greedy algorithm achieves expected value ≤ c₂n^{1/3}.
- **Theorem 5.** The d-dimensional analogue with (d−1)-wise decomposable costs scales as n^{1/d}. The four-sector pairwise-sum model is explicitly "not considered."
- **Consequence for us (my deduction).** Flipping a ↦ 1 − a gives E[3q − T] = Θ(q^{1/3}) for the triangle with U[0,1] values.
- **Leaves open.** The constant; the LP relaxation and its duality gap (that is, core existence); dual variables; K ≥ 4 pairwise models. The paper contains nothing on LP or duality.

### Secondary items

| Source | Read | Relevance |
|---|---|---|
| Kanoria, Saban, Sethuraman, "The size of the core in assignment markets," [arXiv:1407.2576](https://arxiv.org/abs/1407.2576); Operations Research 66 (2018) 620–636 | Full text (arXiv v1, through the proof overview) | Typed bilateral model; expected core size O*(n^{-1/max(K,Q)}) under their Assumption 1 (no balanced sub-market). Remark 1: a balanced one-type market has core size ≥ u with probability 1. Equal sector sizes sit exactly in their excluded case, so core shrinkage does not transfer to P1/P7. |
| Mordant, "A central limit theorem for the random assignment problem," [arXiv:2608.05123](https://arxiv.org/abs/2608.05123) (Aug 2026) | Abstract only | Gives an explicit law for an optimal dual potential (ordered potential gaps become independent exponentials). A likely tool for exact dispersion rates at K = 2 and on trees. |
| Yokote, "Constrained optimal transport with applications to large matching markets," [arXiv:2604.02559](https://arxiv.org/abs/2604.02559) (2026) | Abstract only | TU matching with couples: core may be empty in finite markets but dual optimizers support a nonempty core in the continuum. |
| Kaneko and Wooders, "Cores of partitioning games," Math. Soc. Sci. 3 (1982) 313–327 | Abstract only | Approximate cores are nonempty in large games with few types and small essential coalitions. Deterministic precedent for a vanishing per-capita cost of stability. |
| Vazirani, [arXiv:2402.11437](https://arxiv.org/abs/2402.11437); Gangam et al., AAMAS 2026 | Abstracts only | Leximin and equitable core imputations; worst-case and algorithmic, no random instances. |
| Martin, Mézard, Rivoire (cond-mat/0407623, 0507180); Frieze and Sorkin, RSA 46 (2015) | Abstracts only | iid costs on d-tuples; nothing on decomposable costs or duals. |

Forward-citation lists of Kanoria–Saban–Sethuraman, Hassidim–Romm, Frieze–Pegden–Tkocz, Lucas, Quint, Stuart and Atay–Núñez contain no multi-sided or average-case core follow-up.

## 3. Exact definitions

**P1** (`acyclicity_ijgt/submission/manuscript.tex`, Section 2).
- G = (V, E) is a connected simple graph with K = |V| ≥ 2 sectors, each with m players (v, j).
- Each allowed edge uv carries a nonnegative matrix A_uv; zero matrices are allowed.
- A full team t = (j_v)_{v∈V} has worth w(t) = Σ_{uv∈E} A_uv(j_u, j_v). Only full teams produce.
- v(S) is the maximum total worth of player-disjoint full teams inside S. T = v(N) > 0, and the core is required to be nonempty.
- 𝒞 = {p ≥ 0 : Σ_i p_i = T, Σ_{i∈t} p_i ≥ w(t) for every full team t}, and H_i = max_{y∈𝒞} y_i.
- x is the full ascending-sorted leximin point of 𝒞.
- c(G) = inf_M min_{i:H_i>0} x_i/H_i; **ρ(M) = max_{p∈𝒞} min_{i:H_i>0} p_i/H_i**; R(G) = inf_M ρ(M), over markets on fixed G with m free.
- Theorem 2.1: c(G) = R(G) = 1/K for trees and 0 for cyclic G; the triangle family has ρ(M_q) = 1/(2q+3).

Notation clash: in P1, q indexes the family and the sector size is m = 3q + 1. Below, q means agents per sector, as in the brief.

**P7** (`static_subsidies_ijgt/submission/manuscript.tex`, Section 2).
- Sectors N₁, N₂, N₃ of size n; nonnegative A on N₁×N₂, B on N₁×N₃, C on N₂×N₃.
- v(e) = A_{i₁i₂} + B_{i₁i₃} + C_{i₂i₃} for a full team e; v(Q) is the maximum over player-disjoint full teams in Q.
- T = v(N) > 0; the original core 𝒞 = {x ∈ ℝ^N_+ : x(N) = T, x(e) ≥ v(e) for all e} is assumed nonempty; H_i = max_{x∈𝒞} x_i.
- **s(M) = min{ y(N) − T : y_i ≥ H_i/3 for all i, y(e) ≥ v(e) for every full team e }.**
- Dual: s(M) = max_{λ∈𝒫}(W + S/3 − T) over fractional partial matchings.
- Theorem 2.2: s(M) < T/3 and sup s/T = 1/3.

**Link.** For three sectors, s(M) = 0 if and only if ρ(M) ≥ 1/3. On a path (C ≡ 0), s ≡ 0 by P1's tree theorem.

**Both quantities need a nonempty core.** A random triangle market with continuous values has one exactly when the LP relaxation of the axial 3-D assignment problem has an integral optimum.

## 4. What the sweep implies (derived here, not from the literature, unverified)

1. **Sector-shift bound.** With equal sector sizes, adding constants c_u to each sector with Σc_u = 0 preserves all team sums and the total, so it preserves the core while payoffs stay nonnegative. Hence for any core point p and any player (r, j) with optimal-matching team value d_j:
   H_{r,j} ≥ d_j − Σ_{u≠r} [p_u(j) − min_i p_u(i)].
   Also ρ(M) ≤ T/Σ_i H_i always, and Σ_i H_i ≤ K·T.
2. **Trees.** P1's cut inequalities U_b(i) − U_b(j) ≥ A_ab(j,i) − d_ab(j) are the bilateral stability inequalities on edge ab. Hassidim–Romm's path argument uses only these, so within-sector dispersion is O(deg·log²q/q) with high probability. Then H_i ≥ d_j − O(K·log²q/q), so for fixed K, ρ ≤ 1/K + O(log²q/q), against ρ ≥ 1/K from P1. Needs a check that their balanced-case proof uses no sign constraints; I read it as using only α ≤ 1 and stability.
3. **Any G, given a nonempty core and small dispersion ε.** The same bound gives ρ = 1/K ± O(ε). For the triangle, replacing one player in a tight team suggests ε = Õ(q^{-1/2}), and more generally Õ(q^{-1/deg(v)}). The matching is adversarial, as in Hassidim–Romm, so this needs a real proof.
4. **Subsidy.** y ≡ 1 covers every team (values ≤ 3) and meets every floor (H_i ≤ d_j ≤ 3), so s(M) ≤ 3q − T.
5. **Two per sector.** For any nonnegative pair-additive 2×2×2 triangle market the core is nonempty: the two fractional vertices of the 2×2×2 axial polytope each have value S/2, where S is the sum of all twelve pair entries, and S/2 is also the average of the four perfect matchings. So π(2) = 1, against 209/240 for iid triples. The repository's September pilot agrees: 190 of 190 nonzero m = 2 instances had a core, and 77 of 78 at m = 3 (small integer alphabets; the one empty case is float-flagged, not exactly certified).

## 5. Sweep log

- **Crossref by ISSN, from 2024-01-01, rows = 50.** ISSNs looked up through the Crossref journals API.
  - Requested eight: IJGT 0020-7276; GEB 0899-8256; JET 0022-0531; MOR 0364-765X; JME 0304-4068; SCW 0176-1714; MSS 0165-4896; Economic Theory 0938-2259.
  - Added: Annals of OR 0254-5330; Discrete Applied Mathematics 0166-218X; Operations Research 0030-364X; Theoretical Economics 1933-6837; Random Structures and Algorithms 1042-9832; JAIR 1076-9757; Management Science 0025-1909; Econometrica 0012-9682.
  - 16 queries per journal (256 calls, no failures, 1,856 records title-screened): assignment game core; multi-sided assignment; random market core; cost of stability; large market core matching transfers; three-sided; least core; egalitarian core assignment; assignment; core; nonempty core probability; hypergraph matching; supplier firm buyer; buyer-optimal seller-optimal; asymptotic core convergence; matching transfers random.
  - Nearest hits were not relevant (NTU random markets; Martínez-de-Albéniz–Rafels–Ybern 2026, already cited in P1/P7).
- **arXiv API, 28 queries, all dates plus 2023/2024–2026 windows.** Random assignment game and core; multi-/three-/m-sided assignment; cost of stability with random/average/expected; multidimensional and three-dimensional assignment with random; multi-index matching; random cooperative games; least core random; matching games core random; random assignment problem dual; supplier-firm-buyer; matching for teams.
- **OpenAlex.** Title-and-abstract search from 2023-01-01, 27 queries, 146 unique records across all venues. Per-venue title sweep from 2023 for EC, STOC, FOCS, AAMAS and AAAI (terms: core, assignment, matching, coalition, cooperative, stability, transferable). IJCAI and one AAAI term were rate-limited there; SODA returned nothing under its OpenAlex source IDs.
- **Crossref by proceedings title, from 2023.** EC, SODA, STOC, FOCS, AAAI, AAMAS, IJCAI with 7 queries each (687 records); WINE and SAGT as book chapters (198 chapters, including near-complete SAGT 2023–2026 contents and three WINE volumes). No average-case TU-core paper.
- **Semantic Scholar forward citations.** Kanoria–Saban–Sethuraman (arXiv and OR), Hassidim–Romm (arXiv and EC), Frieze–Pegden–Tkocz, Frieze–Sorkin, Lucas 1995, Quint 1991, Stuart 1997, Atay–Núñez 2019, Martínez-de-Albéniz et al. 2020. Its keyword search was rate-limited (2 of 7 queries returned).
- **Web search.** About 12 queries on: Marchi–Oviedo probabilities; random three-sided core; Crowley–Graham; law of one price multi-sided; cost of stability average-case; random 3-D assignment LP integrality; random TU game core probability; multi-marginal transport and matching for teams; Kaneko–Wooders; Hassidim–Romm follow-ups; LP tightness for multi-way matching.
- **Gaps.**
  - DBLP is behind an anti-bot challenge on all three mirrors, so proceedings coverage is by title keywords through OpenAlex and Crossref, not DBLP.
  - Not searched: Google Scholar, SSRN full text, ProQuest (the Crowley and Graham dissertations).
  - Venue sweeps are keyword-driven, so a paper with an unusual title could be missed.

## 6. Local code notes

Repository: `/Users/maxf/projects/economics-research/game-theory/coopgpt`, read only. Interpreter `/Users/maxf/miniforge3/bin/python3` has scipy 1.15.3, numpy 2.2.6, z3 5.0.0, ortools 9.15 and pulp; highspy is absent.

- **s(M), general input.** `structural_classification_2026-09-06/sharp_subsidy/discovery/probe.py`, function `lp(matrices, m, known_T=None, known_H=None)`.
  - scipy `linprog(method='highs')`, floating point. Builds all m³ team rows, tests core nonemptiness by comparing the fractional cover value with T, computes H_i with 3m LPs, then s with one LP.
  - T is brute force over pairs of permutations, (m!)², unusable beyond m ≈ 4–5.
  - It imports `certificates` from a `/private/tmp/coopgt-simultaneous-20260907` path that no longer exists; the same module is at `structural_classification_2026-09-06/simultaneous_protection/certificates.py`.
  - `discovery/search.py` is a seeded random search for the worst case, not a census.
- **P1 quantities, general graph.** `paper5/graph_core_protection/verification/tree_extension/raw_core.py`, function `solve(edges, matrices, m, k=4)`.
  - Exact rationals with Z3 `Optimize`; T by brute force over permutation tuples. Returns T, the leximin point x, H and `ratio = min x_i/H_i`.
  - That ratio is the leximin quantity behind c(G), not ρ(M). **No general ρ(M) solver exists in the repository**; it needs one extra LP (maximize t subject to p ∈ 𝒞, p_i ≥ t·H_i).
  - The paper's exact ρ values are closed-form certificates in `verification/simultaneous_protection/{certificates.py, verify.py, family.py}`; the package entry point is `verification/run_checks.py` via `build.sh [--smt]`.
- **Triangle pilot.** `pair_additive_campaign_2026-09-05/experiment/pilot.py`, function `float_check(A, B, C, label)`.
  - scipy HiGHS; returns status `empty_core_float` or `nonempty_core_float`, both leximin selectors, H and `min_fraction`. The exact cross-check is `verification_dependency/exact_three.py` (`core`, `check`, `lucas_test`; Z3).
  - **Naming hazard:** `rho` in this file is ‖Lmin − Lmax‖₁/T, not P1's ρ(M).
  - Existing outputs `pilot_summary.json` and `screen_records.json` hold the counts quoted in section 4.
- **P7 submission package.** `paper5/sharp_subsidy/verification/{verify.py, verify_sharp_aggregation.py, verify_endpoint.py, tree_certificate.py}` are certificate checkers only, standard library, with no general solver.
- **For a simulation gate.** Reuse the row construction and the H and s LPs from `probe.py::lp`. Replace the brute-force T with one primal LP (maximize Σλ_e·v(e) over fractional matchings): with continuous values the optimum is almost surely unique, so the core is nonempty exactly when it is integral, and then T equals the LP value. Confirm with scipy `milp` or ortools at small q. Add the one-line ρ LP.

## 7. Sharpened target and first experiment

**Target.** Take a connected sector graph G on K sectors, q agents per sector, pair values iid from a bounded law with continuous density positive at its supremum (benchmark U[0,1]), strict full-team production, and coalition worth by integral matching. The claim to establish is that the worst-case tree/cycle divide of P1 does not survive in ρ but reappears in core existence:
- (i) with high probability, within-sector payoff dispersion on the LP-dual optimal face is small, at most Õ(q^{-1/deg_G(v)}) for sector v;
- (ii) whenever the core is nonempty, ρ(M) = 1/K + O(dispersion) for every G, so the sharp worst-case tree constant is the typical value on trees and cycles alike;
- (iii) π_G(q) = P(core nonempty) is 1 for trees, 1 for the triangle at q = 2, apparently below 1 from q = 3 (one float-flagged pilot instance), with its large-q behaviour and the relative gap (T_LP − T)/T = O(q^{-2/3}) to be determined;
- (iv) the P7 subsidy satisfies s/T = O(q^{-2/3}) trivially, so the question worth asking is whether P(s(M) = 0 | core nonempty) tends to 1.

Unbounded values (exponential), where Hassidim–Romm show one price fails, are the natural non-degenerate contrast.

**Cheapest decisive experiment.** Triangle, iid U[0,1], q ∈ {2, 3, 4, 5, 6, 8, 10, 12, 16, 20}, about 2,000 instances per q (fewer at the top).

1. Solve the primal LP and record integrality, giving π̂(q), and the gap (T_LP − T)/T at small q using MILP for T.
2. On nonempty-core instances compute H_i, ρ (one extra LP), s, and the maximum within-sector spread at the ρ-optimal point.
3. Controls: the path (C ≡ 0), which must give π = 1, ρ ≥ 1/3, s = 0; q = 2, which must give π = 1; exponential values as contrast.

Cost: LPs with at most 8,000 variables or rows; minutes on a laptop.

Decision rule:
- π̂(q) below about 0.05 by q ≈ 10–12: the manuscripts' definitions are vacuous for typical cyclic markets. Restate on the extended core or drop the lane.
- π̂ bounded away from 0, with ρ̂ → 1/3 and spread → 0: the paper is "multi-sided law of one price plus core-existence probability," and protection is a corollary.
- ρ̂ staying bounded below 1/3: my deduction in section 4 is wrong and the lane opens as originally proposed.
