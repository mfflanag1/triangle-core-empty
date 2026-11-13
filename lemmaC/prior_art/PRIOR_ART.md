# Prior art for Lemma C / Conjecture V (random pair-additive 3-index assignment, LP integrality)

Written 2026-11-11 by a Claude Code subagent in a 45-minute time box. No API spending, no outreach.
Complements GATE0_REPORT.md section 2 (Hassidim–Romm, Lucas 1995, Frieze–Pegden–Tkocz, Kanoria et al.,
Martin–Mézard–Rivoire, Frieze–Sorkin, Yokote, Mordant); those are not repeated except where a new
point bears on Lemma C.

## Verdict

**OPEN (not known, not partly known in the required form). Confidence about 80%.**

- No source found states P(LP integral) -> 0 or -> 1, or P(core nonempty) -> 0 or -> 1, for random
  axial 3-index assignment, with or without pair-additive weights. Frieze–Pegden–Tkocz (the only paper on
  the exact pair-additive model) has no LP or duality content (GATE0, full text read). Its two forward
  citations in OpenAlex (Cheng–Liu–Tkocz–Xu 2023; a 2022 IEEE TSMC role-assignment paper) do not touch it.
- **Not even the cleanest analogue is settled rigorously.** For min/max-weight perfect matching on the
  complete graph K_n with iid edge weights, I found no result, rigorous or physical, on the probability that
  the fractional matching LP optimum is integral. What exists: (i) the LP value and the integer value have
  the same limit, π²/12; (ii) for exponential weights, the expected difference is O(log n / n), smaller than
  the standard deviation n^{-1/2} of either (Hessler–Wästlund, rigorous); (iii) physics (replica/cavity)
  computes the 1/N finite-size difference and attributes it to odd cycles in the optimal solution
  (Lucibello–Malatesta–Parisi–Sicuro 2018; Parisi–Perrupato–Sicuro 2020).
- **Inference (mine, not in any source).** On K_n the fractional/integer cost gap is Θ(1/N) on a per-edge
  cost scale of 1/N. That is consistent with O(1) half-integral odd cycles in the LP optimum, and so with
  P(LP integral) tending to a constant strictly between 0 and 1, not to 0. This contrasts with the lane's
  numerics for the 3-index problem: the LP-optimal vertex disrupts about 0.8q of the M-teams (F2), and
  P(empty) rises towards 1 (F3). The 2-D non-bipartite analogue is therefore **not** expected to supply the
  answer by transfer; the 3-index phenomenon looks qualitatively different (mesoscopic, not local).
- **Search gaps.** Google Scholar was not queried directly. The arXiv API was rate-limited ("Rate
  exceeded"), so arXiv was reached through web search and OpenAlex only. ProQuest (Crowley and Graham
  dissertations) was not searched. Physics papers were full-text grepped, not read in full.

## Three most useful sources

1. **Hessler and Wästlund, "Concentration of the cost of a random matching problem"** (preprint dated
   20 August 2008, `math.chalmers.se/~wastlund/martingale.pdf`; OpenAlex lists 2009 with no venue). This
   is the only rigorous concentration result found for an LP relaxation of a random non-bipartite
   matching problem. Var(C_n) = O(1/n) for the LP-relaxed ("friendly model") matching cost uses an
   oracle/exploration process plus a correlation inequality (Lemma 5.1, Lemma 7.1). It also shows
   E|M_n − π²/12| < 1.706/√n, and E(M_n − C_n) = O(log n/n), smaller than either standard deviation.
   This is the template warning for Lemma C: LP* and T can each concentrate while their gap stays
   unresolved, so the signed margin mu* is the right scalar.
2. **Lucibello, Malatesta, Parisi and Sicuro, "The random fractional matching problem"**, J. Stat. Mech.
   (2018) 053301, doi:10.1088/1742-5468/aabbc8, arXiv:1802.02997. This is the explicit 2-D non-bipartite
   analogue of LP versus IP: the same π²/12 limit, with different 1/N corrections. Read it together with
   **Parisi, Perrupato and Sicuro, "Random-link matching problems on random regular graphs"**, J. Stat.
   Mech. (2020) 033301, doi:10.1088/1742-5468/ab7127, arXiv:1910.03504. That paper gives a cycle
   expansion: the integer-minus-fractional difference is carried by odd cycles, each contributing
   positively. These papers supply the cavity language for "fractional structures relative to M". They
   say nothing on P(integral).
3. **Haber and Krivelevich, "On fractional K-factors of random graphs"**, Random Structures & Algorithms
   30(4) (2007) 441–463, doi:10.1002/rsa.20144. Read it with **Johansson, Kahn and Vu, "Factors in random
   graphs"**, RSA (2008), doi:10.1002/rsa.20224. For strictly K1-balanced K (K3 included), a fractional
   K-factor appears w.h.p. at the moment every vertex is covered by a copy of K; JKV give the integral
   triangle-factor threshold at the same scale, with a log factor. This rigorously backs the F1 heuristic
   that LP* and T lose the same leading q^{1/3} term. It is a first-order statement only and does not
   decide integrality.

## Source list

### A. Non-bipartite random matching: LP relaxation versus integer optimum (the cleanest analogue)

**A1. M. Hessler and J. Wästlund, "Concentration of the cost of a random matching problem."** Preprint,
title page dated 20 August 2008; PDF at https://www.math.chalmers.se/~wastlund/martingale.pdf; OpenAlex
lists 2009 with no venue. Full text extracted and grepped; sections 1–3 and 6–8 read.
- *Proves.* Exponential(1) costs on K_n. Theorem 1.2: E|M_n − π²/12| < 1.706/√n. Theorem 1.3: for every x,
  P(|M_n − x| ≥ 1.114/√n) ≥ 1/8, so n^{-1/2} is the right order. Theorem 8.1: var(C_n) = O(1/n), where C_n
  is the cost of an LP-relaxed matching in Wästlund's "friendly model". That model is K_n with loops (loop
  cost mean 2), coefficients σ ∈ [0,1], and a minimum relaxed matching with coefficients in {0, 1/2, 1}
  (Lemma 6.1, from Wästlund 2010). E C_n = 1 − 1/4 + 1/9 − … ± 1/n² = π²/12 + O(1/n²), so
  E(M_n − C_n) = O(log n/n).
- *Bearing.* (i) Tool. The oracle process (reveal costs in order, use exchangeability of new edges and a
  correlation inequality) is the only rigorous variance method found for the value of a random
  non-bipartite matching LP. It does not control the probability of integrality. (ii) Warning. The LP/IP
  gap is far smaller than the fluctuations of either value, the same phenomenon as F3: LP* − T is about
  0.05–0.07, while sd(T) is presumably much larger. Concentration of LP* or T alone cannot give Lemma C,
  which supports the choice of mu*. (iii) The relaxation includes loops, so it is not exactly the
  fractional matching LP on K_n.

**A2. J. Wästlund, "The mean field traveling salesman and related problems,"** Acta Mathematica 204
(2010) 91–150 (citation as given in A1, A3 and A4; not opened). Source of the friendly model, the
exact E C_{k,n} formula, and the half-integrality and monotonicity lemmas used in A1. *Bearing.* Exact
formulas for an LP relaxation's expected value exist only in exponential/friendly models with a
recursive "add a vertex" structure. Nothing comparable exists for 3-index or pair-additive costs.

**A3. C. Lucibello, E. M. Malatesta, G. Parisi and G. Sicuro, "The random fractional matching
problem,"** J. Stat. Mech. (2018) 053301, doi:10.1088/1742-5468/aabbc8, arXiv:1802.02997. Full text
extracted; introduction, results table and conclusions read.
- *Shows (replica, non-rigorous, plus numerics).* The fractional matching LP on K_{2N} has half-integral
  optima made of disjoint odd cycles plus edges, and has the same asymptotic average optimal cost π²/12 as
  integer matching. The finite-size corrections differ. From Table 1, the uniform-weight (µ = 0) 1/(2N)
  coefficient is a ≈ −0.049 (fit) for integer matching and −0.3796 (theory) for the loop-free fractional
  problem. The anomalous N^{-3/2} correction of integer matching is absent in the fractional problem. The
  conclusions say "the contribution of cycles plays an important role" in that anomaly. Numerics average
  over at least 3·10⁶ instances per N.
- *Does not report* the fraction of instances whose fractional optimum is integral, or cycle counts.
- *Bearing.* This gives the closest analogue's LP/IP gap scale, E(T − LP) = Θ(1/N), on an edge-cost
  scale of 1/N. That points to O(1) odd cycles (my inference) and so to a non-degenerate limit for
  P(integral) on K_n. If correct, the 3-index answer (P -> 0, per the lane numerics) does not come from
  "odd cycles exist w.h.p." but from the mesoscopic structure seen in F2. Neither a scoop nor a tool for
  proof; it is cavity-method framing.

**A4. G. Parisi, G. Perrupato and G. Sicuro, "Random-link matching problems on random regular
graphs,"** J. Stat. Mech. (2020) 033301, doi:10.1088/1742-5468/ab7127, arXiv:1910.03504. Full text
extracted and grepped; section on cycle contributions read.
- *Shows (cavity, non-rigorous).* Finite-size corrections come from rare topological structures (cycles).
  Each odd cycle's cost shift φ_ℓ is positive; even cycles contribute little. The sum of odd-cycle
  contributions at z = 3 is about 1.065. The cutoff at cycle length ~√N explains the anomalous term. For
  K_N it re-interprets the Mézard–Parisi 1/N series as a cycle expansion.
- *Bearing.* Supplies the "single cycle embedded in an infinite tree" cavity computation. The 3-index
  analogue would be the cost shift of a fractional structure on k M-teams embedded in the local limit. F2
  says the relevant k grows with q, so the single-structure computation is not enough. Not a scoop.

**A5. E. M. Malatesta, G. Parisi and G. Sicuro, "Fluctuations in the random-link matching problem,"**
Phys. Rev. E 100 (2019) 032102, doi:10.1103/PhysRevE.100.032102 (abstract only via OpenAlex). Replica
derivation of the variance and large-deviation function of the optimal cost on K_N. *Bearing.* A possible
source of sd(T)-type benchmarks; nothing on the LP or the gap.

**A6. S. Sanghavi and D. Shah, "Tightness of LP via max-product belief propagation,"**
arXiv:cs/0508097 (v2, 12 April 2008). Full text extracted; Theorem 1.1 read.
- *Proves.* For MWIS with iid Exp(1) node weights on G(n, c/n), c ≤ 2e, or on random r-regular graphs,
  r ≤ 4: for a uniformly random node i, P[x*_i ∈ {0,1}] ≥ 1 − ε for large n. "Asymptotically tight"
  means integral on a 1 − o(1) fraction of nodes, not globally integral.
- *Bearing.* This is the methodological ceiling of local weak convergence and objective-method arguments:
  they decide integrality at a typical vertex, never global integrality. Our question needs the opposite,
  a global (or mesoscopic) statement. Note that F2 shows the 3-index LP optimum is fractional on most of
  M (k* about 0.8q), so local integrality fails here. If a local-weak-limit description of the LP optimum
  showed fractional mass on a positive fraction of M-teams, P(LP integral) -> 0 would follow at once.
  That is the most promising use of the objective method for Conjecture V (my inference). It is not
  attempted in any source.

**A7. M. Bayati, C. Borgs, J. Chayes and R. Zecchina, "On the exactness of the cavity method for
weighted b-matchings on arbitrary graphs and its relation to linear programs,"** J. Stat. Mech. (2008)
L06001, arXiv:0807.3159 (abstract). Also **S. Sanghavi, D. Malioutov and A. Willsky, "Linear programming
analysis of loopy belief propagation for weighted matching"** (NIPS 2007; arXiv:0705.0760 "Equivalence of
LP relaxation and max-product for weighted matching in general graphs"; search snippets only).
Deterministic: BP/cavity converges to the optimum iff the LP relaxation is tight (no fractional optimum).
*Bearing.* For the 3-index problem, "core nonempty" = "LP tight" = (by analogy) "max-product BP on the
3-index factor graph is exact". So P(BP exact) -> 0 is equivalent to Conjecture V. This is an equivalent
reformulation, not a tool. No random-instance results.

**A8. D. Gamarnik, T. Nowicki and G. Swirszcz, "Maximum weight independent sets and matchings in sparse
random graphs. Exact results using the local weak convergence method,"** Random Structures & Algorithms
(2006), doi:10.1002/rsa.20072, arXiv:math/0309441 (abstract only). Exact limits of maximum-weight
matching and independent set via a local optimality property (exponential weights, sparse graphs).
*Bearing.* As A6: first-order values only; no LP-integrality probability statement was found in it.

### B. Fractional versus integral factors in random (hyper)graphs: threshold level

**B1. S. Haber and M. Krivelevich, "On fractional K-factors of random graphs,"** RSA 30(4) (2007)
441–463, doi:10.1002/rsa.20144 (abstract via the TAU research portal). For strictly K1-balanced K, w.h.p.
a fractional K-factor appears exactly when every vertex is covered by a copy of K; for imbalanced K there
is a gap. It also gives expansion criteria for perfect fractional matchings in hypergraphs. *Bearing.*
Rigorous support that fractional and integral triangle factors share the threshold scale, which explains
F1 (LP* and T both lose the q^{1/3} term). Applying it to the weighted problem needs the level-set
reduction sketched in F1. Silent on the second-order question.

**B2. A. Johansson, J. Kahn and V. Vu, "Factors in random graphs,"** RSA (2008), doi:10.1002/rsa.20224
(found via OpenAlex; content from memory and from the lane's citation, not re-read). Integral
K3-factor threshold n^{-2/3}(log n)^{1/3}. *Bearing.* As B1.

**B3. M. Krivelevich, "Perfect fractional matchings in random hypergraphs,"** RSA 9(3) (1996) 317–334
per the DOI (OpenAlex title match, doi:10.1002/(SICI)1098-2418(199610)9:3<317::AID-RSA4>3.0.CO;2-#;
page range not verified). Also **P. Devlin and J. Kahn, "Perfect fractional matchings in k-out
hypergraphs,"** Electron. J. Combin. 24(3) (2017) P3.60, arXiv:1703.03513 (search snippet). Their
stopping-time result: a random r-graph process stopped when there are no isolated vertices has a perfect
fractional matching w.h.p. *Bearing.* Threshold-level only.

**B4. "On fractional triangle decompositions of random graphs,"** arXiv:2511.15877 (2025; title only via
OpenAlex, not read). Probably irrelevant (decompositions, not factors); listed so it is not re-checked.

### C. Random integer programs: LP–IP gap in general

**C1. S. Borst, D. Dadush and D. Mikulincer, "Integrality gaps for random integer programs via
discrepancy,"** SODA 2023, pp. 1692–1733, doi:10.1137/1.9781611977554.ch65, arXiv:2203.11863 (abstract
via OpenAlex). For max cᵀx, Ax ≤ b, x ∈ {0,1}ⁿ with m constraints and iid or logconcave random A, the
additive gap is O_m(log²n/n) w.h.p. (extending Dyer–Frieze 1989 and Borst et al. 2022). *Bearing.* Wrong
regime: a fixed number m of constraints with n → ∞ variables, where the gap → 0. Our LP has 3q equality
constraints on q³ variables, with the assignment structure. It shows that "LP integral w.h.p." is not
the generic behaviour even when the gap is small. Not a tool.

### D. Three-index polytope structure (for certifying gap 3 of THEOREM_VANISHING)

**D1. M. K. Kravtsov, V. M. Kravtsov and E. V. Lukshin, "On noninteger vertices of the polytope of the
three-index axial assignment problem,"** Diskretnaya Matematika 13(2) (2001) 120–143 (abstract via
mathnet.ru). Lower bounds and counting formulas for vertices with exactly r noninteger components. The
abstract gives no n = 3 count. *Bearing.* A possible independent check of the 1,386 / 1,080 counts in F2
(gap 3); full text not obtained. Companion: Kravtsov, "Description of the types of maximum noninteger
vertices…", Zh. Vychisl. Mat. Mat. Fiz. 46(10) (2006) 1908–1912 (covers n ≥ 4 only).

**D2. D. Vögler, "Kantorovich vs. Monge: a numerical classification of extremal multi-marginal mass
transports on finite state spaces,"** arXiv:1901.04568 (v3 2021; per search listing, journal version
"Geometry of Kantorovich polytopes and support of optimizers for repulsive multi-marginal optimal
transport on finite state spaces," J. Math. Anal. Appl. 502(1) (2021), not verified). Enumerates extreme
points of small multi-marginal (axial) Kantorovich polytopes and asks which are Monge (integral)
optimizers for pairwise costs. *Bearing.* The same object as LP integrality for pair-additive 3-index
costs, but symmetric marginals, worst case and small sizes; no random costs. A possible cross-check for
the vertex enumeration, and an example that non-Monge (fractional) optima for pairwise costs are generic
enough to be classified.

### E. Searched, nothing relevant (no need to repeat)

- OpenAlex full-text search (10 queries): "three-dimensional assignment random", "multidimensional
  assignment random decomposable", "random fractional matching", "LP relaxation integral random weights
  matching", "weighted triangle factor random", "random three-sided assignment game core", "integrality
  gap random assignment problem", "random multi-index assignment linear programming relaxation", "planted
  matching recovery", "fractional triangle factor random graph".
- OpenAlex forward citations of Frieze–Pegden–Tkocz (2 hits), Lucibello et al. (4 hits) and Frieze–Sorkin
  (8 hits): no LP-integrality or core result for random 3-index problems. One hit, Cheng–Liu–Tkocz–Xu,
  "Typical values of extremal-weight combinatorial structures with independent symmetric weights,"
  Electron. J. Combin. 30(1) (2023), doi:10.37236/10237, concerns first-order values on K_n with symmetric
  weights and is irrelevant here.
- Web search on the following returned only deterministic results:
  - LP tightness for random non-bipartite matching. The only hits were blossom-BP papers (arXiv:1509.06849)
    and Gelfand–Sanghavi–Willsky, which are deterministic or algorithmic.
  - Random hypergraph matching-game or multi-sided assignment-game cores from 2024 to 2026. The hits were
    the deterministic core ⇔ LP-integrality characterisation and complexity surveys. "Finite size scaling
    for the core of large random hypergraphs" (AAP 2008) is about the 2-core, not the game core.

## How this bears on Lemma C (summary)

- **Scoop risk: low.** Nothing proves or predicts P(core nonempty) -> 0 for the pair-additive model, or for
  any random axial 3-index model with growing q. Lucas 1995 (in GATE0) covers only q ≤ 3, or one side fixed
  (where P -> 1).
- **Closing risk from an analogue: none found.** The non-bipartite K_n matching LP is the natural analogue,
  but there is no theorem on its integrality probability. Physics suggests O(1) fractional odd cycles there,
  which would be the opposite regime from the 3-index numerics.
- **Usable tools.**
  - (a) The Hessler–Wästlund oracle process and correlation inequality, for variance bounds of LP-type
    values. It may adapt to Efron–Stein for mu*, but only with the exponential/friendly structure, which the
    pair-additive model lacks.
  - (b) The objective method / local weak convergence (A6, A8). If the local limit of (M, LP optimum)
    around a uniform M-team puts positive mass on fractional configurations, then
    P(LP integral) <= P(all q M-teams integral) -> 0 follows directly, without Lemma C. This needs a local
    limit at degree scale q^{1/3} → ∞, which does not exist in the literature for 3-index problems.
  - (c) Haber–Krivelevich and JKV, for the first-order (threshold) layer only.
- **Methodological warning (from A1 and A3).** LP and IP values can share the limit and concentrate at
  rate √n while the gap is O(log n/n). Any proof must work with a signed, gap-level quantity (mu*), not
  with concentration of LP* or T separately.
