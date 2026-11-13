# Literature check 1 (random three-dimensional assignment, LP versus integer optimum)

Returned to the orchestrator by a lookup agent whose parent could not be reached.
Legend: [O] full text opened by that agent; [S] search snippet only; [I] its inference.
Nothing here was re-opened by the orchestrator.

## Martin–Mézard–Rivoire [O: arXiv cond-mat/0407623, cond-mat/0507180]
Costs are iid per hyperedge, not sums of pair costs.
- Replica-symmetric equations are "inconsistent when d > 2". RS entropy is negative beyond
  beta_s = 0.412 ± 0.001.
- d = 3, tripartite normalisation: frozen one-step value 3.126 ± 0.002. Branch-and-bound at
  M ≤ 22: 3.06 ± 0.03 (quadratic fit), 3.09 ± 0.03 (power-law fit). No RS zero-temperature
  value is given for d = 3; read off a figure it lies roughly in [2.72, 2.9] [I].
- No LP relaxation, fractional matching or belief propagation is discussed.
- Annealed bound C(d) ≥ e^(d−2), i.e. 2.718 for d = 3, against 3.126.

## LP versus integer optimum
No source found that states agreement or difference (two searches).
- Frieze–Pegden–Tkocz [O: arXiv 1901.07167]: the pair-additive model. Optimum
  Theta(n^(1/3)); greedy is (3+o(1))-approximate; "the real aim is to replace 3 by 1". No LP.
- Frieze–Sorkin [O: arXiv 1008.0390]: iid Exp(1) costs. Optimum between Omega(1/n) and
  O(log n/n) (the upper bound via Johansson–Kahn–Vu, non-constructive); an algorithm gives
  O(n^(−1+o(1))). "We do not even know how the cost scales with n." No LP comparison.

## Thresholds
- Johansson–Kahn–Vu [quoted in Kahn, O]: M > C_r n log n, order only.
- Kahn, arXiv 1909.06834 [O]: M > (1+eps)(n/r) log n, sharp constant.
- Kahn, arXiv 2008.01605 [O, first page]: a perfect matching appears exactly when the edges
  first cover all vertices.
- Krivelevich 1996 [not opened; restated as Devlin–Kahn Theorem 4, O]: at the first time
  without isolated vertices the hypergraph has a perfect fractional matching whp.
- Devlin–Kahn [O: arXiv 1703.03513]: k-out model, fractional only, including r-partite.
  Integer needs k "at least exponential in r"; for fractional "we don't even know that they
  must be larger than 2".
So in the non-partite random r-graph both hitting times equal the disappearance of isolated
vertices. At bounded degree the two can separate [I].

## Matching on K_n
Hessler–Wästlund [O: math.chalmers.se/~wastlund/martingale.pdf]: an LP relaxation with
loops has expected cost pi^2/12 + O(1/n^2), "hence the expected difference between M_n and
C_n is O(log n/n)". The plain LP is sandwiched between them [I].
