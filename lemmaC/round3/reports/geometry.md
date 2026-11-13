# Round 1 report: geometry line (summary)

Saved by the orchestrator from the agent's message. Scripts, data and tables:
`../geometry/` (tables.txt is the output of summ.py).

Terms: M = identity, e_t = (a_t, b_t, c_t). An off-M team (a_i, b_j, c_k) is pair-keeping
if exactly two of i, j, k coincide (it keeps one pair of M), and pair-disjoint if i, j, k
are distinct.

Claimed complete (unverified, off the critical path):
- P1 (cone): feasible directions at M are y - sum_t z_t e_t with y >= 0 on off-M teams and
  equal loads z_t on a_t, b_t, c_t. Edges are the support-minimal such y, with at most
  2k-1 teams when k M-teams are touched.
- P2 (blow-up): with T(z) the best integer value when the agents of e_t carry load z_t,
  D >= [T(z) - sum z_t v(e_t)]/max z, and D = 0 iff T(z) = sum z_t v(e_t) for every z.
- P3 (integral neighbours): M' = (sigma, tau) is adjacent to M iff <sigma, tau> is
  transitive on the set of moved indices; counts per k-set 3, 26, 426.
- P4 (gain): v·d = G_A(U) + G_B(V) + G_C(W), where U, V, W are the off-diagonal pair loads,
  each a circulation on the M-teams. At q = 2 the fractional gain is half the sum of three
  integral gains, so D = 0 always.
- P5 (pair-keeping family): a potential/shortest-path characterisation of when no
  pair-keeping direction improves.

Enumeration of extreme rays at q = 2, 3, 4: 5, 422, 356,000 (integral 3, 35, 548). q <= 3
exact; q = 4 complete up to float sign tests.

Numerics (floating point):
- Restricted LP gaps: with pair-keeping teams only, P(gap > 0) is 1-4%; with pair-disjoint
  teams only it is about 0; the full LP has 21% (q = 8) to 100% (q = 48). The gap at these
  sizes needs both kinds together.
- Thinning the pair-keeping teams to 1/2, 1/4, 1/8 cuts the mean gap from 0.105 to 0.013,
  0.003, 0.0003 at q = 32.
- The double cover z = 2 improves in about 90% of instances with D > 0 and captures 82-91%
  of E D.
- Control with independent team values of the same marginal law: P(D > 0) = 0.75 at q = 8
  and 0.99 at q = 16; E D = 0.134 and 0.217. This agrees with the sister line's numbers.

Orchestrator's reading: at accessible sizes the gap in the pair-additive market runs
through teams that keep a pair of M. In the limit those teams are a vanishing fraction of
the near-tight ones, so the finite-size mechanism is not the limiting one. This is one more
sign that q <= 128 is pre-asymptotic.
