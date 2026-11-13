# Round 1 report: dualcount line (summary)

Saved by the orchestrator from the agent's summary message. Scripts and data: `../dualcount/`.

Verdict: the dual first moment does not give P(core nonempty) -> 0. The best version
overshoots by about exp(1.26·q). Line stopped.

Proved (claimed complete):
- Relaxed first moment: if p certifies the identity, every p - eta with eta >= 0 and
  per-team sum <= tau' keeps all off-diagonal constraints, with diagonal sums in
  [d_i - tau', d_i]. So P(core nonempty) <= (q!)^2·(6/tau'^3)^q·∫ prod_i w_tau'(sigma_i)·
  P(no fully off-diagonal team blocks at p) dp.
- The relaxed event implies D <= tau·q^(1/3) (tau' = tau·q^(-2/3)).
- In rescaled prices mu(p) = (1/6q^2)·sum (s_i + t_j + u_k)_+^3 + O(q^(2/3)), and Harris
  gives P(no block) >= exp(-mu(1+o(1))).
- A nonempty core is a.s. a bounded polytope of dimension 2q-2 modulo shifts, so it has at
  least 2q-1 vertices.
- One-sector elimination (X2 in CLAIMS.md).

Numerics (evidence only):
- Annealed rate r(tau) of the bound: 4.81, 3.57, 2.16, 1.45, 1.26, 1.45, 1.70 at
  tau = 0.15, 0.3, 0.75, 1.5, 2.7, 4.5, 6. Positive at every admissible constraint.
- Exact vertex enumeration of nonempty cores: mean log count 3.34, 5.04, 6.70, 7.92, 9.59
  for q = 3..7, about exp(1.56·q). Coordinate widths 1.3 to 2.0 times q^(-2/3).
- Why it fails: the bound is P(nonempty)·E[multiplicity]. It needs cores thinner than 0.04
  per dimension; measured effective widths are 0.15 to 0.5, and P(nonempty) decays only
  like exp(-0.13·q).

Gaps named by the agent: the Janson-type upper bound exp(-mu + Delta) for general monotone
events is from memory (Riordan-Warnke), not opened; Laplace asymptotics not written out; a
canonical-vertex count is not evaluated; a truncated first moment is untested.
