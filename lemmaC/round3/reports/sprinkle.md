# Round 1 report: sprinkle line (summary)

Saved by the orchestrator from the agent's message. Scripts and data: `../sprinkle/`
(regret.py, balance.py, repair.py, drot.py, brute.py and their outputs).

Verdict: none of the three targets (E D -> infinity; E D >= c; limsup P(core nonempty) < 1)
is proved. The natural two-round schemes are exact affine maps of the instance, so alone
they give identities that also hold for D ≡ 0 and cannot give E D > 0.

Notation: pairs i (3q^2 of them), d_i = 1 - value. I^{-i} is I with pair i's value set
to 0. rho_IP(i) = T(I) - T(I^{-i}), rho_LP(i) = LP*(I) - LP*(I^{-i}).

Proved (claimed complete, unverified, not on the critical path):
- S0 (affine covariance): deficits alpha·d + beta with alpha > 0 leave the optimal integer
  and fractional solutions unchanged and scale D by alpha.
- S1 (one pair exposed last): with all other values fixed and t = d_i,
  T(t) = T_inf + (tau_i - t)^+, so i is in M iff d_i < tau_i; LP*(t) is convex piecewise
  linear with slope -load_i(t). Consequences: (a) sum_i E rho_IP(i) = E(3q - T) and
  sum_i E rho_LP(i) = E(3q - LP*), so E D = sum_i E[D(I^{-i}) - D]; (c) a.s. D > 0 iff some
  pair has d_i in its window (tau_i, sigma_i); (d) given the others and i in M, d_i is
  uniform on (0, tau_i ∧ 1) and rho_IP(i) = (tau_i ∧ 1) - d_i, so the regret and the deficit
  of a pair of M are exchangeable.
- S2 (deletion is dilation): replacing each pair independently with probability eps by a
  fresh deficit uniform on [1, 1/(1-eps)] gives iid U[0, 1/(1-eps)] deficits, so
  P(D(I^R) > 0) = P(D(I) > 0) exactly.
- S3 (deletion flow balance): E[N_kill; D > 0] = E[N_create; D = 0], where N_kill counts
  pairs whose deletion closes the gap and N_create pairs of M whose deletion opens one.
- S4: P(D > 0) >= P(D(I^{-i}) > 0)/5 for i uniform on M's pairs; E D <= E D(I^{-i}) <= 4·E D.

Gaps stated by the agent: G1, P(D(I^{-i}) > 0) >= c (would give a uniform lower bound);
G2, E[N_kill | D > 0] = o(E[N_create | D = 0]) (would give the theorem); G3, a bound
sum over i not in M of rho_LP(i) <= K·D with K = O(1).

Numerics (evidence): the identities check at q = 3..20. Repairing one deleted pair changes
7.6 of 16 teams at q = 16, so the repair is not local. The LP restricted to teams sharing a
pair with M captures 2-4% of D: the gap lives on teams whose three agents come from three
different teams of M.
