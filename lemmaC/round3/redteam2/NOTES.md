# redteam2 notes (round 3): attack on the two-constants architecture

Floating point is evidence only. Scripts here: `lp_big.py`, `endo.py`, `endo2.py`, `ip_iid.py`.

## Verdict
No contradiction found between (I) and (II). The architecture is coherent; the risk is in certification.

## 1. Where (II) uses fractionality, and why it does not produce an integral assignment
- The only step is x = E[m | tree, costs] and its finite-depth version y^v_e = E[m_e | r-ball of v].
  Cost is c(m) by the tower property (xi_e is ball-measurable); sum over v's teams is 1; Levy upward gives slack -> 0.
  x cannot be degenerate: sum_e x_e = 1 a.s. forces dependence on the tree (rho_RS(xi_e) alone does not sum to 1).
- The same argument gives an INTEGRAL near-perfect matching (keep e iff y^u_e > 1/2 at all three agents) iff the
  recursive tree process is endogenous. So (I) + (II-a) imply non-endogeny. Theorem K does not prove it: the bipartite
  logistic RDE has a non-decaying shift multiplier (-1) and is endogenous (Bandyopadhyay).
- Test (`endo2.py`, K = 4 truncated tree, 2e5 pairs, 40 generations, both marginals re-centred): two copies on the same
  tree started independent settle at corr(X, X') = 0.51, P(same root team) = 0.407, P(|X - X'| < 0.01) = 0.008, cost 2.277.
  No drift to the diagonal. `endo.py` (explicit trees to depth 5, iid boundary) agrees: P(same team) 0.21, 0.30, 0.35,
  0.39, 0.40. Heuristic reason: on a common argmin team Z_parent = -(Z_j + Z_l), variance doubles; the diagonal is
  linearly unstable. So x is genuinely fractional (sum_e x_e^2 about 0.41) and no integral rule comes out.

## 2. Severity-ordered weak points
1. (II-a) must be an EXACT existence theorem with an enclosure of alpha. A finite-depth numerical rule cannot replace
   it: with Aldous patching the loss is min over delta of alpha*delta + c_2*eps^(1/3)/delta = 2 sqrt(c_2 alpha) eps^(1/6)
   = 6.6 eps^(1/6); to fit in beta - alpha = 0.053 needs slack eps < 2e-13. alpha = E(X1+X2+X3)_+^4/8 moves by 9c under a
   shift c of G, so G must be located to about 0.005 in the shift direction. No narrow order interval is invariant
   (T antitone and T(G+c) = T(G) - 2c leaves [G-c, G+c]); Schauder on the re-centred map gives existence but no enclosure.
   A Newton-Kantorovich computer-assisted proof is needed (DT is an integral operator with smooth kernel; I - DT has
   eigenvalue 3 in the shift direction). Useful pin: any exact fixed point has E(X1+X2+X3)_+^3 = 6.
2. Margin: beta - alpha = 0.053 (2.3%), both numbers floating point; concavity of Sigma unreviewed. Physics value for
   a_IP is max_beta f_RS = 2.32 +- 0.01, only 0.02 above e_1*: the required order e_1* <= 2.32 holds at about 2 sigma.
   In the iid model IP/thr = 1.0735, 1.0715, 1.083 at q = 8, 12, 16 against e_1*/2.1496 = 1.0711: consistent, thin.
3. (I-b): the sketch's Markov step covers sharings inside H only. Sharings between M's teams above K and H depend on M.
   Fix: count pairs (t in H_L, t' in H_K) sharing a pair, O(L^2 K^3 q^(2/3)) by Markov, and bound teams of M above L by
   a*q/L on the event being bounded; L = log^2 q gives o(q/log q) sharings, enough. Typical stars must be defined in the
   Poissonised sense. Janson for up-sets is used in the right direction.
4. (II-c): needs E cost <= c_2 n^(1/3) for all n in expectation (or failure probability O(n^(-2/3))); a whp-only form of
   Frieze-Pegden-Tkocz leaves 3*eps*q*P(fail), which is not o(q^(1/3)).
5. Order of limits is fine: q, then r and K, then delta. (III) needs only E (3q-T) from the whp bound.

## 3. Finite-size data: pre-asymptotic claim tested
`lp_big.py`: LP restricted to teams with xi <= K (upper bound on LP deficit; K = 7 vs 6 differs by 1e-4 at q = 1000).
See out_pair_*.jsonl, out_iid_*.jsonl. (3q - LP*)/q^(1/3):
- pair-additive: 2.47 (q=160, round 1), 2.39 (500, n=1), 2.352 +- 0.006 (1000, n=4; a fifth seed 2.325),
  2.313 +- 0.011 (2000, n=3). (value - 2.197)/q^(-1/3) = 1.1 to 1.5 from q = 16 to 2000: compatible with one limit
  near 2.2 and a correction about 1.4 q^(-1/3). A straight line through q <= 160 alone had predicted 2.40 at q = 1000.
- iid teams: 2.214 +- 0.013 (250, n=6), 2.184 +- 0.012 (1000, n=4), 2.202 (2000, n=2), 2.196 (4000, n=1). Flat at the
  fractional value 2.197, below 2.2495 as (II) requires.
Per-instance sd of the ratio at q = 1000 is 0.012 (total sd 0.12, down from 0.45 at q <= 64), as A4 predicts.
- Integer optimum, iid teams, q = 32, K = 5.5 (`ip_iid.py`, out_ip_iid.jsonl; 5 of 6 solved to optimality, the sixth
  bracketed 2.288 to 2.313): IP 2.32 +- 0.03 against LP 2.21 on the same instances, gap 0.11 q^(1/3). The finite-q
  analogue of the star bound, 2.3025 * thr(32)/2.1496 = 2.25, is respected. q = 64 did not solve in 40 s.
Not tested: the integer optimum at large q in either model; pair-additive q = 4000 was started and stopped at the time box.
