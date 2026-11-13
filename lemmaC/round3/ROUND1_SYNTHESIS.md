# Round 1 synthesis and the round 2 plan

Written by the orchestrating session after seven of the eight round-1 lines and three
adversarial checks had reported. Statuses are in `CLAIMS.md`. "Verified" means a fresh
agent of the same model family tried to break the proof and failed; nothing has been read
by a different model family yet.

## What round 1 established

1. **The gap is asymptotically deterministic (verified, A1-A4).** Whp the optimal
   assignment M uses no pair of deficit above C·q^(-1/2+eps), and every team in the LP
   support has deficit below C·q^(-1/2+eps). With the Poincaré bound this gives
   Var T, Var LP*, Var D = O(q^(-1/6+eps)). So

   **P(core nonempty) -> 0 as soon as E D >= q^(-1/12+eps) for large q; in particular if
   liminf E D > 0.**

   This replaces the earlier requirement E D >> q^(1/6). File: `reports/pairdeficit.md`.

2. **The conjecture is not a statement about the limit tree alone (verified, T1-T2).** On
   the Poisson hypertree every invariant fractional perfect matching is the marginal of an
   invariant random integral perfect matching with the same cost. So the tree has one
   optimal cost a_tree for both problems. The difference between the LP and the integer
   optimum is that the finite market cannot realise the cheap invariant integral matchings.
   File: `reports/tree.md` §2.

3. **A first moment conditioned on local stars beats the tree value (verified for the
   Poissonised model, T5).** In the model with independent team costs (rescaled costs at an
   agent form a Poisson process of intensity xi^2/2), conditioning on every agent's binned
   star turns the near-tight hypergraph into a configuration model, and the expected number
   of assignments with usage profile rho has rate
   Sigma[rho] = 3H - 2∫h(rho)dLambda. The rate is concave (checker's proof, unreviewed) and
   vanishes at cost e_1* = 2.30248. With truncation at K = 5 the value is already 2.30246.
   So in that model a_IP >= 2.30. Files: `reports/tree.md` §5, `v2/`.

4. **Candidate LP constants on the tree (numerical).** Replica-symmetric integral fixed
   point: 2.2495 (two independent solvers agree: `../cavity/`, `redteam/`). Fractional
   fixed point: 2.197. Both are below 2.30.

5. **Finite sizes are pre-asymptotic.** At cost level K a near-tight team shares a pair
   with another with probability of order K^2·q^(-1/3), which is not small for q <= 10^4.
   The pair-additive constants measured at q <= 160 (LP 2.53 falling to 2.47; integer 2.57)
   are not the limits. The model without pair-sharing (`reports/sister.md`, model a3) has
   LP 2.12-2.19 and integer 2.18-2.24 at q <= 40, close to the tree numbers.

6. **Routes closed in round 1.** Counting price certificates (overshoots by exp(1.26·q)).
   Two-round exposure as a direct route (all its identities also hold for D ≡ 0). Rounding
   the LP optimum (it is fractional on 90% of agents). Local or block arguments (again).

## The architecture for a proof

Write xi = q^(2/3)·deficit for rescaled costs.

- **(I) Integer lower bound.** Whp 3q - T >= (beta - o(1))·q^(1/3) with beta = 2.30.
  - (I-a) Poissonised model: the star-conditioned first moment, written out with binning,
    truncation and a certified evaluation of the variational value.
  - (I-b) Transfer to the pair-additive market.
- **(II) LP upper bound.** E(3q - LP*) <= (alpha + o(1))·q^(1/3) with alpha < beta.
  - (II-a) An invariant solution on the tree with certified cost alpha < 2.30. Candidate:
    the replica-symmetric recursive tree process for X = min_k(xi_k - X_j - X_l), cost
    about 2.2495.
  - (II-b) Finite-depth approximation and local weak convergence: a nearly feasible
    fractional solution in the finite market with cost close to alpha·q^(1/3).
  - (II-c) Patching the slack at small cost.
- **(III) Conclusion.** E D >= (beta - alpha - o(1))·q^(1/3), far more than item 1 needs.

## Orchestrator's sketches (suggested routes; unverified; check every step)

**Sketch for (I-b): pointwise domination.** Let H be the hypergraph of teams with xi <= K,
with binned costs, in the pair-additive market, and H' the same object in the Poissonised
model. For a configuration h in which no two teams share a pair:
- the teams of h use disjoint pairs, so P(H contains h with those costs) is the product of
  the single-team probabilities, exactly as for H';
- "no other team is near-tight" is at most the probability that no team pair-disjoint from
  h is near-tight, which does not depend on h's values and is at most exp(-mu + Delta) by a
  Janson-type inequality for up-sets (Riordan-Warnke), with mu the Poissonised mean up to
  O(K^3) and Delta = O(K^5·q^(2/3)) = o(q).
So P(H = h) <= exp(o(q))·P(H' = h). For h with s pair-sharings each sharing costs at most a
factor C(K, bins)·q^(2/3), which is exp(o(q)) in total when s <= q^(2/3+eta); and
P(s > q^(2/3+eta)) -> 0 by Markov. Then
E[#cheap assignments; stars typical, s small] <= exp(o(q))·(Poissonised value), and the
events "stars atypical" and "s large" have vanishing probability. Teams above K are handled
by the plain union bound as in the checker's truncation fix; disjoint teams have
independent costs in both models.

**Sketch for (II-b), (II-c): Aldous's split.** Write each pair deficit as d = min(d', d'')
with independent d', d'', P(d' > t) = (1-t)^(1-delta), P(d'' > t) = (1-t)^delta.
- Build the main solution from the d'-market by a finite-depth local rule: for an invariant
  solution m on the tree, y^v_e = P(e is v's team | r-ball of v) sums to 1 over v's teams;
  x_e = min over e's three agents of y^v_e is a packing; by martingale convergence the
  slack s_v = 1 - load_v has mean eps_r -> 0 as r -> infinity. Local weak convergence in
  probability carries cost and slack to the finite market. The d'-market is a U[0,1] market
  with deficits scaled by about 1/(1-delta), and d <= d'.
- Patch with d'', which is independent of the slack vector s. Write s as the mean of the
  indicator of a random balanced triple of agent sets R (dependent rounding, sizes floor or
  ceiling of the sector sum S), chosen independently of d''. On each R take the optimal
  assignment of the d''-sub-market. Since (1-t)^delta <= 1 - delta·t, d'' is
  stochastically below U[0,1]/delta, so by the Frieze-Pegden-Tkocz upper bound the expected
  cost is at most (c_2/delta)·|R|^(1/3). The average over R is a fractional solution with
  marginals exactly s and expected cost at most (c_2/delta)·(eps_r·q)^(1/3).
- Then limsup E(3q - LP*)/q^(1/3) <= alpha/(1-delta) + c_2·eps_r^(1/3)/delta; let r -> infinity,
  then delta -> 0.

**Remark for (II-a).** The integral recursion maps a shift c of its input law to a shift
-2c of its output (T4, R2). So if T(G) is G shifted by c, then G shifted by c/3 is a fixed
point. Existence can therefore be sought for the map "T followed by re-centring", where the
unstable direction is quotiented out.

## Round 2 lines

| Line | Task |
|---|---|
| localconv (running) | Local weak convergence of the pair-additive near-tight hypergraph, with error bounds |
| sc | (I-b): transfer of the star-conditioned first moment to the pair-additive market |
| firstmoment | (I-a): full write-up for the Poissonised model with a certified constant |
| lplocal | (II-b), (II-c): the LP locality upper bound |
| rde | (II-a): existence and certified cost of the tree solution |
| pscout | Large-q numerics in the Poissonised model as a falsification test of alpha < beta |
| redteam2 | Attack the architecture |
