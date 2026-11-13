# Claims register, round 3

Kept by the orchestrating session. Status values: conjectured, claimed, verified, refuted.
"Verified" means the statement survived a fresh agent that saw only the statement and the
proof and was told to break it. A full proof also needs a Codex read.

Agent runs used: 38 of 40 (round 1: 8 lines; checks v1, v2, v3; round 2: localconv, sc, firstmoment, lplocal, rde, pscout, redteam2; checks v4, v5, v6, v7; round 3: interface, rdeA, rdeB; check v8; round 4: cert; checks v9, v10; write-outs fixA, fixB, fixC; second reads v11-v16). Codex reads: 4. Corrections owed: final/ERRATA.md.

## Carried in from earlier rounds

| ID | Statement | Status | File | Checked by |
|---|---|---|---|---|
| P1 | Var T, Var LP* ≤ 4·E(3q − T) | claimed | concentration/CONCENTRATION.md (A2) | same-family check, check/CHECK.md |
| P2 | Proposition R: P(core nonempty) ≤ 8·E(3q − T)/(E D)² | claimed | SYNTHESIS_2026-11-11.md | same-family check, check/CHECK.md |
| P3 | E(3q − T) = Θ(q^(1/3)) (Frieze–Pegden–Tkocz Thm 2) | claimed (published source read) | check/CHECK.md | same-family read of arXiv source |
| P4 | 3q − T ≥ 2.15·q^(1/3) whp | claimed | proof/PROOF_ATTEMPT.md (Thm 1) | none |
| P5 | E(3q − LP*) ≥ 1.36·q^(1/3) | claimed | proof/PROOF_ATTEMPT.md (Thm 2) | none |
| P6 | D > 0 forces an odd Berge cycle in the LP support | claimed | proof/PROOF_ATTEMPT.md (Lemma 5) | none |
| P7 | Ordinal limit: LP optimum always integral | claimed | dual/DUAL.md (Theorem O) | simulation only |
| P8 | Harris: conditioning on M lowers LP* | claimed | proof/PROOF_ATTEMPT.md (Lemma 4) | none |
| P9 | Localisation: D ≥ sum of block gaps | claimed | ../THEOREM_VANISHING.md (Lemma 1) | coordinator read |
| P10 | P(T = S) = 1/q! | claimed | dual/DUAL.md (D2) | coordinator read |
| C0 | Two-constants conjecture: E D ≥ c·q^(1/3) | conjectured | SYNTHESIS_2026-11-11.md | cavity (non-rigorous), data to q = 64 |

## Round 1 lines

| Line | Formulation | Sees favoured approach? |
|---|---|---|
| tree | Statement on the limiting Poisson hypertree | yes |
| pairdeficit | Whp no used pair has deficit above q^(−2/3+ε) | lemma statement only |
| dualcount | First moment over dual certificates | no |
| sprinkle | Two-round exposure | no |
| sister | Easier sister models and transfer | no |
| geometry | Edges of the assignment polytope at an integral vertex | no |
| scout | Targeted numerics | no |
| redteam | Argue E D = o(q^(1/3)) | no |

## New claims

| ID | Statement | Status | File | Checked by |
|---|---|---|---|---|
| S1 | Theorem S: random r-regular 3-partite 3-uniform hypergraph, r <= 6: fractional perfect matching exists, integral one does not; gap >= (1 - mu_r)n | claimed | reports/sister.md §1 | not yet |
| S2 | Proposition A: iid team costs, n·IP >= e - o(1) whp | claimed | reports/sister.md §1 | not yet |
| T1 | Lemma B0: finite hypertree has at most one perfect matching | verified | reports/tree.md §2 | v1 (fresh adversarial agent): holds |
| T2 | Theorem B: every invariant fractional perfect matching on the tree lifts to an invariant random integral one with the same cost | verified | reports/tree.md §2 | v1: holds; exact check on 700 finite pieces, 7,134 roots; add the convention for x_f = 1 |
| T3 | Lemma D: weak duality on the tree | verified with fix | reports/tree.md §3 | v1: needs (x,u) jointly unimodular (take u a factor, or couple conditionally independently); E|u_o| finite suffices |
| T4 | Theorem K: perturbations of the cavity recursion grow by exactly 2 per generation | verified with fix | reports/tree.md §4 | v1: (ii), (iii) hold; (i) 'only if' fails for n >= 3 when the top cells tie (counterexample), add g^{k*}_{n-1} < g^{k*}_n |
| T5 | Star-conditioned first moment: a_IP >= 2.3025 in the Poissonised sister model | verified (formula, concavity, value); write-up and certified evaluation owed | reports/tree.md §5; v2/ | v2 (fresh adversarial agent): re-derived the formula, exact enumeration agrees, e_1* = 2.30248 on three grids, 3,000 profiles fail to beat it; proved Sigma concave (unreviewed); truncation needs a 'high option' fix, K = 5 gives 2.30246 |
| T6 | (S-c) transfer of T5 to pair-additive costs | conjectured | reports/tree.md §5 | none |
| N1 | Sister model without pair-sharing: LP deficit 2.12-2.19, E gap 0.12-0.23 at q <= 16 (numerics) | evidence | reports/sister.md §2 | none |

| R1 | Exchange count: expected number of ways to swap s teams of a fixed assignment for s teams of deficit <= c·q^(-2/3) is at most (s·c³/6q)^s | claimed | redteam/REPORT.md | not yet |
| R2 | Integral tree recursion maps a uniform shift c to -2c (independent of T4(iii); two lines agree) | verified | redteam/REPORT.md; reports/tree.md §4 | v1: K(iii) holds |
| N2 | Tree fixed point (integral, replica-symmetric) costs 2.2496; Bethe entropy negative at beta = 3, 4; free energy rises to about 2.32 at beta = 2 (Monte Carlo) | evidence | redteam/REPORT.md | agrees with cavity 2.2495 and with T5's 2.3025 in order of magnitude |
| N3 | LP optimum is fractional on about 90% of agents for q >= 48, so rounding needs order-q repairs | evidence | redteam/REPORT.md | none |
| N4 | Pair-additive LP deficit ratio drifts down slowly: about 2.53 (q <= 48), 2.50 (64-128), about 2.46 (160, n = 8) | evidence | redteam/lpfrac.jsonl | none |
| A1 | Lemma U: uniform expansion of s-close neighbourhoods over all sets of disjoint positions | verified with cosmetic fix | reports/pairdeficit.md §2 | v3: holds; write 0.0208 for 0.0209 in (b) |
| A2 | Theorem A: whp M uses no pair of deficit >= C·q^(-1/2+eps) | verified | reports/pairdeficit.md §4 | v3: holds; brute force at q = 3, 4, 5 (127,200 rotations, 0 violations); orchestrator also read the argument |
| A3 | Theorem B(LP): whp every LP-support team has deficit < C·q^(-1/2+eps) | verified | reports/pairdeficit.md §3 | v3: holds; 391,852 price checks at q = 8..16, 0 violations |
| A4 | Corollary: Var T, Var LP*, Var D = O(q^(-1/6+eps)); hence P(core nonempty) -> 0 if E D >= q^(-1/12+eps), e.g. if liminf E D > 0 | verified, conditional on P3 (FPT Theorem 2) | reports/pairdeficit.md §1, §5 | v3: holds |
| X1 | Dual first moment (count of price certificates) bounds P(core nonempty) | refuted as a route | dualcount/ (agent's message; summary in reports/dualcount.md) | agent's own numerics: best version overshoots by about exp(1.26·q); a nonempty core has about exp(1.56·q) vertices while P(nonempty) decays only like exp(-0.13·q) |
| X2 | One-sector elimination: a certificate exists iff some a in R^q makes each diagonal team the maximiser of v(i,k,k) - a_i and the identity a maximum matching of w_a(j,k) = max_i(v(i,j,k) - a_i) | claimed | dualcount/ | not checked; not on the critical path |
| K1 | Sprinkle identities S0-S4 (affine covariance; one-pair exposure; deletion flow balance; P(D>0) >= P(D(I^{-i})>0)/5) | claimed | reports/sprinkle.md | not checked; off the critical path |
| X3 | Two-round exposure alone gives E D > 0 | refuted as a route | reports/sprinkle.md | the schemes are affine maps; every identity also holds for D ≡ 0 |
| G1 | Geometry propositions P1-P5 (tangent cone at M, blow-up criterion, integral neighbours, gain as three circulations, pair-keeping family) | claimed | reports/geometry.md | not checked; off the critical path |
| N5 | At q <= 48 the gap needs pair-keeping and pair-disjoint teams together; pair-keeping teams alone give a gap in 1-4% of instances | evidence | reports/geometry.md | none |
| L1 | Local convergence, one root: total-variation distance between the r-ball in H_q(K) and in the Poisson hypertree is at most eps_1, dominated by (9/10)·K²·q^(-1/3) times the expected number of teams in the ball | verified | localconv/LOCALCONV.md §1 | v4: holds; simulated at q up to 262,144, constant near-sharp |
| L2 | Two roots: joint law within 4·eps_1 of the product; variance of neighbourhood-type counts O(q^(5/3)) | verified with wording fix | localconv/LOCALCONV.md §2 | v4: holds; the fluctuation rate is O_P(q^(-1/6)), only the bias is O(q^(-1/3)) |
| L3 | Uniform integrability: E[(1/q)·sum over M of (xi_e - K)_+] <= C·K^(-1/6) for K >= 64, all q, and the same for the LP optimum; repair corollary | verified with fixes, conditional on P3 | localconv/LOCALCONV.md §3 | v4: holds; C_0 must be sup_n E(3n-T_n)/n^(1/3) (finite, since the ratio is bounded for each n and asymptotically; orchestrator's remark); LP constant differs; repair corollary needs the sub-assignment to be independent of the second half of the split; the explicit bound is vacuous at simulable sizes |
| N6 | Pair-sharing partners per near-tight team: (9/10)·K²·q^(-1/3), matched by simulation at q = 32..128; at K = 6 it stays >= 1 until q is about 3·10^4 | evidence + exact first moment | localconv/out_pairshare.txt | none |
| N7 | Scout: largest used pair deficit × q^(2/3) grows like log q (3.1 to 3.8 for M at q = 16..48); P(core nonempty) decays about like exp(-0.12·q); shortest odd cycle in the LP support has length 3 | evidence | scout/REPORT.md | none |
| F1 | LP locality upper bound: for any invariant random integral perfect matching m on the tree, limsup E(3q - LP*)/q^(1/3) <= c(m), given local convergence (L1) and the FPT bound (P3) | verified with a trivial fix | lplocal/LP_LOCALITY.md | v5: holds (93%); one l1/l2 tolerance fix; fractionality enters only through Levy's theorem (E[m_e | r-ball] -> x_e, not m_e) and the patch average, so no integral assignment results unless m is a factor of the tree; residual risk is the unwritten measurable framework and the exact form of the FPT bound (two agents who opened the arXiv source read it as U[0,1] entries, in expectation) |
| E1 | Tree recursion: G(x) = E[(x+X'+X'')_+^3]/6; a solution gives an invariant perfect matching with cost c = E S_+^4/8 = 9·E X (S = sum of three iid copies) | verified (formulae); process construction to be written out | rde/RDE.md | v9: re-derived, E S_+^3 = 6 and 9·E X = E S_+^4/8 confirmed two ways; Codex: formulae check out, the two-directional message process on the infinite-degree tree needs a cylinder-consistency argument |
| E2 | Proposition E: if a convex compact set K of centred laws is invariant under "T then re-centre", a fixed point exists with c <= 9·sup_K(E X + kappa/3) | claimed | rde/RDE.md | not yet |
| E3 | EXISTENCE of the fixed point | verified (E6-E8); a certified bound c_n < 2.302 for some n is the ONLY remaining step (E9) | rdeB/ROUTES.md | v9; independent values c_1..c_6 = 3.41528, 2.66022, 2.39809, 2.30170, 2.26786, 2.25596 |
| N8 | Fixed point numerically: c = 2.2495456100, E X = 0.2499495122 (not 9/4); modulo shift the map is a local contraction with factor about 0.35 (spectrum of P: 1, -0.175, 0.042, ...) | evidence | rde/ (rde_solve.py, rde_spectrum.py) | agrees with two earlier solvers |
| I1 | Theorem I (iid team deficits, each a sum of three uniforms): whp every perfect assignment has total deficit >= 2.302·q^(1/3); interval-arithmetic certificate gives 2.302245 at K = 5, 50 bins | verified | firstmoment/THEOREM_I.md, certify.py | v6 (93%): every lemma re-derived by hand, no floating-point leak in the certificate, INDEPENDENT interval enclosure gives beta >= 2.302028 (so the exponent at 2.30 is <= -0.00339), exact enumeration of 77,527 small cases shows 0 violations; cosmetic: define 'simple' to exclude repeats across bins |
| I2 | Theorem T (transfer): E[Z_eps; H in G, few pair-sharings] <= exp(o(q))·(E'[Z_eps; H' in G^(2theta)] + exp(-(theta/4)·q·ln q)·E'[Z_eps]), with lemmas L1 (factor per pair-sharing), L2 (Janson-type void bound, proved directly), P5 (typicality in the pair-additive market) | verified | sc/TRANSFER.md | v7: holds (93%); L2 is the Boppana-Spencer proof (union bound, exact independence, Harris once in the right direction); numerics at q = 2..32 show no violation; constraints on the interface: theta small, G built from bounded-increment statistics, all of Z_eps covered; purely asymptotic |
| I3 | Theorem I-pa: in the pair-additive market, for every beta <= 2.30224, P(some perfect assignment has total deficit <= beta·q^(1/3)) -> 0. Interface (I-a*) with explicit constants; new Lemma Q replaces the compactness step of Theorem I | verified with cosmetic fixes | interface/THEOREM_I_PA.md, constants.py | v8 (95%): all seven points hold; exponent re-evaluated in exact rationals: <= -0.00059298 at beta = 2.30; fixes: use the true ell_max = 8.2718 in the condition on A; the 'size of q' remark understates how large q must be (ln q > 1.6e5), asymptotic only |
| RT2 | Second red team: no contradiction between the integer bound (I) and the LP bound (II); fractionality is used only in x = E[m | tree, costs]; the same argument would give an integral assignment only if the tree process were endogenous, and a two-copy test shows it is not (correlation 0.51) | evidence + argument | redteam2/NOTES.md | none |
| N9 | Pair-additive LP constant keeps falling: 2.39 (q = 500), 2.352 (1000), 2.313 (2000); iid-team model flat at 2.19-2.20 up to q = 4000; iid integer optimum at q = 32 is 2.32 ± 0.03 | evidence | redteam2/NOTES.md | supports the pre-asymptotic reading and alpha < beta |
| N10 | Poissonised model, K = 6: LP cost per triple 2.165 (q = 16) rising to 2.200 (q = 512), extrapolating to 2.20-2.21; exact integer optimum 2.253, 2.277, 2.285, 2.311 at q = 16, 24, 32, 40 (all closed to optimality), fits give a limit of 2.32 or more; gap 0.09-0.115 and rising | evidence | pscout/results.txt | falsification test passed: consistent with alpha <= 2.25 and beta >= 2.30 |
| E4 | rdeA Lemma 1: exact nonlinear identity in quantile coordinates, T(Q)(u) - T(Q_0)(u) = -2∫p~(u,v)(Q - Q_0)(v)dv - (tail term), p~ a probability kernel; Lemma 2: cell-wise envelope of p~ for a piecewise-constant reference density | claimed | rdeA/EXISTENCE.md | not yet |
| E5 | rdeA Theorem A (reduction): hypotheses H1 (defect), H2 (invariance), H3 (envelope room), H5 (tail envelopes) imply a fixed point with an explicit enclosure of E X | claimed, with gaps (H5 and the Schauder wrapper not written in full) | rdeA/EXISTENCE.md | not yet |
| N11 | rdeA measurements (float64 only, nothing interval-verified): contraction constant 0.571 in a weighted sup norm on quantile functions (spectral value 0.351), defect 1.9e-4, radius 6.8e-5, conditional enclosure c in [2.2465, 2.2526]; H3 has 6-14% room | evidence | rdeA/measure_variants.log | none |
| E6 | rdeB Theorem M: on increasing-failure-rate laws, T is monotone in the dispersive order (which ignores shifts); every output of T is IFR | verified | rdeB/ROUTES.md | v9: holds, no fix; needs only IFR of one variable (Prekopa twice); about 230,000 level comparisons, no grid-stable violation; a designed non-IFR control does violate, so the test has power |
| E7 | rdeB Corollary A: Q_n = T^n(delta_0) is below every fixed point in dispersive order, hence c* <= c_n = 3(E T(Q_n) + 2 E Q_n) for every n | verified | rdeB/ROUTES.md | v9: holds; the 'hence' is Lemma 3 (centre both laws; dispersive order with equal means gives convex order) |
| E8 | rdeB Corollary B: existence of a fixed point (monotone limit of Q_n under an explicit super-solution from the scaling law T(aQ)(theta) = a·T(Q)(theta/a^3)) | verified with cosmetic fixes | rdeB/ROUTES.md | v9 (97% for existence plus cost <= c_n): drop one from-memory citation (not needed); the limit's cost bound needs no continuity |
| E9 | Certified bound: every fixed point has cost 9·E X <= 2.2519 (run B, stage 10, prints 2.251887; run A gives 2.274755 at stage 5), by iterating T from delta_0 with one-sided rounding to piecewise-exponential IFR laws and interval arithmetic | verified | cert/CERTIFICATE.md, cert.py | v10 (97%): rounding safety, bounds on g and g', final bound and arithmetic all hold; own code recomputes g and g' at all 4,000 points of all 7 transitions of run A (theta below g, s above g' everywhere) and gives c = 2.25188597 for the last law of run B; both runs reproduce; run B's laws regenerated in v10/. Codex source audit: scheme sound; quote 2.2519 because the printed digit is round-to-nearest |
| W1 | Theorem 1 written out in full (closes the two Codex gaps; new Lemma G for the 1/p cut-off; explicit union bounds; Poincare inequality proved in the text; Lemma V on uniqueness of pair loads) | verified with one fix applied | final/THEOREM_1.md | v11 (independent second read, 93%): all points hold; one arithmetic slip in Theorem A Step 6 (a dropped factor p; the constant kappa too large); the orchestrator corrected the display and set kappa = 0.3^{r_A}/(8·4^{r_A+1}), a weaker constant than the referee's suggestion; nothing downstream changes |
| W2 | Theorem I certificate hardened: verified exp, log and Euler constant; outward-rounded decimals; exported constant asserted against an exact rational. Psi_0 <= -3.84751328322; -Psi_0/theta >= 2.30224586; perturbed exponent <= -0.000592983 at beta = 2.30. All original conclusions stand | written | final/certI/ (NOTES.md, certify_rigorous.py, constants_rigorous.py) | none yet. The library defect was reproduced at 80 bits (iv.exp: 21 of 100,000; iv.log: 1 of 1,000,000); the original log line for Psi_0 was 2e-12 below its own endpoint, unused downstream |
| W3 | Tree process and locality written out: construction of the recursive tree process on the infinite-degree tree (moments from the fixed-point equation, cylinder consistency, root-independence, invariance, perfect matching, cost and normalisation); decomposition theorem with conventions; Lemma D with hypotheses; LP locality step 3 replaced by total-variation convergence; the unused repair corollary corrected | written; second reads launched | final/TREE_AND_LOCALITY.md | v14 (locality), v15 (tree process, existence, certificate) launched |
| R2a | Second independent read of the transfer and the join (I2, I3) | holds (95%) | sc/TRANSFER.md, interface/THEOREM_I_PA.md | v13: every link of the chain re-derived; the planted comparison's factor contains no dependence on the candidate assignment's pair values; Janson bound re-derived and simulated; exponent recomputed in exact rationals (<= -0.0005929831 at beta = 2.30); cosmetic: prefactor inside the supremum, the coefficient in Lemma Q is <= 0 rather than = 0, no table row for beta = 2.30224 (one supplied), and the 'size of q' remark is understated by many orders of magnitude (q above about 10^52,900 at beta = 2.30) |
| R2b | Second independent read of the LP locality theorem (F1, with the total-variation replacement of its step 3) | holds (92%) | lplocal/LP_LOCALITY.md, final/TREE_AND_LOCALITY.md (d), localconv/LOCALCONV.md §1-2 | v14: all seven points hold; the construction needs one ball of radius r+1 with K = r and bounded functionals, so Theorem A suffices; brute force of the patch at q = 4 (40 of 40); building the rule on the auxiliary uniform market is essential and is what the proof does; opened arXiv:1901.07167v1 and confirmed Theorem 2(b) in expectation. Flag: with Theorem 2, the theorem forces the cheap tree matching not to be a function of the marked tree; no file claims that it is (checked by search) |
| R2c | Second independent read of Theorem I and its hardened certificate (I1, W2) | holds (93%) | firstmoment/THEOREM_I.md, final/certI/ | v12: worked backwards from the final inequality; re-derived the conditional first moment (exact on 240 cases; 0 violations in 59,584), the concavity (adversarial search finds no violation; one apparent violation was the referee's own artefact), the limit step; read certify_rigorous.py line by line (verified exp, log, Euler constant, exact rational inputs, export) and re-ran it; recomputed Psi_0 = -3.847516394 by two other methods. Fixes: keep the prefactor inside the maximum; cite certify_rigorous.py in section 5. Note: no slack beyond 2.30224 |
| R2d | Second independent read of the tree process, existence and the cost certificate (W3 (a)-(c), E6-E9) | holds: (A) 93%, (B) 92%, (C) 90% | final/TREE_AND_LOCALITY.md, rdeB/ROUTES.md, cert/ | v15: no error; invariance proved is the team-level mass-transport identity, which is what Theorem 3 needs; cost formula re-derived without the tree; own c_1..c_10 agree with v9 and v10 to all digits; certificate runs reproduce digit for digit; one citation in Lemma 4 to replace by a direct argument; one numerical loose end on non-IFR pairs (off the proof path) |
| R2e | Further independent read of the tree-process construction (W3 (a)-(c)) | holds (92%) | final/TREE_AND_LOCALITY.md | v16: all five claims and parts (b), (c) hold; G finite and the law atomless are derived, not assumed; outward messages are functions of the inward family; simulation on truncated trees (about 148 million agent checks) finds 0 inconsistencies and cost 2.2452 ± 0.0021 against a truncated prediction 2.2482; one clarifying sentence owed in Theorem A6 Step 1. The orchestrator checked that hypothesis (H0) of LP_LOCALITY.md is word for word the identity (MT) proved there |

## Orchestrator notes after three reports (sister, tree, redteam)

- The conjecture is not a statement about the tree alone (T2). The working architecture is
  now: (I) an integer lower bound beta from a first moment conditioned on local stars (T5);
  (II) an LP upper bound alpha from a local fractional rule on the tree, carried to finite
  q; then E D >= (beta - alpha - o(1))·q^(1/3). Candidate numbers: beta = 2.3025,
  alpha = 2.2495 (replica-symmetric matching) or 2.197 (fractional fixed point).
- Finite-q data cannot see the limit: at cost level K the chance that a near-tight team
  shares a pair with another is about K²·q^(-1/3), which is not small until q is far above
  10^4. So the pair-additive numbers (2.53, 2.57) are pre-asymptotic, and the sister model
  without pair-sharing (N1) is the better guide to the limit constants.
- Step (II) needs a patching argument. Aldous's device for the assignment problem (split
  each cost into two independent parts, use one for the main solution and the other for
  patching) is the natural tool; the sprinkle line may already be on it.

Round 2 plan and the proof architecture: see ROUND1_SYNTHESIS.md.

## Reads by a different model family (Codex, via cross-model; 2026-11-11)

Full texts: `codex/CODEX_READ_integer_bound.md`, `codex/CODEX_READ_theorems_1_3_4.md`.

- **Theorem 2 (integer bound; I1, I2, I3):** no substantive mathematical error; confidence
  95% conditional on the interval library. Findings: (gap) the certificate prints interval
  endpoints with a routine that does not guarantee outward rounding; Codex re-ran the
  certificate and confirmed by direct interval comparison that Psi_0 < -3.84751328, so the
  imported decimal is safe; fix the reporting. (gap) the prefactor's dependence on the
  perturbation parameters is stated too narrowly; harmless. (cosmetic) "equivalently" in
  MAIN_THEOREM.md should be "in particular" (fixed). Codex also verified an error budget at
  beta = 2.30224.
- **Theorem 1 (A1-A4):** appears correct after two pieces of bookkeeping are written out:
  the exceptional events (Lemma U failing, |L| < np/2) in the rotation estimate, and the
  passage of the expansion through the 1/p cut-off before U(b) applies. The pair-level
  independence in the rotation argument is confirmed.
- **Theorem 3 (F1, L1-L3):** appears correct; the patching survives. Gaps: the "repair
  corollary" of LOCALCONV.md is not justified for an arbitrary rule (LP_LOCALITY.md does not
  use it in that form); the l1/l2 approximation step needs the per-shape argument (same
  point as v5).
- **Theorem 4 (T1-T2, E1):** the decomposition is correct (add the convention at x_f = 1).
  The construction of the two-directional message process on the infinite-degree tree in
  RDE.md §1.3 is compressed and needs a cylinder-consistency argument; given it, the
  matching and the cost formula c = E S_+^4/8 = 9·E X check out. Codex names this the
  single weakest step of the second packet.
- **Import P3:** Codex opened the primary PDF: Frieze-Pegden-Tkocz Theorem 2(b) is for
  independent uniform pair costs and bounds the EXPECTED cost of a completed greedy
  assignment by c_2·n^(1/3). This is the form the proofs use.
- **Derivation of the main theorem:** correct given the components; apply Chebyshev with
  t·E D for every fixed t to get D/E D -> 1 (fixed in MAIN_THEOREM.md).
- **Existence and the cost certificate (E1, E6-E9); third Codex read, a source audit without
  re-running:** no reversed inequality found. Verdicts: invariant matching and cost
  identities sound (95%); Theorem M, the cost comparison and existence sound after a minor
  moment clarification (94%); the certificate's enclosure scheme appears sound for runs A
  and B (90%). Findings: (gap) the printed "2.251887" comes from a round-to-nearest format,
  so the last digit is unsupported; quote 2.2519 (any internal bound that prints as
  2.251887 is below 2.2518875), or assert against an exact rational; the margin to 2.30 is
  unaffected. (gap) ROUTES.md needs all absolute moments of a fixed point; add
  E X_-^k <= E Y_+^(k+3)/((k+1)(k+2)(k+3)). (gap) cert.py's exact-argument justification
  holds for the advertised lattice parameters only; restrict the accepted parameters.
  (cosmetic) "all exponential moments" should say positive parameters; the super-solution's
  multiplier is existential, not numerical. Full text:
  `codex/CODEX_READ_existence_and_certificate.md`.
- **Fourth Codex read (the final write-outs: final/THEOREM_1.md, final/TREE_AND_LOCALITY.md,
  final/certI/):** every gap from the earlier reads is now closed: the rotation's
  exceptional events, the expansion through 1/p (Lemma G), the corrected Step 6 and its
  constant kappa, the two-directional cylinder consistency of the tree process, the
  decomposition conventions and weak duality, and the LP approximation step (by the
  total-variation argument). Verdicts: Theorem 1 96%; tree process, matching and cost 94%;
  parts (b)-(e) 94%; certificate logic 93% (not re-run by Codex; re-run by v12). New
  findings, all in parts not used on the main line and now corrected in the file: a false
  inequality and a division by zero in the optional §D.3; an insufficient independence
  hypothesis in the side corollary §E. Full text: `codex/CODEX_READ_final_writeouts.md`.

## Final status (12 November 2026)

Stopped under the rule "a proof survives verification". 38 of 40 agent runs used, plus four
Codex reads.

| Part of the proof | Author | Same-family adversarial reads | Codex |
|---|---|---|---|
| Theorem 1 (variance) | r1-pairdeficit; written out by fixA | v3, v11 (one constant fixed) | read 2 (gaps), read 4 (closed, 96%) |
| Theorem I (independent teams) + certificate | r2-firstmoment; hardened by fixC | v6, v12 | read 1 (95%), read 4 (logic 93%) |
| Transfer and join (Theorem 2) | r2-sc, r3-interface | v7, v8, v13 | read 1 |
| Local convergence (Theorem A) | r2-localconv | v4, v14 (its five steps) | read 2 |
| LP locality (Theorem 3) | r2-lplocal; step 3 replaced by fixB | v5, v14 | read 2, read 4 |
| Tree: decomposition, process, matching, cost (Theorem 4) | r1-tree, r2-rde; written out by fixB | v1, v9, v15, v16 | read 2, read 3, read 4 (94%) |
| Existence in the dispersive order | r3-rdeB | v9, v15 | read 3 (94%) |
| Cost certificate (9·E X <= 2.2519) | r4-cert | v10, v15 | read 3 (source audit, 90%) |

No read found an error that affects a theorem. Fixes applied or owed are in
`final/ERRATA.md`. Trusted base: the interval library's four basic operations, IEEE
arithmetic in the BLAS sums of the cost certificate, and Frieze-Pegden-Tkocz Theorem 2(b).
No person has read any part.
