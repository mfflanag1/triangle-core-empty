# Cavity computation of a_IP and a_LP in the local limit (2026-11-11)

Lane agent (Claude Code), 60-minute time box. **Non-rigorous throughout**: replica-symmetric cavity
plus population dynamics, with floating-point numerics. Not reviewed by another model family.
Interpreter `/Users/maxf/miniforge3/bin/python3`, at most 4 single-threaded workers, no API spend.

## Summary

| quantity | value | status |
|---|---|---|
| a_RS: integral RS fixed point, cost per team | **2.2495 +- 0.0005** (E X = 0.2500, so a = 9 E X = 9/4 to 4 digits) | deterministic grid solution, two energy formulas agree |
| a_LP(tree): fractional fixed point (n-fold relaxation, n -> inf) | **2.197 +- 0.002** (n = 2, 4, 8: 2.2106, 2.2003, 2.1975) | population dynamics, drift-corrected |
| a_RS - a_LP(tree) | **0.052 +- 0.003** | |
| integral RS fixed point stable against fractional perturbations? | **No.** Kinks grow by about x2.2 per generation and saturate at the fractional fixed point | control: for bipartite matching (d = 2) kinks decay and the energy returns to pi^2/6 |
| finite-q E(LP* - T)/q^(1/3), fit a q^(1/3) + b + c q^(-1/3) | 0.156 +- 0.026 (2-term fit: 0.097 +- 0.003). With the slope fixed at 0.052 the fit is poor (chi2 19.7 on 4 dof) | data from STRUCTURE.md section 1 |

Reading (heuristic). The LP is convex, so I take its local limit to be the fractional fixed point:
a_LP = 2.197. The integer problem's RS fixed point (2.2495) is unstable, both in relaxation space
(this note) and to bug proliferation (section 4). So the integer problem is probably
replica-symmetry-broken, with a_IP >= a_RS = 2.2495. That gives **a_IP - a_LP >= 0.052**, which
supports the two-constants conjecture. The finite-q gap data prefer a larger constant, about 0.1-0.16,
which fits an RSB correction that raises a_IP above a_RS. I did not compute the 1RSB value of a_IP.

Absolute levels do **not** match the naive finite-q ratios. 3q - T and 3q - LP* divided by q^(1/3)
sit at 2.58 and 2.53, which is 13-15% above the cavity numbers. A fit with a q^(-1/3) correction
moves both to the cavity values (3q - LP*: a = 2.205 +- 0.255; 3q - T: a = 2.36 +- 0.26). But
those fits are ill-conditioned, so the scaling validation is weak (section 5).

## 1. Local limit

**Scale.** Pair deficits d = 1 - (pair value) are iid U[0,1], so the team deficit
D = d_ab + d_bc + d_ac has P(D < t) = t^3/6 exactly for t <= 1. An agent is in q^2 teams, so the
expected number of its teams with D < t is q^2 t^3 / 6. This is O(1) exactly when t ~ q^(-2/3).
Rescale xi = q^(2/3) D. An agent's incident team costs then form, in the limit, a Poisson process
of intensity xi^2/2 d xi on [0, inf). Since 3q - T = sum over M's teams of D = q^(1/3) x
(mean xi of a matched team), a_IP is the expected rescaled cost of the root agent's team, and
likewise for a_LP. (Check: STRUCTURE.md's per-pair deficit x q^(2/3) of 0.86 gives a team cost of
2.58, which is the same number as 3q - T over q^(1/3).)

**Reconciling with STRUCTURE.md's "4.5 q eps^(3/2)".** That count is in *reduced cost* relative
to the LP dual p_lp, not in raw deficit. Write p(f) = 3 - (u_a + u_b + u_c). Near-tight then means
D_f < u_a + u_b + u_c + eps, with u's of order q^(-2/3) and a sum of about 2.5 in rescaled units.
The count is E[(U + s)^3]/6 with s = eps q^(2/3) and U the rescaled price sum. Over s in
[0.5, 3] this is locally well fitted by s^(3/2): (2.55 + 1)^3/6 = 7.4 against 6.6 observed at
s = 1, and (2.55 + 2.64)^3/6 = 23 against 22.3 observed at eps = 0.2, q = 48. So "q eps^(3/2)" is
(eps q^(2/3))^(3/2), a function of the same scaling variable. The q^(2/3) scale and the s^3/6
raw-cost intensity are consistent.

**Pair sharing does not survive the limit, but it dies slowly.** Fix a near-tight team (a,b,c), so
each of its pair deficits is O(q^(-2/3)). Another team sharing pair (b,c), say (a',b,c), is near-tight
only if d_a'b + d_a'c = O(q^(-2/3)). That has probability O(q^(-4/3)), and there are q choices of a'.
So each near-tight team shares a pair with O(q^(-1/3)) other near-tight teams. Within one agent's
star, the indicators for (b,c) and (b,c') are dependent through d_ab. The extra second moment is
q^3 t^5 = O(q^(-1/3)), so the star process is asymptotically Poisson. The limit object is therefore
a 3-uniform Poisson-weighted hypertree with independent team costs. Every agent has the same
offspring law, and the three sides are exchangeable. **The leading finite-size correction is
relative order q^(-1/3)** (0.28 at q = 48). That is large, and it is why absolute levels are hard to
check against q <= 48.

## 2. Recursions

**Integer (n = 1).** For an agent j at the top of its subtree, let X_j = C(j must be matched inside)
- C(j is taken by the parent team). Root i has child teams k with cost xi_k and other agents
j_k, l_k:

    X_i = min_k ( xi_k - X_{j_k} - X_{l_k} ),    so    P(X > x) = exp( -E[(x + X' + X'')_+^3] / 6 ).

A team with cost xi and incoming messages X1, X2, X3 is matched iff xi < S = X1 + X2 + X3. By
Mecke/Campbell this gives two energy formulas, a = E[S_+^4]/8 and the Bethe form
a = 3 E X + E[S_+^4]/12. The normalisation is E[S_+^3]/6 = 1, since every agent is in exactly one
team. Together these force a = 9 E X at the fixed point.

**Fractional (LP).** The message becomes a function. h_j(x) is the cost of j's subtree when the parent
team takes weight x in [0,1] of j, minus the cost at x = 0. It is convex with h_j(0) = 0. In the
integral case h_j(x) = -X_j x. The recursion is

    G_k(x) = xi_k x + h_{j_k}(x) + h_{l_k}(x),   Phi = inf-convolution of the G_k,   h_i(t) = Phi(1 - t) - Phi(1).

The edge marginal is x*_f = argmin_x [xi x + h_a + h_b + h_c], and the energy is
a = int (xi^2/2) xi E x*(xi) d xi. Numerically, weights are restricted to {0, 1/n, ..., 1}: this is
the n-fold relaxation, with each h stored as n nondecreasing cell slopes. The inf-convolution is then
"take the n smallest cell slopes over all children, sorted", which is exact for the n-fold problem.
Letting n -> inf approaches the LP. n = 1 is the integer RS recursion.

On a *finite* hypertree the LP is integral (balanced matrix), and linear h's map to linear h's, so
the integral solution is always a fixed point. The question is whether it is the *attracting* one
on the infinite tree. The same mechanism makes the independent-set LP on random regular graphs
(all-1/2) differ from the tree's integral recursion (Takabe-Hukushima-type analyses of LP
relaxations).

## 3. Results (`out/`)

**RS integral** (`rs_fixed_point.py`, grid h -> 0, Richardson): E X = 0.24995, sd X = 0.8632,
normalisation 1.000000, a_direct = 2.24954, a_bethe = 2.24954. Population dynamics at N = 200k
(`precise_n1_linear_N200k_partial.log`, drift-corrected by `regress.py`) gives 2.2489 +- 0.0010. This
validates the population-dynamics estimator.

**Fractional fixed point** (`precise.py`, N = 80k, 45 sweeps, 15 burn-in; `regress.py` regresses the
per-sweep energy on (norm - 1) and reports the intercept, with SE inflated for autocorrelation):

| n | a(norm = 1) | mean kink s_max - s_min | teams with x* > 0 per agent | teams with x* = 1 per agent |
|---|---|---|---|---|
| 1 | 2.2489 +- 0.0010 | 0 | 1.00 | 1.00 |
| 2 | 2.2106 +- 0.0009 | 0.52 | 1.69 | 0.31 |
| 4 | 2.2003 +- 0.0007 | 0.78 | 2.22 | 0.15 |
| 8 | 2.1975 +- 0.0007 | 0.90 | 2.53 | 0.09 |
| 16 | about 2.206 +- 0.010 (N = 20k, Bethe estimator via `estimate.py`, not drift-corrected) | 0.97 | | |

The successive differences 0.040, 0.010, 0.003 shrink roughly like n^(-2). Extrapolating gives
**a_LP(tree) = 2.197 +- 0.002** (statistical). Finite-N systematics are not separately quantified;
the n = 1 check suggests they are at or below 0.001.

**Instability** (`out/run_n4_tiny.log`): starting from the integral fixed point plus kinks of size
1e-4 (n = 4), the mean kink grows 4.7e-4, 1.1e-3, 2.5e-3, 5.7e-3, 1.3e-2, ... (about x2.2 per sweep,
with one sweep about one generation). It saturates at 0.77, the same fixed point reached from
O(1) random kinks. So the fractional fixed point looks unique and attracting, and the integral one
is linearly unstable.

**Control, d = 2** (`control_d2.py`): bipartite matching on the PWIT, where the LP is integral and the
answer is pi^2/6 = 1.6449. From O(1) kinks the mean kink decays (0.34 -> 0.14 over 25 sweeps,
marginal and slow) and the energy is 1.6403 +- 0.0029. The method therefore does not manufacture
fractional fixed points where none should exist.

**Comparison with the finite-q LP structure.** LP-support teams per agent (counting the M-team):
the data give about 1 + 1.73 = 2.7 at q = 48, still rising. The cavity gives 2.53 at n = 8, also
still rising with n. The data's non-half-integral vertices with large denominators at large q fit
the cavity's increasing kink and support with n, and the n = 2 (half-integral) value is clearly not
the limit.

## 4. Why the integer RS is unstable, and the RSB caveat

At zero temperature X_i depends only on the argmin team's two messages, with derivatives -1 and -1.
A perturbation at depth g therefore reaches the root along 2^g argmin paths, each with derivative
+-1, and the spin-glass susceptibility grows like 2^g. For d = 2 there is one path (1^g, marginal),
which is the known bipartite case where RS is exact. For d >= 3 this is the bug-proliferation
instability. It is consistent with Martin-Mezard-Rivoire's frozen glassy (RSB) phase for
multi-index matching with iid d-tuple costs (cond-mat/0407623, 0507180; that is a different cost
density, r = 0, so not the same constants). The relaxation-space growth of about x2.2 per
generation above is the same instability seen from the LP side. Consequences:
- a_IP is **probably not** 2.2495. With RSB the ground-state energy is expected to be at or above the
  RS value (as for SK), so I use a_IP >= 2.2495 as a heuristic bound, not a theorem. A 1RSB
  (energetic, population-of-populations) computation is the next step. My first attempt had a
  non-stationary energy estimator, with noise of about +-1.5, and was discarded.
- The LP is convex. The working assumption, standard but unproved here, is that the attracting
  fractional fixed point gives its exact local limit. Strictly, the tree computation only shows that
  some fractional object on the limiting tree costs 2.197 per agent. Whether finite-q LP optima
  converge to it is the usual sofic/objective-method gap (Aldous and Wastlund close it for bipartite
  matching; nothing closes it here).

## 5. Finite-q comparison (`fit_finite_q.py`; data from STRUCTURE.md section 1, q = 12..48)

| series | a q^(1/3) + b | a q^(1/3) + b + c q^(-1/3) | a fixed at the cavity value |
|---|---|---|---|
| 3q - T | a = 2.63 +- 0.03 (chi2 7.6/4) | a = 2.36 +- 0.26 | a = 2.2495: b = 2.09, c = -3.18, chi2 6.7/4 (acceptable) |
| 3q - LP* | a = 2.53 +- 0.03 (chi2 8.8/4) | a = 2.21 +- 0.26 | a = 2.197: b = 2.01, c = -2.77, chi2 7.2/4 (acceptable) |
| LP* - T | 0.097 +- 0.003 (chi2 8.8/4) | 0.156 +- 0.026 (chi2 3.5/3) | 0.052: chi2 19.7/4 (poor) |

- The level data are consistent with the cavity levels only if the O(1) correction is large
  (b of about +2, c of about -3), so that it cancels against q^(-1/3) and makes the ratio look flat
  at about 2.55 over q in [12, 48]. That is possible (the relative correction is q^(-1/3)), but it
  is not a confirmation. **The check of a_IP against 2.55 does not validate the scaling; it is
  merely not inconsistent with it.**
- The gap is measured much more precisely, because T and LP* are strongly correlated. It points to a
  limit constant above 0.052, in the range 0.1-0.16. That is what one expects if a_LP = 2.197 and
  a_IP = a_RS + (RSB correction of about 0.05-0.1).

## 6. Answer to the lane question

- a_LP (tree, fractional fixed point) = 2.197 +- 0.002. Solid as a number. Its identification with
  lim (3q - LP*)/q^(1/3) is a standard but unproved convexity/RS assumption.
- a_IP: the RS value 2.2495 is exact for the RS ansatz but probably not the true value, since RS is
  unstable. Best heuristic: a_IP >= 2.2495, and the finite-q gap fit suggests about 2.30-2.35.
- **They differ:** a_IP - a_LP >= 0.052 +- 0.003 under the RS-lower-bound assumption, and about
  0.1-0.16 from the finite-q gap fit. Both are positive. This supports Lemma G(ii) and the
  two-constants conjecture.
- How solid: the fractional-vs-integral split on the tree is robust. It is the same from tiny and
  large initial kinks, monotone in n, and absent in the d = 2 control. The absolute connection to
  q <= 48 is weak because of q^(-1/3) corrections. a_IP itself needs a 1RSB computation.

## Next steps

1. Energetic 1RSB for the integer recursion, with a Bethe (stationary) estimator. The zero-complexity
   y gives a_IP(1RSB). Prediction from the gap data: 2.30-2.35.
2. Data at larger q (64, 96) for the gap only, which is cheap with the reduced-cost MILP. If the
   cavity picture is right, E(LP* - T)/q^(1/3) keeps rising toward 0.1-0.16.
3. Proof direction suggested here: show that the integral fixed point is not the minimiser among
   unimodular fractional perfect matchings of the Poisson hypertree. An explicit improving local
   fractional perturbation of positive density gives a_LP < a_RS. Separately, a lower bound
   a_IP >= a_RS (objective method / Wastlund-type) would finish Lemma G(ii) heuristically. The first
   half is a statement about a single tree object and may be provable.

## Files

- `rs_fixed_point.py`: deterministic RS fixed point and energies. `out_rs_density.npy` holds the X density.
- `lp_popdyn.py`: n-fold relaxation population dynamics (`update`, `energy`), used for the instability run.
- `precise.py`: long runs with per-sweep energies. Outputs `out/precise_n{2,4,8}.log` and `.json`, plus final populations.
- `regress.py`: drift-corrected energy (intercept at norm = 1).
- `estimate.py`: Bethe/stationary estimator and fractional diagnostics from a saved population.
- `control_d2.py`: bipartite PWIT control. Outputs `out/control_d2_n8.log` and `.json`.
- `fit_finite_q.py`: finite-q fits.
