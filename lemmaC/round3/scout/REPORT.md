# Scout report (round 3): steering numerics

Status: final for this time box (11 Nov 2026). Everything here is floating point and is evidence only, except the two
remarks marked "argument", which are short derivations I would defend but which no one else has checked.

Scripts (all in this folder): `scout.py` (items 1, 3, 5), `sister.py` and `extra.py --u` (item 2), `core.py` (item 4),
`extra.py` (deficit tails, triangles), `analyze.py` and `fits.py` (tables). Records are in `out/` (3.4 MB).
Seeds: `cg.draw(q, idx, 770031)` for the main runs; tags 770032 to 770034 for sister, core and U-model runs.
LP* and T come from `lemmaC/scaling/cg.py` (column generation, reduced-cost-fixing MILP). The q = 48 run was stopped at
the first 211 indices in order, so there is no selection on run time.

## Readings

1. **Used deficits grow slowly, like a power of log q; the squared-deficit sum is 3.9·q^(−1/3).**
   Largest pair deficit used by M, times q^(2/3): mean 3.13, 3.37, 3.52, 3.76 at q = 16, 24, 32, 48. Largest pair deficit
   with positive LP load: 3.30, 3.67, 3.92, 4.20. Fits: 1.53 + 0.58·ln q for M and 0.96 + 0.85·ln q for the LP; a power
   law would need exponent 0.17 (M) or 0.23 (LP). The fits alone do not separate log from a small power on this range.
   What argues for log-type growth is that the per-pair law collapses: the share of M's 3q pairs with rescaled deficit
   above t = 1, 2, 3, 4, 5 is 0.34, 0.085, 0.016, 0.002, 0.0002 at q = 16, 24 and 32 alike (tail close to
   exp(−0.8·t^1.5)), so the maximum is the extreme of 3q draws from one law. Not bounded, not a visible power.
   Sum of squared pair deficits of M, times q^(1/3): 3.82, 3.87, 3.92, 3.89 (flat). For the LP: sum of load·d² is 1.49,
   1.31, 1.20, 1.02 (about the same as M), and sum of load²·d² times q^(1/3) is 2.83, 2.48, 2.25, 2.06.
   The exact squared gradient of D, the sum of (load − 1_M)²·d², is 0.53, 0.61, 0.64, 0.61, against the 8·E(3q − T) =
   52 to 75 used in fact 2.

2. **The sister model loses integrality much faster.** Sister = each team value an independent sum of three U[0,1]
   (model IH; same marginal as the brief's model, so P(deficit < t) = t³/6 for t ≤ 1). P(LP optimum integral), n = 400:
   0.57, 0.44, 0.38, 0.30, 0.25, 0.125, 0.085, 0.030, 0.045, 0.013, 0.005 at q = 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24.
   Pair-additive model on the same code and fresh seeds: 0.978, 0.927, 0.877, 0.860, 0.777, 0.677, 0.552, 0.477, 0.367,
   0.203, 0.125 (agrees with the brief at 16 and 24). Independent U[0,1] team values (full LP, n = 200): 0.63, 0.345,
   0.215, 0.135, 0.045, 0.015 at q = 4, 6, 8, 10, 12, 16. In all three, −ln P(integral) is roughly linear in q for
   q ≥ 10, with slope 0.12 (pair), 0.22 (IH), 0.35 (U) per unit q. E D in IH is 0.068 at q = 4 and 0.176 at q = 12,
   against 0.0006 and 0.026 in the pair model.

3. **The shortest odd cycle has length 3 and does not grow.** Among instances with D > 0, the shortest odd Berge cycle
   in the fractional part of the LP support has length 3 in 75%, 89%, 97%, 99% of cases at q = 16, 24, 32, 48; the
   remainder are 5 or 7. For cycles in which no team contains a third cycle agent (the ones that witness
   non-balancedness) the shares are 68%, 78%, 89%, 96%. The number of Berge triangles per fractional support is 2.4,
   5.4, 8.8 at q = 16, 24, 32, so roughly linear in q. Only about 10% of triangles are centred on a team of M.
   The fractional part covers 71%, 77%, 83%, 87% of agents; all-weights-1/2 solutions fall from 35% to 1%.
   The LP keeps about 57% of its weight on M's own teams, about 20% on teams keeping one pair of M and about 23% on
   teams sharing no pair with M.

4. **When the core is nonempty it is full-dimensional, with about 6.6q facets and exponentially many vertices.**
   In all 645 instances with D = 0 (q = 4 to 16) the centre point has no tight off-M team and positive margin, so the
   core has dimension 2q (2q − 2 after removing the two sector shifts). Facets (off-M teams tight somewhere on the
   core): mean 19, 27, 33, 48, 78, 96 at q = 4, 5, 6, 8, 12, 16, that is 3.1 to 3.5 per dimension. Teams keeping one
   pair of M are 77% of the facets at q = 4 and 46% at q = 16. Vertices of the core modulo shifts: median 160, 924,
   4076 at q = 4, 5, 6 (mean 202, 1418, 7853; largest 66,674); not attempted beyond q = 6. The largest uniform margin
   eps*, times q^(2/3), has mean 0.161, 0.143, 0.142, 0.102, 0.065, 0.045: the core thins faster than q^(−2/3).

5. **T is sharply concentrated: sd T ≈ 0.85·q^(−1/6).** Y = (3q − T)/q^(1/3) has mean 2.572, 2.586, 2.598, 2.579 and
   sd 0.221, 0.171, 0.155, 0.118 at q = 16, 24, 32, 48 (n = 1000, 1000, 400, 211). Quantiles 90% / 95% / 99%:
   2.85 / 2.94 / 3.12, 2.81 / 2.89 / 3.00, 2.81 / 2.86 / 2.96, 2.75 / 2.80 / 2.83. Skewness 0.14 to 0.32, so close to
   normal with a slightly heavier upper side. sd T times q^(1/6) is 0.88, 0.84, 0.88, 0.82.

## Surprises and points against the brief

- **The brief's "sd of 3q − T is about 0.45" hides a trend.** It falls from 0.556 (q = 16) to 0.428 (q = 48), in line
  with q^(−1/6). This is what the Poincaré inequality predicts once the squared-deficit sum is used instead of 3q − T:
  the bound on sd T is 2·sqrt(3.9)·q^(−1/6).
- **Facts 1 and 2 are loose by a factor of order q^(2/3), and the target can be lowered (argument).** The gradient of D
  in the exponential coordinates is (load_i − 1_M(i))·d_i, so Var D ≤ 4·E Σ(load_i − 1_M(i))²·d_i² ≤ 8·E[d_max·(3q − T)],
  where d_max is the largest pair deficit used by M or by the LP optimum. If a lemma gives d_max ≤ q^(−2/3+o(1)) outside
  an event of probability o(q^(−4/3)), then Var D ≤ q^(−1/3+o(1)) and P(core nonempty) → 0 follows from
  E D ≥ q^(−1/6+ε); a constant lower bound on E D would be enough. **Gap:** the lemma on d_max. A two-team swap only
  gives conditions on off-M pairs, which fact 8 says are biased the wrong way.
  At present sizes the sharper Chebyshev bound is still vacuous: it needs E D above about 1.55, against 0.16 at q = 48.
- **P(core nonempty) falls exponentially in q in the data** (−ln P grows by about 0.12 per unit q; the brief's rows at
  32 and 40 fit the same slope). Any Chebyshev route proves only a polynomial rate.
- **The sister model is not a faithful stand-in.** Pair-additivity protects integrality: a result for independent team
  values will not transfer without using the pair structure, although the direction (the pair model has the smaller
  gap) is the same at every size measured.
- **Uniqueness (argument).** LP* depends on the team weights only through the pair loads. Two different fractional
  vertices can have the same pair loads: on a 2×2×2 block, weight 1/2 on the four teams of one parity class and weight
  1/2 on the four of the other. So when D > 0 the optimal team weights need not be unique, although the pair loads
  are almost surely. When D = 0 the optimum is unique almost surely, because an assignment is determined by its pair
  loads; strict complementarity then gives item 4's first two findings (no off-M team tight in the relative interior,
  dimension 2q) as a fact rather than an observation.
- Nothing else contradicts the brief. E D = 0.041, 0.072, 0.104, 0.161 with fit 0.109·q^(1/3) − 0.238; P(D > 0) =
  0.647, 0.866, 0.963, 1.000; sd D = 0.048 to 0.055; 3q − LP* ≈ 2.54 to 2.57 times q^(1/3).

## Limits

- Odd-cycle search is exact up to length 15; no instance hit the cap or the node budget.
- The LP optimum reported by the solver is one vertex; see the uniqueness remark.
- Facet counts use an LP per candidate team with tolerance 1e-8; vertex counts use cdd in floating point.
- Sister-model D was computed for q ≤ 12 only (the MILP is slow there).

## Tables

### 1. Used deficits (pair deficit x q^(2/3); S2 = sum of squared pair deficits)
| q | n | max d, M: mean / med / 95% / max | max d, LP: mean / med / 95% / max | max team def M (mean) | S2(M) | S2(M) q^(1/3) | sum load d^2 (LP) | sum load^2 d^2 (LP) | sum (load-1_M)^2 d^2 | 3q-T |
|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 1000 | 3.13 / 3.04 / 4.30 / 5.86 | 3.30 / 3.23 / 4.49 / 5.95 | 4.94 | 1.514 | 3.82 | 1.486 | 1.121 | 0.529 | 6.48 |
| 24 | 1000 | 3.37 / 3.28 / 4.51 / 6.12 | 3.67 / 3.61 / 4.74 / 5.65 | 5.20 | 1.342 | 3.87 | 1.310 | 0.858 | 0.606 | 7.46 |
| 32 | 400 | 3.52 / 3.43 / 4.68 / 6.22 | 3.92 / 3.85 / 4.98 / 6.09 | 5.37 | 1.236 | 3.92 | 1.196 | 0.709 | 0.637 | 8.25 |
| 48 | 211 | 3.76 / 3.75 / 4.80 / 7.05 | 4.20 / 4.14 / 5.17 / 5.99 | 5.59 | 1.069 | 3.89 | 1.023 | 0.566 | 0.611 | 9.37 |

Variance check (Poincare, constant 4): bound on sd = 2 sqrt(E S2).
| q | sd T observed | 2 sqrt(E S2(M)) | 2 sqrt(E(3q-T)) (fact 1) | sd D observed | 2 sqrt(E sum (load-1_M)^2 d^2) | sqrt(8 E(3q-T)) (fact 2) | E D | Chebyshev bound on P(D=0) with grad D | with fact 2 |
|---|---|---|---|---|---|---|---|---|---|
| 16 | 0.556 | 2.46 | 5.09 | 0.048 | 1.45 | 7.20 | 0.041 | 1280.6 | 31396 |
| 24 | 0.495 | 2.32 | 5.46 | 0.055 | 1.56 | 7.73 | 0.072 | 468.4 | 11531 |
| 32 | 0.492 | 2.22 | 5.74 | 0.054 | 1.60 | 8.12 | 0.104 | 234.4 | 6074 |
| 48 | 0.428 | 2.07 | 6.12 | 0.051 | 1.56 | 8.66 | 0.161 | 93.8 | 2880 |

### 3. Odd cycles in the fractional part of the LP support (instances with D > 0)
| q | n (D>0) | weak: share 3 / 5 / 7 / >=9 | weak mean | strong: share 3 / 5 / 7 / >=9 | strong mean | frac teams (mean) | frac agents / 3q | all weights 1/2 | min weight (median) |
|---|---|---|---|---|---|---|---|---|---|
| 16 | 647 | 0.75 / 0.22 / 0.03 / 0.00 | 3.56 | 0.68 / 0.29 / 0.03 / 0.00 | 3.71 | 26.1 | 0.71 | 0.35 | 0.3333 |
| 24 | 866 | 0.89 / 0.09 / 0.01 / 0.00 | 3.23 | 0.78 / 0.21 / 0.01 / 0.00 | 3.47 | 46.1 | 0.77 | 0.13 | 0.1171 |
| 32 | 385 | 0.97 / 0.03 / 0.01 / 0.00 | 3.08 | 0.89 / 0.10 / 0.01 / 0.00 | 3.23 | 71.3 | 0.83 | 0.04 | 0.0345 |
| 48 | 211 | 0.99 / 0.01 / 0.00 / 0.00 | 3.02 | 0.96 / 0.04 / 0.00 / 0.00 | 3.08 | 117.3 | 0.87 | 0.01 | 0.0060 |

### 5. Upper tail of Y = (3q - T)/q^(1/3)
| q | n | mean | sd | sd of T | skew | 50% | 90% | 95% | 99% | max | min | E D | sd D | P(D>0) | (3q-LP*)/q^(1/3) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 1000 | 2.572 | 0.221 | 0.556 | 0.14 | 2.56 | 2.85 | 2.94 | 3.12 | 3.27 | 1.90 | 0.041 | 0.048 | 0.647 | 2.556 |
| 24 | 1000 | 2.586 | 0.171 | 0.495 | 0.21 | 2.58 | 2.81 | 2.89 | 3.00 | 3.14 | 2.11 | 0.072 | 0.055 | 0.866 | 2.561 |
| 32 | 400 | 2.598 | 0.155 | 0.492 | 0.21 | 2.59 | 2.81 | 2.86 | 2.96 | 3.06 | 2.23 | 0.104 | 0.054 | 0.963 | 2.565 |
| 48 | 211 | 2.579 | 0.118 | 0.428 | 0.32 | 2.57 | 2.75 | 2.80 | 2.83 | 2.93 | 2.29 | 0.161 | 0.051 | 1.000 | 2.535 |

### 2. Sister model: P(LP optimum integral)
| q | 4 | 5 | 6 | 7 | 8 | 10 | 12 | 14 | 16 | 20 | 24 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| P: P(integral) | 0.978 | 0.927 | 0.877 | 0.860 | 0.777 | 0.677 | 0.552 | 0.477 | 0.367 | 0.203 | 0.125 |
| P: n | 400 | 400 | 400 | 400 | 400 | 400 | 400 | 400 | 400 | 400 | 400 |
| P: E D | 0.0006 | 0.0024 | 0.0068 | 0.0060 | 0.0096 | 0.0154 | 0.0261 |  |  |  |  |
| IH: P(integral) | 0.570 | 0.438 | 0.380 | 0.300 | 0.247 | 0.125 | 0.085 | 0.030 | 0.045 | 0.013 | 0.005 |
| IH: n | 400 | 400 | 400 | 400 | 400 | 400 | 400 | 400 | 400 | 400 | 400 |
| IH: E D | 0.0684 | 0.0828 | 0.0913 | 0.1124 | 0.1239 | 0.1556 | 0.1755 |  |  |  |  |
| U: P(integral) | 0.630 |  | 0.345 |  | 0.215 | 0.135 | 0.045 |  | 0.015 |  |  |
| U: n | 200 |  | 200 |  | 200 | 200 | 200 |  | 200 |  |  |
| U: E D |  |  |  |  |  |  |  |  |  |  |  |

### 4. Core geometry when D = 0 (essential core Q, dimension 2q-2; brief's core has dimension 2q)
| q | n | core nonempty | off-M teams tight at centre (max) | eps* q^(2/3): mean / med | eps*>0 share | facets: mean / med / max | facets sharing 1 pair with M / 0 pairs (mean) | facets/(2q-2) | vertices: med / mean / max | mean width x q^(2/3) |
|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 200 | 192 | 0 | 0.161 / 0.120 | 1.00 | 18.7 / 18 / 44 | 14.4 / 4.3 | 3.12 | 160 / 202 / 1047 | 1.56 |
| 5 | 200 | 194 | 0 | 0.143 / 0.110 | 1.00 | 26.8 / 25 / 53 | 18.9 / 8.0 | 3.35 | 924 / 1418 / 6656 | 1.71 |
| 6 | 50 | 47 | 0 | 0.142 / 0.098 | 1.00 | 32.6 / 31 / 73 | 21.6 / 11.0 | 3.26 | 4076 / 7853 / 66674 | 1.68 |
| 8 | 150 | 125 | 0 | 0.102 / 0.072 | 1.00 | 48.1 / 43 / 88 | 28.7 / 19.4 | 3.44 | n/a | 1.92 |
| 12 | 100 | 58 | 0 | 0.065 / 0.043 | 1.00 | 77.5 / 68 / 173 | 39.7 / 37.8 | 3.52 | n/a | 2.01 |
| 16 | 60 | 29 | 0 | 0.045 / 0.036 | 1.00 | 96.2 / 89 / 170 | 44.5 / 51.7 | 3.21 | n/a | 1.71 |

### 1b. Tails of used pair deficits (same instances): share of the 3q used pairs with d q^(2/3) > t
| q | n | M: t=1 | 2 | 3 | 4 | 5 | 6 | LP load: t=1 | 2 | 3 | 4 | 5 | 6 | share of S2(M) from d q^(2/3) > 2 | > 3 | > 4 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 300 | 0.3383 | 0.0856 | 0.0154 | 0.0019 | 0.0001 | 0.0000 | 0.3357 | 0.0826 | 0.0133 | 0.0017 | 0.0001 | 0.0000 | 0.461 | 0.149 | 0.029 |
| 24 | 200 | 0.3357 | 0.0812 | 0.0156 | 0.0023 | 0.0002 | 0.0000 | 0.3327 | 0.0797 | 0.0142 | 0.0017 | 0.0000 | 0.0000 | 0.454 | 0.155 | 0.036 |
| 32 | 80 | 0.3402 | 0.0880 | 0.0173 | 0.0025 | 0.0003 | 0.0000 | 0.3378 | 0.0838 | 0.0144 | 0.0013 | 0.0001 | 0.0000 | 0.478 | 0.163 | 0.037 |

Extreme-value collapse test from the main runs: -ln P(max_M <= t)/(3q) (q-independent if the used deficits behave like 3q iid draws from one rescaled law)
| q | t=3 | 3.5 | 4 | 4.5 | 5 |
|---|---|---|---|---|---|
| 16 | 0.0153 | 0.0056 | 0.0019 | 0.0007 | 0.0002 |
| 24 | 0.0166 | 0.0061 | 0.0021 | 0.0007 | 0.0003 |
| 32 | 0.0169 | 0.0064 | 0.0025 | 0.0008 | 0.0003 |
| 48 | 0.0205 | 0.0064 | 0.0026 | 0.0007 | 0.0003 |

### 3b. LP weight by pairs shared with M, and Berge triangles in the fractional support (D > 0 instances)
| q | n | LP weight/q on teams sharing 0 pairs with M | 1 pair | 3 pairs (M teams) | n (D>0) | triangles per instance | strong triangles | share of strong triangles whose centre (a,b,c) is an M team | centre in LP support |
|---|---|---|---|---|---|---|---|---|---|
| 16 | 300 | 0.187 | 0.232 | 0.581 | 205 | 2.43 | 1.37 | 0.114 | 0.046 |
| 24 | 200 | 0.233 | 0.200 | 0.567 | 177 | 5.37 | 2.03 | 0.097 | 0.075 |
| 32 | 80 | 0.247 | 0.189 | 0.564 | 77 | 8.75 | 3.00 | 0.078 | 0.065 |

### Fits

```
maxdM means [3.128 3.367 3.522 3.763] se [0.019 0.02  0.032 0.042] | log fit 1.53 + 0.58 ln q, resid 0.007 | power fit exponent 0.170, resid 0.013 | (ln 3q)^(2/3) fit -0.34 + 1.41 x, resid 0.008
maxdL means [3.298 3.67  3.922 4.199] se [0.02  0.018 0.029 0.038] | log fit 0.96 + 0.85 ln q, resid 0.046 | power fit exponent 0.227, resid 0.067 | (ln 3q)^(2/3) fit -1.81 + 2.08 x, resid 0.038
q [16 24 32 48]
sd T [0.556 0.495 0.492 0.428] x q^(1/6) [0.882 0.84  0.876 0.816] se of sd [0.012 0.011 0.017 0.021] loglog slope -0.221
E S2(M) q^(1/3) [3.815 3.872 3.923 3.886] loglog slope -0.315
E gradD2 [0.529 0.606 0.637 0.611]
E sum load^2 d^2 x q^(1/3) [2.825 2.475 2.251 2.058]
E D [0.0406 0.0719 0.1042 0.1614] se [0.0015 0.0017 0.0027 0.0035] ED/q^(1/3) [0.0161 0.0249 0.0328 0.0444]
E D fit 0.109 q^(1/3) + -0.238
corr(D, defT) [np.float64(0.12), np.float64(0.16), np.float64(0.17), np.float64(0.11)]
P q [ 4  5  6  7  8 10 12 14 16 20 24] -ln P(int) [0.02 0.08 0.13 0.15 0.25 0.39 0.59 0.74 1.   1.6  2.08] slope per unit q over q>=10: 0.124
IH q [ 4  5  6  7  8 10 12 14 16 20 24] -ln P(int) [0.56 0.83 0.97 1.2  1.4  2.08 2.47 3.51 3.1  4.38 5.3 ] slope per unit q over q>=10: 0.225
U q [ 4  6  8 10 12 16] -ln P(int) [0.46 1.06 1.54 2.   3.1  4.2 ] slope per unit q over q>=10: 0.353
```
