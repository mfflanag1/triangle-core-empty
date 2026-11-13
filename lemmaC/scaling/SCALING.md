# Scaling of the gap D = LP* - T (random triangle market, uniform pair values), 2026-11-11

Lane agent (Claude Code), 60-minute time box, machine shared (4 workers max). Not reviewed. **Floating point only
(scipy HiGHS dual simplex + scipy milp). Evidence, not certificates.**

## Files
- `cg.py`: LP* by column generation (restricted equality LP, exact dense pricing of all q^3 reduced costs,
  stop at max r <= 1e-9), then T by reduced-cost fixing as in `../structure/mustruct.py` (MILP on {r_f >= -delta},
  delta0 = 0.3 for q > 48, doubled until T_inc >= LP* - delta). Seeds for new runs:
  `SeedSequence([20261111, 0, 0, q, idx, 7])` (tag 7, a sixth entry, so no collision with gate-1 5-entry seeds).
  `--gate1` drops the tag and reproduces gate-1 instances for validation.
- `analyse.py` -> `out/analysis.txt`: the table and fits below. q = 16..48 come from `../structure/out/struct_q*.jsonl`
  (gate-1 seeds, full LP); q = 56, 64 from `out/scal_q*.jsonl` (tag 7).
- `out/val_full_q032.jsonl`, `out/val_q040.jsonl`, `out/val_q048.jsonl`: validation records.

## Validation of the speed-up
- Column generation vs the full q^3 LP, q = 32 (4 gate-1 instances): |LP_cg - LP_full| <= 3e-14.
- Against `structure/` records (full LP + MILP), gate-1 idx 0..11 at q = 40 and 48: max |dLP*| = 1.1e-13,
  max |dT| = 9e-11.
- In every run (q = 32..80) the initial column set (top-3 third member for every pair, in all three directions,
  plus the diagonal; about 9q^2 columns) was already LP-optimal: one CG iteration, no priced-in columns.
  The LP is not the bottleneck (<= 4 s at q = 80).
- **The MILP is the bottleneck.** Median wall time per instance (machine loaded): 19 s at q = 56, 35 s at q = 64,
  about 240 s at q = 80 (candidate set ~7,200 teams at delta = 0.3; no delta doubling was needed at q >= 56).
  A q = 96 test instance had not finished after 9 minutes and was stopped. **q = 96 and 128 were not reached.**

## Table

Counts: q = 16..40: 200; q = 48: 120 (gate-1 seeds, from `structure/`); q = 56: 70; q = 64: 23 (idx 0..22, tag 7;
runs stopped at the time box, instances returned in index order so no selection on run time); q = 80: 2 (listed
separately, not in fits).

| q | n | E D | sd D | se E D | P(D>0) | E D/q^(1/3) | E D/q^(1/4) | E D/q^(1/6) | E(3q-T) | /q^(1/3) | E(3q-LP*) | (3q-LP*)/(3q-T) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 200 | 0.0389 | 0.0462 | 0.0033 | 0.645 | 0.0154 | 0.0195 | 0.0245 | 6.438 | 2.555 | 6.400 | 0.9940 |
| 24 | 200 | 0.0694 | 0.0534 | 0.0038 | 0.875 | 0.0241 | 0.0314 | 0.0409 | 7.477 | 2.592 | 7.408 | 0.9907 |
| 32 | 200 | 0.1082 | 0.0544 | 0.0038 | 0.965 | 0.0341 | 0.0455 | 0.0607 | 8.207 | 2.585 | 8.099 | 0.9868 |
| 40 | 200 | 0.1273 | 0.0528 | 0.0037 | 0.985 | 0.0372 | 0.0506 | 0.0689 | 8.897 | 2.601 | 8.769 | 0.9857 |
| 48 | 120 | 0.1564 | 0.0537 | 0.0049 | 1.000 | 0.0430 | 0.0594 | 0.0820 | 9.346 | 2.572 | 9.190 | 0.9833 |
| 56 | 70 | 0.1879 | 0.0428 | 0.0051 | 1.000 | 0.0491 | 0.0687 | 0.0961 | 9.860 | 2.577 | 9.672 | 0.9809 |
| 64 | 23 | 0.2117 | 0.0350 | 0.0073 | 1.000 | 0.0529 | 0.0749 | 0.1059 | 10.290 | 2.572 | 10.078 | 0.9794 |

q = 80 (n = 2): D = 0.249, 0.201; 3q - T = 11.34, 11.07; 3q - LP* = 11.09, 10.87. Consistent with the trend
(mean 0.225, ratio about 0.980), but two points carry no weight.

## Fits (bootstrap over instances within each q, 4000 draws)

Power law E D ~ q^b (OLS on log means):

| q range | b | 90% | 95% |
|---|---|---|---|
| 16-64 | 1.21 | [1.13, 1.30] | [1.11, 1.32] |
| 24-64 | 1.11 | [1.02, 1.20] | [1.01, 1.22] |
| 32-64 | 1.00 | [0.90, 1.10] | [0.88, 1.12] |
| 40-64 | 1.10 | [0.95, 1.24] | [0.93, 1.27] |

Local slopes: 1.43, 1.54, 0.73, 1.13, 1.19, 0.89 (16->24 ... 56->64). Weighted pure power over 16-64: b = 1.13.

- **D/q^(1/3) is still rising**: 0.015, 0.024, 0.034, 0.037, 0.043, 0.049, 0.053. No plateau in range.
- Offset fits (weighted by se): E D = 0.112 q^(1/3) - 0.248 (chi2 13.6 / 5 dof); 0.199 q^(1/4) - 0.364
  (chi2 17.2); 0.116 log q - 0.290 (chi2 31.6, poor). Pure power fits best (chi2 10.4 / 5 dof).
- 3q - T ~ q^0.336 [0.327, 0.345]; 3q - LP* ~ q^0.325 [0.316, 0.334]. 3q - T = 2.57 q^(1/3) flat (2.56-2.60).
- **Ratio (3q - LP*)/(3q - T)**: 0.9940 -> 0.9794, i.e. 1 - ratio = E D / E(3q - T) = 0.0060, 0.0093, 0.0132,
  0.0143, 0.0167, 0.0191, 0.0206: still falling (ratio) / rising (1 - ratio), roughly like q^0.9 over this range,
  so the exponent of 3q - LP* sits ~0.01 below that of 3q - T. No sign of convergence of 1 - a_LP/a_IP yet; the
  q^(1/3)-plus-offset reading gives a_IP - a_LP of about 0.11 (about 4% of a_IP = 2.57) as an asymptote.
- sd D: 0.046-0.054 for q = 16-48, then 0.043 (n = 70) and 0.035 (n = 23). Flat or slightly falling, never rising.
  E D / sd D = 4.4 at q = 56 and 6.0 at q = 64.

## Reading

1. **Over q = 16-64, E D grows like q^1.0-1.2, and every bootstrap interval excludes b <= 1/4 and b <= 1/6 by a
   wide margin** (lowest 95% endpoint 0.88). E D/q^(1/4) rises 0.020 -> 0.075 monotonically.
2. This exponent **cannot be asymptotic**: D <= S - T = O(q^(1/3)) a priori, and 3q - T itself is only 2.57 q^(1/3).
   So the data are pre-asymptotic in the sense that the gap is still filling in toward its eventual scale; the
   observed b ~ 1.1 says nothing directly about the limiting exponent except that the approach has not stalled.
3. What the data do support: no saturation or turnover up to q = 64 (the last local slope, 0.89, is within noise of
   the others; se at q = 64 is 0.007). The offset model E D ~ 0.11 q^(1/3) - 0.25 fits acceptably and is the reading
   consistent with the two-constants conjecture (a_LP < a_IP by about 4%). The range q <= 64 cannot rule out a later
   slowdown to an exponent below 1/4 (e.g. a gap growing like log q eventually); log q fits clearly worse here.

## Implication for Proposition R

Proposition R needs E D / q^(1/4) -> infinity (or / q^(1/6) with the q^(1/3) bound on E(3q - T)). Numerically:
- E D / q^(1/4) and E D / q^(1/6) rise steadily through q = 64 with no flattening; the trend supports b > 1/4 and
  b > 1/6 over the accessible range.
- But the measured b (~1.1) is a transient on the way to an exponent that is at most 1/3, so the data cannot certify
  the asymptotic hypothesis of Proposition R. The best-fitting asymptotic model (c q^(1/3) - d, c ~ 0.11) is exactly
  the two-constants conjecture, which would give the q^(1/3) rate, comfortably above q^(1/4).
- The rigorous variance budget is the weak point, not the mean: empirical sd D is about 0.04-0.05 (flat), versus the
  O(q^(1/6)) Chebyshev scale used in Prop. R. At q = 64 the actual ratio E D / q^(1/4) is only 0.075, so
  Proposition R as stated becomes numerically effective only at astronomically large q; the proof route must rest on
  the a_LP < a_IP limit (cavity / objective method), not on extrapolating these finite-q exponents.

## Not done
- q = 96, 128: MILP too slow under the time box (> 9 min per instance at q = 96 on the loaded machine). LP* alone is
  cheap there (one CG pass). A faster T (warm start from LP rounding, tighter delta from an incumbent, or a
  Lagrangian upper bound on T giving a lower bound on D) would be needed.
- q = 64 has only 23 instances and q = 80 only 2.
