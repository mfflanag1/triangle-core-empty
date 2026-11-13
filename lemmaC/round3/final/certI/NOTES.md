# Theorem I certificate, hardened (Fix C)

Written 2026-11-12. Scope: the interval certificate of `firstmoment/THEOREM_I.md` section 5
and the error budget of `interface/THEOREM_I_PA.md` section 2. The mathematics is unchanged.
The originals (`firstmoment/certify.py`, `interface/constants.py`) were not edited.

## 1. Verdict

No mathematical error was found. Every number stated in the two theorem files that is used
downstream is confirmed by the hardened scripts:

| Quantity (K = 5, B = 50, theta = 1.6712, weights `w_K5_e10.json`) | Original log | Hardened, certified |
|---|---|---|
| Psi_0 upper bound | −3.84751328323 | **Psi_0 ≤ −3.84751328322** (and ≤ −3.84751328 at 8 decimals) |
| −Psi_0/theta lower bound | 2.30224586119 | **≥ 2.30224586119** |
| theta·beta + Psi_0, beta = 2.30 | −0.0037532832 | ≤ −0.0037532832 |
| beta = 2.302 | −0.00041088323 | ≤ −0.0004108832 |
| beta = 2.3022 | (−0.00007 in the text) | ≤ −0.0000766432 |
| beta = 2.30224 | not printed | ≤ −0.0000097952 |
| Perturbed exponent theta·beta + Psi_0 + Err, beta = 2.30 | −0.00059298308 | **≤ −0.0005929830** (Err ≤ 0.0031602970) |
| beta = 2.302 | −0.00016488483 | **≤ −0.0001648848** (Err ≤ 0.0002459952) |

The perturbed exponents use the 8-decimal constant Psi_0 ≤ −3.84751328, as the original did.
Second certificate (K = 5, B = 20, theta = 1.6671, `w_K5_e25.json`): Psi_0 ≤ −3.83621828482,
−Psi_0/theta ≥ 2.30113267639 (original log: −3.83621828482 and 2.3011326764).

Exact internal upper endpoints of the Psi_0 enclosure (binary rationals):

- original script: −581419768675579087666533/2^77 = −3.8475132832280412…
- hardened script: −1162839537350571633950917/2^78 = −3.8475132832261005…

The hardened bound is weaker by 1.94·10⁻¹². All of that is the enclosure of Euler's constant
(section 3), whose lower end is 6.5·10⁻¹³ below the true value; Phi enters with weight 3.
The replacement of exp and log moves nothing at this scale.

## 2. Weaknesses, reproduced

`test_defect_and_replacements.py`, log `out_test_defect_and_replacements.txt` (mpmath 1.3.0,
pure-Python backend, 80 bits, which is the precision of both scripts; arguments drawn from the
ranges the scripts use; reference values at 600 bits).

- **`iv.exp` is not an enclosure.** Upper endpoint below the true value for 21 of 100,000
  arguments, relative shortfall up to 2.0·10⁻²⁸ (2^−80 = 8.3·10⁻²⁵). One case is proved
  without any reference to mpmath's own functions: x = −0x1.31c2c4b9c3ff0p−2, where the
  200-bit verified exponential has lower bound 0.741860139211600604812658140023 and `iv.exp`
  returns upper endpoint 0.741860139211600604812658139927.
- **`iv.log` is not an enclosure.** Upper endpoint below the true value for 1 of 1,000,000
  arguments (relative shortfall 3.0·10⁻³⁰), proved the same way at x = 0x1.5a6695e45d861p+0.
  The lower endpoints never failed in these samples.
- Cause (read in `libmp/libmpi.py`): `mpi_exp` and `mpi_log` call `mpf_exp` / `mpf_log` with
  directed rounding. Those compute an approximation with a few guard bits and round *that*
  in the requested direction; when the approximation lies on the wrong side of a
  representable number the result is on the wrong side too.
- `iv.euler` contained the reference value at 80 bits, but its routine has the same design
  (a fixed-point value assumed to be rounded down), so it is replaced as well.
- **Printing.** `iv.nstr` rounds to nearest. On 20,000 random endpoints the printed
  "upper bound" was below the endpoint in 10,055 cases. A live instance: the original log
  line "Psi0 <= −3.84751328323" is 2.0·10⁻¹² below the script's own endpoint
  −3.8475132832280. The 8-decimal value −3.84751328 that `constants.py` copied is above the
  endpoint by 3.2·10⁻⁹, so nothing downstream was affected.

The size of the defect (below 10⁻²⁷ relative per call) is far too small to move a margin of
4·10⁻⁴; the point of the fix is that the certificate no longer rests on that estimate.

## 3. What changed

`certify_rigorous.py` (also a module: `certify(wfile, delta)` returns the exact endpoint).

1. Trusted base reduced to: `mpmath.iv` addition, subtraction, multiplication, division,
   comparisons, and conversion of ints and floats (outward rounded); Python integers and
   `fractions.Fraction`. Not used: `iv.exp`, `iv.log`, `iv.euler`, interval `**`, `iv.nstr`.
2. `rexp(x)`: r = x/2^k with |r| ≤ 2^−10 (asserted on the interval endpoints);
   exp(r) ∈ Σ_{j<12} r^j/j! + [−ρ, ρ], ρ = 2^−119/12! ≥ |r|^12/12!·e^|r| (Lagrange
   remainder); then k interval squarings. Adapted from `cert/cert.py` after checking the
   remainder bound and the reduction.
3. `rlog(x)`: for a point p > 0 and any real y_0 (the double `math.log(p)`, a guess that is
   not trusted), ln p = y_0 + ln(1 + d) with d = p/exp(y_0) − 1 enclosed by `rexp`, and
   d − d² ≤ ln(1 + d) ≤ d for |d| ≤ 1/2 (|d| ≤ 2^−20 is asserted). An interval is treated
   endpoint by endpoint, since ln is increasing.
   Proof of the inequality: ln(1 + d) ≤ d always. For 0 ≤ d: ln(1 + d) ≥ d − d²/2. For
   −1/2 ≤ d < 0, with s = −d: −ln(1 − s) = s + s²/2 + s³/3 + … ≤ s + (s²/2)/(1 − s) ≤ s + s².
4. Euler's constant: γ = H_n − ln(n + 1/2) − S_n with 0 ≤ S_n ≤ 2/(45n²), n = 2^16, H_n summed
   in interval arithmetic. Proof: c_k = H_k − ln(k + 1/2) → γ, and with y = 1/(2k),
   c_{k−1} − c_k = ln((1 + y)/(1 − y)) − 2y = 2·Σ_{j≥1} y^(2j+1)/(2j + 1), which lies between
   0 and (2y³/3)/(1 − y²) ≤ (16/15)/(12k³) for k ≥ 2; and Σ_{k>n} k^−3 ≤ 1/(2n²).
   Result: γ ∈ [0.57721566490088, 0.57721566491124].
5. Powers are written as products (arguments are nonnegative where it matters).
6. Weights are parsed as exact decimals (no float round trip).
7. Output: every printed bound is produced by `dec_le` / `dec_ge`, which take the exact binary
   endpoint as a `Fraction`, apply floor or ceiling at the requested number of decimals, and
   assert the inequality. The run also writes `psi0_export_<weights>.json` with the exact
   endpoint as a rational, the outward-rounded decimals (asserted), the parameters and the
   SHA-256 of the weights file.
8. Two added assertions: the mesh is increasing; Σ λ_b ρ_b + ε_0 encloses 1.

`constants_rigorous.py`.

1. Same `rexp` / `rlog`; no `iv.exp`, `iv.log`, `**`, `iv.nstr`.
2. The constant −3.84751328 is asserted, in exact rational arithmetic, to be at least the exact
   endpoint. By default the endpoint is read from `psi0_export_w_K5_e10.json` (weights hash and
   K, B, eta, theta checked). With `--fresh` it is recomputed in the same process. Both were
   run and give identical output.
3. The hypothesis c' ≥ 1 + ln^+(λ_max + δ_1 + 2τ) is asserted including the case ln^+ = 0, and
   |ln J| < theta·K is asserted on both sides.

Evidence for the replacements (`test_defect_and_replacements.py`, 80 bits, 600-bit reference):
`rexp` 0 failures in 100,000 points and 60,000 checks on wide intervals (relative width at
most 1.5·10⁻¹⁷); `rlog` the same counts, 0 failures (relative width at most 2.4·10⁻¹⁷); the
four interval operations 0 failures in 468,775 containment checks against exact rationals
(points and wide intervals of either sign); `dec_le` / `dec_ge` 0 failures in 20,000.

Cost: 75 s for the main certificate (original: 32 s).

## 4. Item 3: corrected statements (the two theorem files were not edited)

The equation with the prefactor is (3.1) of `interface/THEOREM_I_PA.md`. Equation (3.1) of
`firstmoment/THEOREM_I.md` is the choice-function bound and has no prefactor; the matching
place there is the last paragraph of its section 3.2.

**(3.1) of THEOREM_I_PA.md, corrected.** For every set 𝒯 of degree arrays and q^(2/3) ≥ K,

  E'[Z_β; D(H') ∈ 𝒯] ≤ c_q·(q + 1)^(B+2) · sup over D ∈ 𝒯 and profiles μ of
            [ Π_b (n_b(D) + 1)² · exp(q·Ξ(μ; D)) ],
  c_q = e²·exp(K⁶/(36q(1 − p)))·(1 − p)^(−q) ≤ e²·exp(K^(9/2)/30 + K^(3/2)/5).

The degree-dependent factor Π_b (n_b + 1)² is inside the supremum. If
N(𝒯) := sup over D ∈ 𝒯 of max_b ν_b(D) is finite, the right side is at most
C_0·(1 + N(𝒯))^(2B)·q^(3B+2)·exp(q·sup Ξ), with C_0 = e²·2^(B+2)·exp(K^(9/2)/30 + K^(3/2)/5)
depending only on K and B. (Bounds used: p ≤ 1/6 and q ≥ K^(3/2) when q^(2/3) ≥ K;
q + 1 ≤ 2q; n_b + 1 ≤ q·(1 + ν_b).)

**(I-a\*)(iii), corrected.** For q^(2/3) ≥ K,
E'[Z_β; H' ∈ 𝒢_q^(2τ)] ≤ C·q^(3B+2)·exp(−(m − Err)·q) with
C = C_0(K, B)·(1 + λ_max + δ_1 + 2τ)^(2B). So C depends on K, B and on δ_1 + 2τ, and not on
q. Under the standing hypothesis c' ≥ 1 + ln^+(λ_max + δ_1 + 2τ) one has
C ≤ C_0(K, B)·(1 + e^(c'−1))^(2B); if δ_1 + 2τ ≤ 1 (true for both certified rows, where it is at
most 2·10⁻⁵) then C ≤ C_0(K, B)·(2 + λ_max)^(2B), which depends only on K and B. All of these
parameters are fixed before q → ∞, so sections 5 and 6 of THEOREM_I_PA.md are unaffected.

**THEOREM_I.md section 3.2, last sentence, corrected.** The polynomial factor
e²·(q + 1)²·Π_b (n_b + 1)² depends on the degree array and stays inside the maximum over D in
(2.2). On 𝒯_δ it is at most e²·(q + 1)²·(1 + (λ_max + δ)·q)^(2B) by 3.3(i), uniformly in D, so
together with the (q + 1)^B profiles the total factor is exp(o(q)) uniformly over 𝒯_δ. That is
what Lemma 3.2 uses.

## 5. Other observations

- Three displayed numbers in THEOREM_I.md section 5.3 are rounded to nearest and overstate
  the computed enclosure by at most one unit in the last digit: the lower end of ε_0
  (0.02081077389; outward value 0.02081077388), the upper end of ρ_1 (0.7565004437; outward
  value 0.7565004438), and the 12-digit Psi_0 line of the log. None is used downstream.
  Outward-rounded values: Phi ∈ [−0.522226313363, −0.522223198517],
  ε_0 ∈ [0.02081077388, 0.02081090331], ρ_1 ∈ [0.7564957391, 0.7565004438],
  Σ λ_b ρ_b + ε_0 ∈ [0.99999792704, 1.00000414596], ln J ∈ [−6.106378003, −6.106378002],
  E W ≤ 0.72310583.
- The table of THEOREM_I_PA.md section 2 is confirmed as printed (m ≥ 0.0037532 and
  0.0004108; Err ≤ 0.0031604 and 0.0002460; m − Err ≥ 0.00059 and 0.00016; A ≤ 44.2374 and
  52.2374).
- What is still trusted: the four interval operations and comparisons of mpmath 1.3.0, Python
  integer and rational arithmetic, and the script logic (one author for the hardening, on top
  of the one author of the original). The interval operations passed the exact-rational test;
  that is evidence, not a proof of the library.
- Not done: an independent re-implementation in another interval library; a second read of
  the hardened scripts.

## 6. Files

- `certify_rigorous.py`, `out_certify_rigorous_K5_e10.txt`, `out_certify_rigorous_K5_e25.txt`
- `psi0_export_w_K5_e10.json`, `psi0_export_w_K5_e25.json` (exact hand-off)
- `constants_rigorous.py`, `out_constants_rigorous.txt` (export mode),
  `out_constants_rigorous_fresh.txt` (`--fresh`)
- `test_defect_and_replacements.py`, `out_test_defect_and_replacements.txt`

Re-run, from this folder, with P = /Users/maxf/miniforge3/bin/python3:

    $P certify_rigorous.py ../../firstmoment/w_K5_e10.json 0.005     # 75 s
    $P certify_rigorous.py ../../firstmoment/w_K5_e25.json 0.005     # 30 s
    $P constants_rigorous.py                                         # 1 s
    $P constants_rigorous.py --fresh                                 # 80 s
    $P test_defect_and_replacements.py 100000                        # about 4 min
