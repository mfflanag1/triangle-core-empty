# Independent enclosure of the tree matching cost

The fine run gives an exact upper bound of **10135505296350745 / 4503599627370496** on the cost of the tenth rounded law. This is below 2.250535 and hence below 2.2519. The coarser grid gives **2534881595504131 / 1125899906842624**, below 2.251428. Lemmas 7.13 and 7.14 in §7.6 make either result an upper bound on the fixed-point cost used in the finite-market fractional upper bound.

## Reproduction

These runs used python-flint 0.9.0 with FLINT 3.6.0. From this directory, run

```sh
python -B cert_fractional.py
python -B cert_fractional.py --grid-bits 11 --precision 224
python -B cert_fractional_checks.py
```

The first two transcripts are `cert_fractional_run.txt` and `cert_fractional_coarse.txt`. The per-stage decimals round upward. The program converts its final Arb upper enclosure to an exact dyadic rational, then checks it against 22519/10000 using Python integers. Every run starts with the unit atom at zero and constructs its own ten rounded laws.

## Enclosures and arithmetic

1. On a fresh uniform mesh of width `h=2^-12` (coarse run `2^-11`), the program retains the first and successive strictly increasing lower bounds on the hazard between about 0.0000001 and 36, including the first one reaching 36. The fine run uses 176-bit Arb real balls, the coarse run 224-bit balls. Slopes are positive running maxima of upper derivative bounds, establishing the shape of the next law.
2. Write `y=-x`. On a rounded piece, the survival probability increases exponentially with rate equal to the piece's slope. Its integral from one endpoint to the other is the increase in the survival probability divided by that slope. Arb encloses the breakpoints, integrals and survival probabilities. Every comparison placing a positive grid sample on a piece must be decided by the balls. An overlap aborts the run. The atom at zero is accounted for separately. The integral of survival minus one over the support gives the mean of the current law.
3. The cell masses and normalized within-cell first moments are rounded downward, and the chord masses upward, to multiples of `2^-70`. Three FLINT integer polynomial products compute the convolution coefficients for the Lemma 7.15. The positive-part powers in the lemma are evaluated by exact prefix moments of orders zero through three. At grid index n, each prefix moment is the sum of the coefficients at j ≤ n times (n−j) to the requested power. Shifting n uses the binomial formula and the new coefficient. Integer floors and ceilings produce the one-sided hazard and derivative values used by the next law.
4. The tangent lower bound for the hazard integrates to an upper bound on the mean of its inverse at an exponential level. The code covers the initial segment by its length, integrates every exponential segment with Arb's exponential and `expm1`, and analytically bounds the entire last tail. Three times the sum of this upper bound and twice the enclosed current-law mean gives the final cost upper endpoint.

The calculation uses the analytic rounding, chain and lattice-enclosure lemmas in §7. The implementation, mesh, exact-polynomial convolutions and exact integer kernel evaluation are separate from the deposited NumPy/mpmath computation. No random samples, saved laws, BLAS operations, or data files enter the calculation. Independent analytic test cases check the prefix recurrence against direct sums, the cell formulas for a truncated exponential law, and both kernel inequalities against exact rational expectations over nine pairs of off-grid atoms.
