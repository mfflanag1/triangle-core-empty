"""Certify the fractional tree bound using Arb balls and FLINT integer polynomials.

Run with Python 3 and python-flint 0.9: python cert_fractional.py
The file is self-contained and uses no stored laws or original certificate routines.
"""

import argparse
from dataclasses import dataclass
from fractions import Fraction
from itertools import chain, repeat

import flint
from flint import arb, ctx, fmpz_poly

# The cost and hazard meshes use the same dyadic width. The other grid rounds
# masses and parameters for exact convolutions and exact stage-to-stage input.
GRID_BITS = 12
FRACTION_BITS = 70
UNIT = 1 << FRACTION_BITS
PRECISION = 176
LEVEL_START = (UNIT + 10**7 - 1) // 10**7
LEVEL_END = 36 * UNIT
LAST_STAGE = 10


def dyadic(integer, bits):
    return arb((integer, -bits))


def endpoint(value, bits, upper):
    """Integer bounding value * 2**bits, rounded in the requested direction."""
    scaled = value * (1 << bits)
    if upper:
        return int(scaled.upper().ceil().unique_fmpz())
    return int(scaled.lower().floor().unique_fmpz())


@dataclass
class MeshPoint:
    index: int        # t = index * 2^-GRID_BITS
    theta: int        # lower bound for g(t), in units of 2^-FRACTION_BITS
    derivative: int   # positive lower bound for g'(t), same units
    slope: int        # monotone upper bound for g'(t), same units


def cells(points, grid_bits=GRID_BITS):
    """Enclose cell masses, normalized first moments, chord weights and EX.

    Work in the reversed coordinate y=-x. The survival S(y)=P(X>-y)
    increases, with S(0) interpreted as the mass of the terminal atom. Each
    segment contributes (S(right)-S(left))/slope to the survival integral.
    """
    if not points:
        return [UNIT], [0], [UNIT], arb(0)

    # Successive breakpoints in y, moving away from the atom at zero.
    levels = [dyadic(p.theta, FRACTION_BITS) for p in reversed(points)] + [arb(0)]
    slopes = [dyadic(p.slope, FRACTION_BITS) for p in reversed(points)]
    breaks = [arb(0)]
    survival = [(-level).exp() for level in levels]
    primitive = [arb(0)]
    for i, s in enumerate(slopes):
        breaks.append(breaks[-1] + (levels[i] - levels[i + 1]) / s)
        primitive.append(primitive[-1] + (survival[i + 1] - survival[i]) / s)

    mass = [max(0, endpoint(survival[0], FRACTION_BITS, False))]
    moment = [0]
    chord = []
    previous_average = arb(0)  # J_0=0 at the atom in the chord formula
    previous_survival = survival[0]
    previous_primitive = arb(0)
    i = 0
    j = 0
    while i < len(slopes):
        j += 1
        y = dyadic(j, grid_bits)
        while i < len(slopes) and breaks[i + 1] <= y:
            i += 1
        if i == len(slopes):
            current_survival = arb(1)
            integral = primitive[i] + y - breaks[i]
        else:
            if not (y < breaks[i + 1] and breaks[i] <= y):
                raise ArithmeticError("Inconclusive comparison of lattice and breakpoints")
            current_survival = (-levels[i] + slopes[i] * (y - breaks[i])).exp()
            integral = primitive[i] + (current_survival - survival[i]) / slopes[i]

        average = (integral - previous_primitive) * (1 << grid_bits)
        mass.append(max(0, endpoint(current_survival - previous_survival, FRACTION_BITS, False)))
        moment.append(max(0, endpoint(average - previous_survival, FRACTION_BITS, False)))
        # The maximum with zero is safe for any upper enclosure as well.
        chord.append(max(0, endpoint(average - previous_average, FRACTION_BITS, True)))
        previous_average, previous_survival, previous_primitive = average, current_survival, integral

    chord.append(max(0, endpoint(1 - previous_average, FRACTION_BITS, True)))
    if sum(mass) > UNIT or sum(chord) < UNIT:
        raise ArithmeticError("Cell probabilities failed the mass checks")
    # E X = -integral_0^y (1-S(z)) dz. The upper endpoint y is beyond the support.
    return mass, moment, chord, integral - y


def power_sums(coefficients, degree):
    """Yield sum_j coefficients[j]*(n-j)_+^k at each n, k=0,...,degree.

    For k=0, sum over j<=n. Pascal's identity shifts n to n+1.
    All coefficients and all intermediate moments are exact integers.
    """
    sums = [0] * (degree + 1)
    for coefficient in chain(coefficients, repeat(0)):
        if degree == 3:
            sums[3] += 3 * sums[2] + 3 * sums[1] + sums[0]
        sums[2] += 2 * sums[1] + sums[0]
        sums[1] += sums[0]
        sums[0] += coefficient
        yield sums


def build_mesh(mass, moment, chord, grid_bits=GRID_BITS):
    p, d, nu = map(fmpz_poly, (mass, moment, chord))
    # d rounds D/h down. The polynomials are indexed by j=-k. FLINT
    # computes convolution exactly in Z[x], including its internal FFTs.
    pp = power_sums(map(int, p * p), 3)
    dp = power_sums(map(int, d * p), 2)
    uu = power_sums(map(int, nu * nu), 2)

    # Numerators of the tangents come from phi_3(u)=u_+^3/6 and
    # phi_2(u)=u_+^2/2. Coefficients each carry two mass factors.
    g_divisor = 6 << (FRACTION_BITS + 3 * grid_bits)
    derivative_divisor = 2 << (FRACTION_BITS + 2 * grid_bits)
    points = []
    last_level, last_slope = 0, 0
    for n in range((64 << grid_bits) + 1):
        m, first, upper = next(pp), next(dp), next(uu)
        theta = (m[3] + 6 * first[2]) // g_divisor
        if theta < LEVEL_START or theta <= last_level:
            continue
        derivative = (m[2] + 4 * first[1]) // derivative_divisor
        slope = (upper[2] + derivative_divisor - 1) // derivative_divisor
        if derivative <= 0 or derivative > slope:
            raise ArithmeticError("Invalid derivative enclosure")
        last_level, last_slope = theta, max(last_slope, slope)
        points.append(MeshPoint(n, theta, derivative, last_slope))
        if theta >= LEVEL_END:
            return points
    raise ArithmeticError("Insufficient t range")


def cost_bound(points, mean_x, grid_bits=GRID_BITS):
    """Integrate exp(-g) using a lower supporting line on each interval."""
    total = dyadic(points[0].index, grid_bits)
    for i, point in enumerate(points):
        theta = dyadic(point.theta, FRACTION_BITS)
        a = dyadic(point.derivative, FRACTION_BITS)
        area = (-theta).exp() / a
        if i + 1 < len(points):
            width = dyadic(points[i + 1].index - point.index, grid_bits)
            area *= -(-a * width).expm1()
        total += area
    return 3 * (total + 2 * mean_x)


def main(grid_bits, precision):
    ctx.prec = precision
    print(f"python-flint {flint.__version__}, FLINT {flint.__FLINT_VERSION__}, "
          f"precision {ctx.prec} bits, grid 2^-{grid_bits}, coefficients 2^-{FRACTION_BITS}", flush=True)
    points = []
    for stage in range(LAST_STAGE + 1):
        mass, moment, chord, mean_x = cells(points, grid_bits)
        points = build_mesh(mass, moment, chord, grid_bits)
        upper = endpoint(cost_bound(points, mean_x, grid_bits), 54, True)
        outward_six_digits = (upper * 10**6 + (1 << 54) - 1) >> 54
        print(f"Stage {stage:2d}: {len(mass)-1:6d} cells, {len(points):6d} mesh points, "
              f"cost <= {outward_six_digits // 10**6}.{outward_six_digits % 10**6:06d}", flush=True)
    result = Fraction(upper, 1 << 54)
    target = Fraction(22519, 10000)
    if result > target:
        raise ArithmeticError(f"Final upper endpoint {result} exceeds {target}")
    print(f"Exact upper endpoint {result} <= {target}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-bits", type=int, default=GRID_BITS)
    parser.add_argument("--precision", type=int, default=PRECISION)
    args = parser.parse_args()
    if args.grid_bits < 1 or args.precision < 2:
        parser.error("grid bits and precision must be positive, with precision at least two")
    main(args.grid_bits, args.precision)
