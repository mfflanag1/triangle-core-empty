"""Focused checks of the independent fractional enclosure on analytic laws."""

from fractions import Fraction

from flint import arb, ctx

from cert_fractional import (GRID_BITS, FRACTION_BITS, PRECISION, UNIT, MeshPoint,
                             build_mesh, cells, dyadic, power_sums)


def check_moment_recurrence():
    coeff = [3, 8, 0, 5, 2]
    for degree in (2, 3):
        for n, row in zip(range(14), power_sums(coeff, degree)):
            for k in range(degree + 1):
                expected = sum(a * (n - j)**k for j, a in enumerate(coeff) if j <= n)
                assert row[k] == expected


def check_survival_lattice():
    # X = min(E-5/4,0). Exact endpoints include a gridpoint at the terminal kink.
    points = [MeshPoint(1, 5 * UNIT // 4, UNIT, UNIT)]
    mass, moment, chord, mean = cells(points)
    a = arb(5) / 4
    analytic_mean = 1 - a - (-a).exp()
    assert (mean - analytic_mean).contains(0)
    h = dyadic(1, GRID_BITS)
    at_zero = (-a).exp()
    assert dyadic(mass[0], FRACTION_BITS) <= at_zero
    for j in (1, 2, len(mass) // 2, len(mass) - 1):
        left, right = at_zero * ((j - 1) * h).exp(), at_zero * (j * h).exp()
        p = right - left
        j_bar = p / h
        # With a unit exponential slope, average_s = (S(right)-S(left))/h.
        assert dyadic(mass[j], FRACTION_BITS) <= p
        assert dyadic(moment[j], FRACTION_BITS) <= j_bar - left
        assert dyadic(chord[j - 1], FRACTION_BITS) >= (
            j_bar if j == 1 else j_bar * (1 - (-h).exp()))
    assert dyadic(chord[-1], FRACTION_BITS) >= 1 - (-(-h).expm1()) / h


def check_kernel_enclosures():
    # Three off-grid/at-grid atoms, with exact dyadic locations and masses.
    # Compute the true expectations directly over the nine ordered pairs.
    atoms = [(0, Fraction(0), Fraction(1, 4)),
             (147, Fraction(1, 4), Fraction(1, 2)),
             (311, Fraction(3, 4), Fraction(1, 4))]
    n = 312
    mass, moment, chord = [0] * n, [0] * n, [0] * n
    for j, fraction, probability in atoms:
        mass[j] += int(UNIT * probability)
        moment[j] += int(UNIT * probability * fraction)
        chord[j] += int(UNIT * probability * (1 - fraction))
        if j:
            chord[j - 1] += int(UNIT * probability * fraction)
    mesh = build_mesh(mass, moment, chord)
    for point in [mesh[0], mesh[100], mesh[235], mesh[500], mesh[-1]]:
        g, derivative = Fraction(0), Fraction(0)
        for j, a, p in atoms:
            for k, b, q in atoms:
                u = max(0, point.index - j + a - k + b) / Fraction(1 << GRID_BITS)
                g += p * q * u**3 / 6
                derivative += p * q * u**2 / 2
        assert Fraction(point.theta, UNIT) <= g
        assert Fraction(point.derivative, UNIT) <= derivative <= Fraction(point.slope, UNIT)


if __name__ == "__main__":
    ctx.prec = PRECISION
    check_moment_recurrence()
    check_survival_lattice()
    check_kernel_enclosures()
    print("Passed: exact moment recurrence, cell enclosures for a truncated exponential, "
          "and kernel bounds for a three-atom off-grid law.")
