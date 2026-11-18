"""Independently enclose the T1 integer first moment and market transfer exponent.

Run with Python 3 and python-flint. The only external input is the paper's vector u.
See cert_integer.md for the quadrature and endpoint proofs.
"""

import argparse
import json
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

import flint
from flint import arb, ctx, fmpq


def rational(value):
    q = Fraction(value)
    return arb(fmpq(q.numerator, q.denominator))


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def widen(value, error):
    """An enclosing ball after an additive error of at most error."""
    return (value - error).union(value + error)


def load_witness():
    path = Path(__file__).resolve().parents[1] / "round3/firstmoment/w_K5_e10.json"
    data = json.loads(path.read_text(), parse_float=Decimal)
    require((data["K"], data["B"], data["eta"], data["theta"]) ==
            (5, 50, Decimal("0.1"), Decimal("1.6712")), "unexpected witness parameters")
    require(len(data["u"]) == data["B"], "unexpected number of coordinates")
    return [rational(u).exp() for u in data["u"]], rational(data["uinf"]).exp()


def expectations(weights, base, rates, panels):
    """Integrate on dyadic intervals; bound the complementary head and tail."""
    def laplace(t):
        decays = [(-w * t).exp() for w in weights]
        log_psi = -base * t + sum(lam * (decay - 1) for lam, decay in zip(rates, decays))
        return log_psi.exp(), decays

    start = arb(1) / 2**32
    stop = arb(2**13)
    mean_w = base + sum(lam * w for lam, w in zip(rates, weights))
    mean_log = arb(0)
    means = [arb(0) for _ in range(len(weights) + 1)]
    all_weights = weights + [base]
    log_error = start * (1 + mean_w)
    errors = [start * w for w in all_weights]
    # Each interval is [a,2a]. Its spacing is a/panels.
    for k in range(-32, 13):
        a = arb(2)**k
        step = a / panels
        local_log = arb(0)
        local_means = [arb(0) for _ in all_weights]
        for j in range(panels + 1):
            t = a + j * step
            coef = 1 if j in (0, panels) else (4 if j % 2 else 2)
            psi, decays = laplace(t)
            local_log += coef * ((-t).exp() - psi) / t
            # Multiply by each prefactor w once per interval.
            weighted_psi = coef * psi
            for i, decay in enumerate(decays):
                local_means[i] += weighted_psi * decay
            local_means[-1] += weighted_psi
        mean_log += step * local_log / 3
        for i, w in enumerate(all_weights):
            means[i] += step * w * local_means[i] / 3
        # Composite Simpson error <= 81*a*f(a/2)/(180*panels**4)
        # for each Laplace transform of a positive measure.
        half = a / 2
        psi, decays = laplace(half)
        scale = 81 * a / (180 * panels**4)
        log_error += scale * ((-half).exp() + psi) / half
        for i, (w, decay) in enumerate(zip(weights, decays)):
            errors[i] += scale * w * psi * decay
        errors[-1] += scale * base * psi

    log_error += ((-stop).exp() + (-base * stop).exp() / base) / stop
    for i, w in enumerate(weights):
        errors[i] += w * (-(w + base) * stop).exp() / (w + base)
    errors[-1] += (-base * stop).exp()
    enclosed = [widen(mean, err) for mean, err in zip(means, errors)]
    return widen(mean_log, log_error), enclosed[:-1], enclosed[-1]


def run(bits, panels):
    require(bits >= 64 and panels > 0 and panels % 2 == 0,
            "use at least 64 bits and an even positive number of panels")
    ctx.prec = bits
    weights, base = load_witness()
    rates = [arb(3*b*b - 3*b + 1) / 6000 for b in range(1, 51)]
    theta = rational("1.6712")
    beta = rational("2.302")
    K = 5
    phi, rho, eps = expectations(weights, base, rates, panels)
    require((eps > 0) and (eps < 1) and all(0 < r < 1 for r in rho),
            "invalid star probabilities")
    require((eps + sum(lam * r for lam, r in zip(rates, rho))).contains(1),
            "star mass identity failed")

    def tail_tilt(x):
        return (-theta * x).exp() * (1 + theta*x + (theta*x)**2 / 2) / theta**3

    tails = [tail_tilt(arb(i) / 10) for i in range(51)]
    high = tails[-1]
    psi0 = 3 * phi - 3 * eps + high * (eps / base)**3
    log_magnitudes = [abs(high.log())]
    for i, (lam, w, r) in enumerate(zip(rates, weights, rho)):
        mean_tilt = (tails[i] - tails[i+1]) / lam
        log_magnitudes.append(abs(mean_tilt.log()))
        psi0 += lam * ((1 - r)**3 + mean_tilt * (r/w)**3).log()
    c = 17
    delta1, delta2, tau = (rational(Fraction(1, n)) for n in (10**6, 50000, 4000000))
    box = (8*c + 4 + 2*theta*K) / 3
    require(all(ell < theta*K for ell in log_magnitudes), "tilt exceeds box hypothesis")
    load_bound = rates[-1] + delta1 + 2*tau
    require(load_bound > 1 and c >= 1 + load_bound.log(), "load hypothesis failed")
    error = 3 * (delta2 + 2*tau*(box + arb(2).log())) + 4*c*(delta1 + 2*tau)
    error += (6*len(weights) + 2) * arb(-1-c).exp()
    iid_exponent = theta * beta + psi0
    exponent = iid_exponent + error
    require(exponent < 0, "corrected first-moment exponent is not strictly negative")

    print(f"python-flint {flint.__version__}, FLINT {flint.__FLINT_VERSION__}, "
          f"{bits} bits, {panels} Simpson panels per octave")
    print("Interval endpoints below are outward-rounded to multiples of 10^-10.")
    for label, value in [("Psi0", psi0), ("unperturbed exponent", iid_exponent),
                         ("transfer error", error), ("market exponent at beta=2.302", exponent),
                         ("fixed-witness unperturbed exponent at beta=2.30225",
                          theta * rational("2.30225") + psi0)]:
        lower = int((value.lower() * 10**10).floor().unique_fmpz())
        upper = int((value.upper() * 10**10).ceil().unique_fmpz())
        print(f"{label}: [{Decimal(lower).scaleb(-10):f}, {Decimal(upper).scaleb(-10):f}]")
    print(f"Exact dyadic market exponent upper endpoint: {exponent.upper().fmpq()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bits", type=int, default=128)
    parser.add_argument("--panels", type=int, default=256)
    args = parser.parse_args()
    run(args.bits, args.panels)
