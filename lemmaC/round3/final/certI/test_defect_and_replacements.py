"""Evidence script (no part of the proof).
 A. mpmath 1.3.0 iv.exp / iv.log are not rigorous enclosures at the precision used by certify.py and
    constants.py (80 bits; 120 bits also shown): endpoint violations against 600-bit mp values, and one
    violation of each kind confirmed rigorously (200-bit rexp / rlog bound vs the 80-bit iv endpoint).
 B. The replacements rexp, rlog (certify_rigorous.py) enclose the 600-bit reference on every random argument.
 C. iv + - * / against exact rational arithmetic (points and wide intervals).
 D. The enclosure of Euler's constant; iv.euler for comparison.
 E. dec_le / dec_ge against exact rationals; how often the old printing routine (iv.nstr) rounds the wrong way.
Argument ranges are those of the two scripts: exp at -t*w in -[1e-8, 1800], at u in [-7, 1], at -theta*x in
[-9, 0]; log at mesh ratios 1.005 / 1.01, at values in (0, 2].
"""
import sys, os, math, random, time
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath
from mpmath import iv, mp
import certify_rigorous as C
from certify_rigorous import rexp, rlog, lo_frac, hi_frac, dec_le, dec_ge

N = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
mp.prec = 600
random.seed(20261112)
print("mpmath %s, backend %s, N = %d per family" % (mpmath.__version__, mpmath.libmp.BACKEND, N))

EXP_GEN = [lambda: -random.random() * 10 ** random.uniform(-8, 3.26), lambda: random.uniform(-7, 1), lambda: -random.uniform(0, 9)]
LOG_GEN = [lambda: 1 + random.random() * 0.01, lambda: random.random() * 2 + 1e-9, lambda: 10 ** random.uniform(-9, 0)]
def ends(e): return mp.mpf(e._mpi_[0]), mp.mpf(e._mpi_[1])


def set_prec(p):
    iv.prec = p
    C._C[:] = [iv.mpf(1) / math.factorial(j) for j in range(12)]     # 1/j! at the new precision; _REM stays a valid bound


def rlog_refined(p):
    """rlog with a refined guess (for the 200-bit confirmation only): any real y is a valid guess."""
    y = iv.mpf(math.log(float(p)))
    for _ in range(2):
        y = (y + (p / rexp(y) - 1)).mid
    d = p / rexp(y) - 1
    return iv.mpf([(y + (d.a - d.a * d.a)).a, (y + d.b).b])


print("\nA. defect of iv.exp / iv.log")
for prec in (80, 120):
    set_prec(prec)
    for name, f, g, gens in (("exp", iv.exp, mp.exp, EXP_GEN), ("log", iv.log, mp.log, LOG_GEN)):
        lo_bad = hi_bad = 0; worst = 0.0; first = None; t0 = time.time(); n = N * (10 if name == "log" else 1)   # the log defect is rarer
        for i in range(n):
            x = gens[i % 3](); a, b = ends(f(iv.mpf(x))); t = g(mp.mpf(x))
            if a > t: lo_bad += 1; worst = max(worst, float((a - t) / abs(t))); first = first or (x, "lower")
            if b < t: hi_bad += 1; worst = max(worst, float((t - b) / abs(t))); first = first or (x, "upper")
        print("  prec %3d iv.%s: lower endpoint above the true value %d times, upper endpoint below it %d times, of %d; "
              "worst relative violation %.2e (2^-prec = %.1e)  [%.0fs]" % (prec, name, lo_bad, hi_bad, n, worst, 2.0 ** -prec, time.time() - t0), flush=True)
        if first and prec == 80:
            x, side = first
            e80 = f(iv.mpf(x))
            set_prec(200)
            ref = rexp(iv.mpf(x)) if name == "exp" else rlog_refined(iv.mpf(x))       # rigorous 200-bit enclosure
            proved = (ref.a > e80.b) if side == "upper" else (ref.b < e80.a)
            print("     witness x = %r (= %s exactly): iv.%s(x) %s endpoint = %s; rigorous 200-bit enclosure of the true value: [%s, %s]; violation proved: %s"
                  % (x, x.hex(), name, side, mp.nstr(ends(e80)[1 if side == "upper" else 0], 32), mp.nstr(ends(ref)[0], 32), mp.nstr(ends(ref)[1], 32), proved is True))
            set_prec(80)
set_prec(80)

print("\nB. replacements at prec 80 (as used by certify_rigorous.py and constants_rigorous.py)")
for name, f, g, gens in (("rexp", rexp, mp.exp, EXP_GEN), ("rlog", rlog, mp.log, LOG_GEN)):
    bad = 0; wmax = 0.0; t0 = time.time()
    for i in range(N):
        x = gens[i % 3](); a, b = ends(f(iv.mpf(x))); t = g(mp.mpf(x))
        bad += not (a <= t <= b); wmax = max(wmax, float((b - a) / abs(t)))
    print("  %s on points: %d failures of %d; largest relative width %.2e; %.0f us per call" % (name, bad, N, wmax, 1e6 * (time.time() - t0) / N), flush=True)
    bad = 0
    for i in range(N // 5):
        x = gens[i % 3](); hi = x + abs(x) * 10 ** random.uniform(-12, -1); a, b = ends(f(iv.mpf([x, hi])))
        for y in (x, hi, min(max(random.uniform(x, hi), x), hi)):          # doubles inside [x, hi]
            t = g(mp.mpf(y)); bad += not (a <= t <= b)
    print("  %s on intervals of relative width 1e-12..1e-1: %d failures of %d point checks" % (name, bad, 3 * (N // 5)), flush=True)

print("\nC. iv + - * / at prec 80 against exact rationals (points and wide intervals of either sign; corners and one interior pair)")
def rnd(): return random.random() * 10 ** random.randint(-8, 3) * random.choice([1, -1])
bad = cnt = 0
for i in range(N // 4):
    x1 = x2 = rnd(); y1 = y2 = rnd()
    if i % 2: x1, x2 = sorted((x1, rnd()))
    if i % 4 >= 2: y1, y2 = sorted((y1, rnd()))
    X, Y = iv.mpf([x1, x2]), iv.mpf([y1, y2]); div = y1 > 0 or y2 < 0
    res = [X + Y, X - Y, X * Y] + ([X / Y] if div else [])
    for p, q in ((x1, y1), (x1, y2), (x2, y1), (x2, y2), (min(max(random.uniform(x1, x2), x1), x2), min(max(random.uniform(y1, y2), y1), y2))):
        P, Q = Fraction(p), Fraction(q)
        for z, t in zip(res, [P + Q, P - Q, P * Q] + ([P / Q] if div else [])):
            cnt += 1; bad += not (lo_frac(z) <= t <= hi_frac(z))
print("  %d failures of %d containment checks" % (bad, cnt))

print("\nD. Euler's constant")
g = C.euler_gamma(); a, b = ends(g); t = mp.euler + 0
print("  rigorous enclosure [%s, %s] contains the 600-bit value: %s; width %.2e" % (dec_le(g, 14), dec_ge(g, 14), a <= t <= b, float(b - a)))
a, b = ends(iv.euler + 0)
print("  iv.euler (80 bits) contains it: %s  (not used by the hardened scripts; its fixed-point routine carries no proof of directed rounding)" % (a <= t <= b))

print("\nE. decimal hand-off")
bad = wrong_hi = wrong_lo = 0; M = N // 5
for i in range(M):
    x = iv.mpf(random.uniform(-4, 4)) / 3; nd = random.randint(6, 12)
    lo_s, hi_s = dec_le(x, nd), dec_ge(x, nd)
    bad += not (Fraction(lo_s) <= lo_frac(x) and Fraction(hi_s) >= hi_frac(x) and Fraction(hi_s) - Fraction(lo_s) <= Fraction(2, 10 ** nd))
    pr = lambda e: Fraction(iv.nstr(e, 12).strip("[]").split(",")[0])     # the original scripts print "[v, v]" for a point
    wrong_hi += pr(x.b) < hi_frac(x); wrong_lo += pr(x.a) > lo_frac(x)
print("  dec_le / dec_ge: %d failures of %d (exact rational comparison)" % (bad, M))
print("  old routine iv.nstr(endpoint, 12): printed upper bound below the upper endpoint in %d of %d cases; printed lower bound above the lower endpoint in %d" % (wrong_hi, M, wrong_lo))
