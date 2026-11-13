"""evidence only: rexp against 400-bit mp.exp; iv +,-,*,/ against exact rationals."""
import sys, random, time
sys.argv = ['cert.py', '4', '10', '10', '1']
sys.path.insert(0, '/Users/maxf/projects/economics-research/game-theory/coopgpt/average_case_protection_2026-11-09/lemmaC/round3/cert')
import cert
from cert import iv, rexp
from mpmath import mp
from fractions import Fraction
mp.prec = 400
random.seed(5); bad = 0; wmax = 0; t0 = time.time(); n = 40000
for i in range(n):
    x = [-random.random() * 40, -random.random() * 1e-3, -random.random() * 1e-9, -random.random()][i % 4]
    e = rexp(iv.mpf(x)); t = mp.exp(mp.mpf(x))
    a = mp.mpf(e._mpi_[0]); b = mp.mpf(e._mpi_[1])
    bad += not (a <= t <= b); wmax = max(wmax, float((b - a) / t))
print("rexp: %d failures of %d; max relative width %.2e; %.0f us per call" % (bad, n, wmax, 1e6 * (time.time() - t0) / n))
def fr(raw):
    s, m, e, bc = raw; v = Fraction(int(m)) * (Fraction(2) ** e); return -v if s else v
bad = 0
for i in range(40000):
    p = random.random() * 10 ** random.randint(-8, 3) * random.choice([1, -1]); q = random.random() * 10 ** random.randint(-8, 3) + 1e-9
    P, Q = Fraction(p), Fraction(q)
    for z, exact in ((iv.mpf(p) + iv.mpf(q), P + Q), (iv.mpf(p) - iv.mpf(q), P - Q), (iv.mpf(p) * iv.mpf(q), P * Q), (iv.mpf(p) / iv.mpf(q), P / Q)):
        bad += not (fr(z._mpi_[0]) <= exact <= fr(z._mpi_[1]))
print("iv + - * /: %d failures of 160000 against exact rational arithmetic" % bad)
