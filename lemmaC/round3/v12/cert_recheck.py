"""V12: re-run the hardened certificate in-process (no export written), compare with the exported rational, and
probe the pieces: exact parsing, midpoint exactness of the mesh, rexp / rlog against 60-digit references on
edge inputs, the gamma enclosure, and my own high-precision Psi_0 against the certified endpoint."""
import sys, json, math, random, time
from fractions import Fraction
LANE = "/Users/maxf/projects/economics-research/game-theory/coopgpt/average_case_protection_2026-11-09/lemmaC/round3"
sys.path.insert(0, LANE + "/final/certI")
import mpmath
from mpmath import iv, mp
import certify_rigorous as C
mp.dps = 60; random.seed(12)
W = LANE + "/firstmoment/w_K5_e10.json"
d = json.loads(open(W).read(), parse_float=str)
print("types after parse_float=str:", {k: type(v).__name__ for k, v in d.items() if k != "u"}, set(type(v).__name__ for v in d["u"]))
print("R('1.6712') encloses 2089/1250:", C.lo_frac(C.R("1.6712")) <= Fraction(2089, 1250) <= C.hi_frac(C.R("1.6712")),
      " width:", float(C.hi_frac(C.R("1.6712")) - C.lo_frac(C.R("1.6712"))))
# mesh: midpoints exact, endpoints exact
tf = 2.0 ** -27; a = iv.mpf(tf); bad = 0; n = 0
while tf < 900:
    tf *= 1.005; b = iv.mpf(tf); m = (a + b) / 2
    bad += not (C.lo_frac(m) == C.hi_frac(m) == (Fraction(float(a.a)) + Fraction(tf)) / 2) or C.lo_frac(b) != Fraction(tf); a = b; n += 1
print("mesh cells checked: %d, midpoints/endpoints not exact: %d" % (n, bad))
# rexp / rlog edge probes
bad = 0; cnt = 0
for _ in range(4000):
    x = random.choice([-random.random() * 10 ** random.uniform(-12, 3.3), random.uniform(-25, 3), -2.0 ** random.randint(-40, 10), 2.0 ** -random.randint(1, 60)])
    e = C.rexp(iv.mpf(x)); t = mp.exp(mp.mpf(x)); cnt += 1; bad += not (mp.mpf(e._mpi_[0]) <= t <= mp.mpf(e._mpi_[1]))
    hi = x + abs(x) * random.random(); e = C.rexp(iv.mpf([x, hi]))
    for y in (x, hi): cnt += 1; bad += not (mp.mpf(e._mpi_[0]) <= mp.exp(mp.mpf(y)) <= mp.mpf(e._mpi_[1]))
print("rexp: %d failures of %d" % (bad, cnt))
bad = 0; cnt = 0
for _ in range(4000):
    x = random.choice([1 + random.random() * 0.01, 10 ** random.uniform(-300, 300), 2.0 ** random.randint(-60, 60), random.random()])
    e = C.rlog(iv.mpf(x)); t = mp.log(mp.mpf(x)); cnt += 1; bad += not (mp.mpf(e._mpi_[0]) <= t <= mp.mpf(e._mpi_[1]))
print("rlog: %d failures of %d" % (bad, cnt))
g = C.euler_gamma(); print("gamma enclosure contains mp.euler:", mp.mpf(g._mpi_[0]) <= mp.euler <= mp.mpf(g._mpi_[1]))
t = time.time(); lines = []
out = C.certify(W, 0.005, log=lines.append)
exp_ = json.load(open(LANE + "/final/certI/psi0_export_w_K5_e10.json"))
print("re-run: %.0fs; Psi0_up = %.16f; equals exported rational: %s; logged lines identical to out_certify_rigorous_K5_e10.txt (except timing): %s"
      % (time.time() - t, float(out["Psi0_up"]), out["Psi0_up"] == Fraction(exp_["Psi0_up"]),
         [l for l in lines][1:] == open(LANE + "/final/certI/out_certify_rigorous_K5_e10.txt").read().splitlines()[1:len(lines)]))
mine = mp.mpf("-3.847516394113")
print("my tanh-sinh Psi_0 = %s <= certified endpoint: %s; slack %.3e; theta*2.302 + endpoint = %.10f" % (mine, mine <= mp.mpf(out["Psi0_up"].numerator) / out["Psi0_up"].denominator,
      float(mp.mpf(out["Psi0_up"].numerator) / out["Psi0_up"].denominator - mine), float(Fraction("1.6712") * Fraction("2.302") + out["Psi0_up"])))
