"""V6 independent interval enclosure of Psi_0 (THEOREM_I.md Def. 1.1), written from the definition.
Different enclosure principle from certify.py: ONLY monotonicity of psi(t)=E exp(-tW) and of exp(-t w_b) psi(t)
(Laplace transforms of nonnegative variables are decreasing); no convexity, no derivative, no midpoint.
Mesh t_j = 2^-30 * (1+2^-LG)^j (exact binary floats), up to t1 >= 1024.
  int_{t0}^{t1} psi dt/t   in [sum psi(t_j) ln(t_j/t_{j-1}),  sum psi(t_{j-1}) ln(t_j/t_{j-1})]
  int_0^{t1} f dt          in [sum f(t_j)(t_j - t_{j-1}) + t0 f(t0),  sum f(t_{j-1})(t_j - t_{j-1}) + t0 f(0)]
  Phi = -gamma - ln t0 + [0, t0*EW] - int_{t0}^{t1} psi dt/t - [0, psi(t1)/(t1 w_inf)]
  tails of rho_b, eps_0 beyond t1 in [0, psi(t1)]."""
import sys, json, time
from fractions import Fraction
from mpmath import iv
iv.prec = 64
F = "/Users/maxf/projects/economics-research/game-theory/coopgpt/average_case_protection_2026-11-09/lemmaC/round3/firstmoment/" + sys.argv[1]
LG = int(sys.argv[2]); d = json.load(open(F))
def Q(x):
    f = Fraction(str(x)); return iv.mpf(f.numerator) / iv.mpf(f.denominator)
B = int(d["B"]); eta = Q(d["eta"]); th = Q(d["theta"]); K = eta * B
u = [Q(v) for v in d["u"]]; uinf = Q(d["uinf"])
w = [iv.exp(v) for v in u]; winf = iv.exp(uinf)
lam = [eta ** 3 * (3 * b * b - 3 * b + 1) / 6 for b in range(1, B + 1)]
def prim(x):  # exp(-th x) * (x^2/(2 th) + x/th^2 + 1/th^3);  int_a^b e^{-th x} x^2/2 dx = prim(a) - prim(b)
    return iv.exp(-th * x) * (x * x / (2 * th) + x / (th * th) + 1 / (th * th * th))
kap = [prim(eta * (b - 1)) - prim(eta * b) for b in range(1, B + 1)]; J = prim(K)
EW = winf + sum(l * x for l, x in zip(lam, w))
def node(t):
    e = [iv.exp(-t * x) for x in w]
    s = t * winf
    for l, x in zip(lam, e): s += l * (1 - x)
    p = iv.exp(-s)
    return p, e
ZERO = iv.mpf(0)
t0f = 2.0 ** -30; r = 1.0 + 2.0 ** -LG
tprev = iv.mpf(t0f); pprev, eprev = node(tprev)
Ilo = ZERO; Ihi = ZERO
# head cell [0,t0] for the dt-integrals: f decreasing, f(0)=1
rlo = [tprev * x * pprev for x in eprev]; rhi = [tprev for _ in range(B)]
elo = tprev * pprev; ehi = tprev
tf = t0f; n = 0; tic = time.time()
while tf < 1024.0:
    tf = tf * r; t = iv.mpf(tf); p, e = node(t)
    lg = iv.log(t / tprev); h = t - tprev
    Ilo += p * lg; Ihi += pprev * lg
    hp = h * p; hpp = h * pprev
    for i in range(B):
        rlo[i] += hp * e[i]; rhi[i] += hpp * eprev[i]
    elo += hp; ehi += hpp
    tprev, pprev, eprev = t, p, e; n += 1
t1 = tprev; tail = pprev.b
Phi_hi = (-iv.euler - iv.log(iv.mpf(t0f)) + EW * t0f - Ilo).b
Phi_lo = (-iv.euler - iv.log(iv.mpf(t0f)) - Ihi - pprev / (t1 * winf)).a
tot = 3 * iv.mpf(Phi_hi); nlo = ZERO; nhi = ZERO; wr = 0.0
rho1 = None
for i in range(B):
    ra = (w[i] * rlo[i]).a; rb = (w[i] * rhi[i] + tail).b
    assert ra > 0 and rb < 1
    if i == 0: rho1 = (ra, rb)
    C = iv.exp(-3 * u[i]) * kap[i] / lam[i]
    va = iv.log((1 - iv.mpf(ra)) ** 3 + C * iv.mpf(ra) ** 3).b; vb = iv.log((1 - iv.mpf(rb)) ** 3 + C * iv.mpf(rb) ** 3).b
    v = va if va > vb else vb          # convex in r (second derivative 6(1-r)+6Cr > 0): max at an endpoint
    tot += lam[i] * iv.mpf(v)
    nlo += lam[i] * iv.mpf(ra); nhi += lam[i] * iv.mpf(rb); wr = max(wr, float(((rb - ra) / ra).b))
ea = (winf * elo).a; eb = (winf * ehi + tail).b
Cinf = iv.exp(-3 * uinf) * J
ca = (-3 * iv.mpf(ea) + Cinf * iv.mpf(ea) ** 3).b; cb = (-3 * iv.mpf(eb) + Cinf * iv.mpf(eb) ** 3).b
tot += iv.mpf(ca if ca > cb else cb)
P = tot.b
print("file=%s ratio=1+2^-%d cells=%d t1=%s time=%.0fs" % (sys.argv[1], LG, n, iv.nstr(t1, 8), time.time() - tic))
print("Phi in [%s, %s]" % (iv.nstr(iv.mpf(Phi_lo), 12), iv.nstr(iv.mpf(Phi_hi), 12)))
print("eps0 in [%s, %s]; rho_1 in [%s, %s]; worst rel width rho %.2e" % (iv.nstr(iv.mpf(ea), 10), iv.nstr(iv.mpf(eb), 10), iv.nstr(iv.mpf(rho1[0]), 10), iv.nstr(iv.mpf(rho1[1]), 10), wr))
print("norm in [%s, %s]" % (iv.nstr(iv.mpf((nlo + ea).a), 10), iv.nstr(iv.mpf((nhi + eb).b), 10)))
print("V6 UPPER BOUND Psi0 <= %s" % iv.nstr(iv.mpf(P), 12))
print("V6 -Psi0/theta >= %s" % iv.nstr(iv.mpf((-iv.mpf(P) / th).a), 12))
for bt in ("2.30", "2.301", "2.302"):
    v = (th * Q(bt) + iv.mpf(P)).b
    print("  beta=%s: theta*beta+Psi0 <= %s %s" % (bt, iv.nstr(iv.mpf(v), 8), "NEGATIVE" if v < 0 else "not certified"))
