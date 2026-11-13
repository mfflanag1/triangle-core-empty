"""Interval-arithmetic certificate (mpmath.iv) for THEOREM_I.md, section 5.
Model: bins I_b = ((b-1)eta, b*eta], b = 1..B, K = B*eta, plus the high option.
Input: w file with theta, u_b = ln w_b, uinf = ln w_inf (decimal strings read as exact rationals).
Output: rigorous upper bound on
  Psi0 = 3 Phi + sum_b lam_b ln((1-rho_b)^3 + C_b rho_b^3) - 3 eps + C_inf eps^3
and beta_cert = -Psi0/theta.   Enclosures:
  psi(t) = E exp(-tW) = exp(-t w_inf - sum_b lam_b (1 - exp(-t w_b)))  is convex in t, and so is exp(-t w_b) psi(t).
  On each cell [a,b]: chord >= function >= tangent at the midpoint.
  Phi = -gamma - ln t0 + [0, m1 t0] - int_{t0}^{t1} psi dt/t - [0, psi(t1)/(t1 w_inf)]
  rho_b = w_b int_0^inf exp(-t w_b) psi dt,  tail beyond t1 in [0, psi(t1)];  eps = w_inf int_0^inf psi dt, same tail."""
import sys, json, time
from fractions import Fraction
from mpmath import iv
iv.prec = 80
d = json.load(open(sys.argv[1])); delta = float(sys.argv[2]) if len(sys.argv) > 2 else 0.01
def R(x): f = Fraction(str(x)); return iv.mpf(f.numerator) / iv.mpf(f.denominator)
B = d["B"]; eta = R(d["eta"]); K = eta * B; th = R(d["theta"])
u = [R(v) for v in d["u"]]; uinf = R(d["uinf"])
w = [iv.exp(v) for v in u]; winf = iv.exp(uinf)
def G(x): return -iv.exp(-th * x) * (x * x / (2 * th) + x / th ** 2 + 1 / th ** 3)
lam = [((eta * b) ** 3 - (eta * (b - 1)) ** 3) / 6 for b in range(1, B + 1)]
kap = [G(eta * b) - G(eta * (b - 1)) for b in range(1, B + 1)]
J = -G(K)
m1 = winf + sum(l * x for l, x in zip(lam, w))
def ev(t):
    """t: iv number. returns psi, psi', list f_b = exp(-t w_b) psi."""
    e = [iv.exp(-t * x) for x in w]
    S = sum(l * (1 - x) for l, x in zip(lam, e))
    psi = iv.exp(-t * winf - S)
    dpsi = -psi * (winf + sum(l * x * y for l, x, y in zip(lam, w, e)))
    return psi, dpsi, [x * psi for x in e]
t0 = 2.0 ** -27; tic = time.time()
# cell [0, t0] for rho_b and eps
a = iv.mpf(0); b = iv.mpf(t0); pa, _, fa = iv.mpf(1), None, [iv.mpf(1)] * B
pb, _, fb = ev(b); m = (a + b) / 2; pm, dpm, fm = ev(m)
rho_lo = [(b - a) * x for x in fm]; rho_hi = [(b - a) * (x + y) / 2 for x, y in zip(fa, fb)]
eps_lo = (b - a) * pm; eps_hi = (b - a) * (pa + pb) / 2
I_lo = iv.mpf(0); I_hi = iv.mpf(0); n = 0; tf = t0
while True:
    a, pa, fa = b, pb, fb
    tf = tf * (1 + delta); b = iv.mpf(tf); pb, _, fb = ev(b); m = (a + b) / 2; pm, dpm, fm = ev(m)
    h = b - a; lg = iv.log(b / a)
    I_lo += (pm - dpm * m) * lg + dpm * h
    bc = (pb - pa) / h; I_hi += pa * lg + bc * (h - a * lg)
    for i in range(B):
        rho_lo[i] += h * fm[i]; rho_hi[i] += h * (fa[i] + fb[i]) / 2
    eps_lo += h * pm; eps_hi += h * (pa + pb) / 2
    n += 1
    if (pb / (b * winf)).b < 1e-12 and pb.b < 1e-12: break
t1 = b; tail = pb.b
I = iv.mpf([I_lo.a, I_hi.b])
Phi = -iv.euler - iv.log(iv.mpf(t0)) + iv.mpf([0, (m1 * t0).b]) - I - iv.mpf([0, (pb / (t1 * winf)).b])
rho = [iv.mpf([(w[i] * rho_lo[i]).a, (w[i] * rho_hi[i]).b + tail]) for i in range(B)]
eps = iv.mpf([(winf * eps_lo).a, (winf * eps_hi).b + tail])
norm = sum(l * r for l, r in zip(lam, rho)) + eps
def up(x): return x.b
tot = 3 * Phi; worst_rel = 0
for i in range(B):
    C = iv.exp(-3 * u[i]) * kap[i] / lam[i]
    assert rho[i].a >= 0 and rho[i].b <= 1
    vals = [iv.log((1 - r) ** 3 + C * r ** 3) for r in (iv.mpf(rho[i].a), iv.mpf(rho[i].b))]
    tot += lam[i] * iv.mpf(max(up(vals[0]), up(vals[1])))
    worst_rel = max(worst_rel, float(((rho[i].b - rho[i].a) / rho[i].a).b))
Cinf = iv.exp(-3 * uinf) * J
chi = [-3 * e + Cinf * e ** 3 for e in (iv.mpf(eps.a), iv.mpf(eps.b))]
tot += iv.mpf(max(up(chi[0]), up(chi[1])))
Psi0_up = tot.b; beta_cert = (-iv.mpf(Psi0_up) / th).a
print("file=%s delta=%g cells=%d t1=%.4g time=%.0fs" % (sys.argv[1], delta, n, float(t1.b), time.time() - tic))
print("K=%s eta=%s B=%d theta=%s" % (d["K"], d["eta"], B, d["theta"]))
print("Phi in [%s, %s]  width %.2e" % (iv.nstr(iv.mpf(Phi.a), 12), iv.nstr(iv.mpf(Phi.b), 12), float((Phi.b - Phi.a).b)))
print("eps in [%s, %s]; rho_1 in [%s, %s]; worst relative width of rho_b: %.2e" % (iv.nstr(iv.mpf(eps.a), 10), iv.nstr(iv.mpf(eps.b), 10), iv.nstr(iv.mpf(rho[0].a), 10), iv.nstr(iv.mpf(rho[0].b), 10), worst_rel))
print("check sum_b lam_b rho_b + eps in [%s, %s] (must contain 1)" % (iv.nstr(iv.mpf(norm.a), 12), iv.nstr(iv.mpf(norm.b), 12)))
print("ln J in [%s, %s]; m1 <= %s" % (iv.nstr(iv.mpf(iv.log(J).a), 10), iv.nstr(iv.mpf(iv.log(J).b), 10), iv.nstr(iv.mpf(m1.b), 8)))
print("Psi0 <= %s" % iv.nstr(iv.mpf(Psi0_up), 12))
print("beta_cert = -Psi0/theta >= %s" % iv.nstr(iv.mpf(beta_cert), 12))
for bt in ("2.30", "2.3", "2.302"):
    v = (th * R(bt) + iv.mpf(Psi0_up)).b
    print("  beta=%s: theta*beta + Psi0 <= %s  -> %s" % (bt, iv.nstr(iv.mpf(v), 8), "NEGATIVE (certified)" if v < 0 else "not certified"))
