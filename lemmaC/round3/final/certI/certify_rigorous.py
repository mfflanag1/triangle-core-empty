"""Hardened copy of ../../firstmoment/certify.py (certificate of THEOREM_I.md, section 5).
The mathematics is unchanged.  Three things differ from the original:

 1. No mpmath library function is trusted beyond interval + - * / , comparisons and the outward-rounding
    conversion of ints and floats.  iv.exp, iv.log, iv.euler and the interval power ** are not used:
      rexp   exp(x) = exp(r)^(2^k), r = x/2^k, |r| <= 2^-10 (asserted on the interval endpoints),
             exp(r) in sum_{j<12} r^j/j! + [-rho, rho], rho = 2^-119/12! >= |r|^12/12! * e^|r|  (Lagrange).
      rlog   for a point p > 0 and ANY real y0 (here the double math.log(p), used only as a guess):
             ln p = y0 + ln(1 + d), d = p/exp(y0) - 1, and d - d^2 <= ln(1 + d) <= d for |d| <= 1/2;
             d is enclosed with rexp, |d| <= 2^-20 is asserted.  A general interval is handled endpoint by
             endpoint (ln is increasing).
      GAMMA  Euler's constant: gamma = H_n - ln(n + 1/2) - S_n with 0 <= S_n <= 2/(45 n^2), n = 2^16.
             Proof: c_k = H_k - ln(k + 1/2) -> gamma and c_{k-1} - c_k = ln((k+1/2)/(k-1/2)) - 1/k
             = 2 sum_{j>=1} y^(2j+1)/(2j+1) with y = 1/(2k), which lies in [0, (2y^3/3)/(1-y^2)]
             and (2y^3/3)/(1-y^2) <= (16/15)/(12 k^3) for k >= 2; sum_{k>n} k^-3 <= 1/(2n^2).
 2. Every printed bound is an exact decimal obtained from the exact binary endpoint by rational floor/ceil
    in the stated direction (dec_le / dec_ge, asserted).  The bounds used downstream are also written, as
    exact rationals and as outward-rounded decimals, to psi0_export_<weights>.json, and constants_rigorous.py
    asserts its constant against that exact rational (or against a fresh run of certify()).
 3. The weights file is read with exact decimal parsing (no float round trip).

Formulas enclosed (unchanged):  psi(t) = E exp(-tW) = exp(-t w_inf - sum_b lam_b (1 - exp(-t w_b))) is convex
in t, and so is exp(-t w_b) psi(t).  On each cell [a,b]: chord >= function >= tangent at the midpoint.
  Phi = -gamma - ln t0 + [0, m1 t0] - int_{t0}^{t1} psi dt/t - [0, psi(t1)/(t1 w_inf)]
  rho_b = w_b int_0^inf exp(-t w_b) psi dt, tail beyond t1 in [0, psi(t1)];  eps = w_inf int_0^inf psi dt, same tail.
  Psi0 = 3 Phi + sum_b lam_b ln((1-rho_b)^3 + C_b rho_b^3) - 3 eps + C_inf eps^3.

Usage: certify_rigorous.py WEIGHTS.json [delta]      (delta = mesh ratio - 1, default 0.01; the certificates use 0.005)
"""
import sys, os, json, time, math, hashlib
from fractions import Fraction
from mpmath import iv
from mpmath.libmp import fzero

iv.prec = 80
HERE = os.path.dirname(os.path.abspath(__file__))


def R(x):
    """interval enclosure of an exact rational given as int, Fraction or decimal string."""
    f = Fraction(x)
    return iv.mpf(f.numerator) / iv.mpf(f.denominator)


# ---------- exact endpoints and outward-rounded decimals ----------
def _frac(raw):
    s, m, e, bc = raw
    assert m != 0 or raw == fzero, "endpoint is not a finite number"
    v = Fraction(int(m)) * Fraction(2) ** e
    return -v if s else v


def lo_frac(x): return _frac(x._mpi_[0])     # exact lower endpoint of the iv number x
def hi_frac(x): return _frac(x._mpi_[1])     # exact upper endpoint


def _fmt(n, nd):
    return "%s%d.%0*d" % ("-" if n < 0 else "", abs(n) // 10 ** nd, nd, abs(n) % 10 ** nd)


def dec_ge(x, nd):
    """exact decimal with nd digits after the point that is >= the upper endpoint of x (x: iv number or Fraction)."""
    f = x if isinstance(x, Fraction) else hi_frac(x)
    n = math.ceil(f * 10 ** nd)
    assert Fraction(n, 10 ** nd) >= f
    return _fmt(n, nd)


def dec_le(x, nd):
    """exact decimal with nd digits after the point that is <= the lower endpoint of x."""
    f = x if isinstance(x, Fraction) else lo_frac(x)
    n = math.floor(f * 10 ** nd)
    assert Fraction(n, 10 ** nd) <= f
    return _fmt(n, nd)


# ---------- verified exp, log, Euler's constant (interval + - * / only) ----------
_C = [iv.mpf(1) / math.factorial(j) for j in range(12)]                  # enclosures of 1/j!
_REM = iv.mpf([-1, 1]) * (iv.mpf(2.0 ** -120) / math.factorial(12) * 2)    # |r|^12/12! * e^|r| <= 2^-119/12! for |r| <= 2^-10
_RMAX = iv.mpf(2.0 ** -10)


def rexp(x):
    """rigorous enclosure of {exp(t): t in x} for an iv number x."""
    xm = max(abs(float(x.a)), abs(float(x.b)))
    k = 0 if xm == 0.0 else max(0, math.frexp(xm)[1] + 11)
    r = x / iv.mpf(2.0 ** k)
    assert (r.a >= -_RMAX) is True and (r.b <= _RMAX) is True
    P = _C[11]
    for j in range(10, -1, -1):
        P = _C[j] + r * P
    E = P + _REM
    for _ in range(k):
        E = E * E
    return E


def _rlog_point(p):
    """rigorous enclosure of ln p for a positive point p (an iv number of zero width)."""
    pf = float(p)
    assert pf > 0.0 and math.isfinite(pf)
    y0 = iv.mpf(math.log(pf))            # a double, exact as an iv number; only a guess
    d = p / rexp(y0) - 1                 # p = exp(y0) * (1 + d)
    assert (d.a > -2.0 ** -20) is True and (d.b < 2.0 ** -20) is True
    return iv.mpf([(y0 + (d.a - d.a * d.a)).a, (y0 + d.b).b])


def rlog(x):
    """rigorous enclosure of {ln t: t in x} for an iv number x > 0."""
    assert (x.a > 0) is True
    return iv.mpf([_rlog_point(x.a).a, _rlog_point(x.b).b])


def euler_gamma(n=1 << 16):
    H = iv.mpf(0)
    for k in range(1, n + 1):
        H += iv.mpf(1) / k
    g = H - rlog(iv.mpf(2 * n + 1) / 2) - iv.mpf([0, 1]) * (iv.mpf(2) / (45 * n * n))
    # consistency with the literature value (not used in the proof)
    assert lo_frac(g) < Fraction("0.5772156649015328606") < hi_frac(g)
    return g


def cube(x): return x * x * x            # used for nonnegative arguments only


# ---------- the certificate ----------
def certify(wfile, delta=0.01, log=print):
    """returns a dict; 'Psi0_up' (Fraction) is the exact upper endpoint of the enclosure of Psi_0."""
    raw = open(wfile, "rb").read()
    d = json.loads(raw.decode(), parse_float=str)            # decimals kept as strings, read as exact rationals by R
    tic = time.time()
    GAMMA = euler_gamma()
    B = d["B"]; eta = R(d["eta"]); K = eta * B; th = R(d["theta"])
    u = [R(v) for v in d["u"]]; uinf = R(d["uinf"])
    assert len(u) == B and Fraction(d["eta"]) * B == Fraction(d["K"])
    w = [rexp(v) for v in u]; winf = rexp(uinf)
    def G(x): return -rexp(-th * x) * (x * x / (2 * th) + x / (th * th) + 1 / (th * th * th))
    lam = [(cube(eta * b) - cube(eta * (b - 1))) / 6 for b in range(1, B + 1)]
    kap = [G(eta * b) - G(eta * (b - 1)) for b in range(1, B + 1)]
    J = -G(K)
    m1 = winf + sum(l * x for l, x in zip(lam, w))
    def ev(t):
        """t: iv number. returns psi, psi', list f_b = exp(-t w_b) psi."""
        e = [rexp(-t * x) for x in w]
        S = sum(l * (1 - x) for l, x in zip(lam, e))
        psi = rexp(-t * winf - S)
        dpsi = -psi * (winf + sum(l * x * y for l, x, y in zip(lam, w, e)))
        return psi, dpsi, [x * psi for x in e]
    t0 = 2.0 ** -27
    # cell [0, t0] for rho_b and eps
    a = iv.mpf(0); b = iv.mpf(t0); pa, fa = iv.mpf(1), [iv.mpf(1)] * B
    pb, _, fb = ev(b); m = (a + b) / 2; pm, dpm, fm = ev(m)
    rho_lo = [(b - a) * x for x in fm]; rho_hi = [(b - a) * (x + y) / 2 for x, y in zip(fa, fb)]
    eps_lo = (b - a) * pm; eps_hi = (b - a) * (pa + pb) / 2
    I_lo = iv.mpf(0); I_hi = iv.mpf(0); n = 0; tf = t0
    while True:
        a, pa, fa = b, pb, fb
        tf = tf * (1 + delta); b = iv.mpf(tf); pb, _, fb = ev(b); m = (a + b) / 2; pm, dpm, fm = ev(m)
        assert (b.a > a.b) is True
        h = b - a; lg = rlog(b / a)
        I_lo += (pm - dpm * m) * lg + dpm * h
        bc = (pb - pa) / h; I_hi += pa * lg + bc * (h - a * lg)
        for i in range(B):
            rho_lo[i] += h * fm[i]; rho_hi[i] += h * (fa[i] + fb[i]) / 2
        eps_lo += h * pm; eps_hi += h * (pa + pb) / 2
        n += 1
        if (pb / (b * winf)).b < 1e-12 and pb.b < 1e-12: break
    t1 = b; tail = pb.b
    I = iv.mpf([I_lo.a, I_hi.b])
    Phi = -GAMMA - rlog(iv.mpf(t0)) + iv.mpf([0, (m1 * t0).b]) - I - iv.mpf([0, (pb / (t1 * winf)).b])
    rho = [iv.mpf([(w[i] * rho_lo[i]).a, (w[i] * rho_hi[i] + tail).b]) for i in range(B)]
    eps = iv.mpf([(winf * eps_lo).a, (winf * eps_hi + tail).b])
    norm = sum(l * r for l, r in zip(lam, rho)) + eps
    def up(x): return x.b
    tot = 3 * Phi; worst_rel = 0
    for i in range(B):
        C = rexp(-3 * u[i]) * kap[i] / lam[i]
        assert (rho[i].a >= 0) is True and (rho[i].b <= 1) is True and (C.a > 0) is True
        vals = [rlog(cube(1 - r) + C * cube(r)) for r in (rho[i].a, rho[i].b)]
        tot += lam[i] * max(up(vals[0]), up(vals[1]))
        worst_rel = max(worst_rel, float(((rho[i].b - rho[i].a) / rho[i].a).b))
    Cinf = rexp(-3 * uinf) * J
    assert (eps.a >= 0) is True and (Cinf.a > 0) is True
    chi = [-3 * e + Cinf * cube(e) for e in (eps.a, eps.b)]
    tot += max(up(chi[0]), up(chi[1]))
    Psi0_up = tot.b; beta_cert = (-Psi0_up / th).a
    lnJ = rlog(J)
    name = os.path.basename(wfile)
    log("file=%s delta=%g cells=%d t1=%.4g prec=%d time=%.0fs" % (name, delta, n, float(t1), iv.prec, time.time() - tic))
    log("K=%s eta=%s B=%d theta=%s" % (d["K"], d["eta"], B, d["theta"]))
    log("gamma in [%s, %s]" % (dec_le(GAMMA, 14), dec_ge(GAMMA, 14)))
    log("Phi in [%s, %s]  width <= %s" % (dec_le(Phi, 12), dec_ge(Phi, 12), dec_ge(Phi.b - Phi.a, 9)))
    log("eps in [%s, %s]; rho_1 in [%s, %s]; worst relative width of rho_b: %.2e" % (dec_le(eps, 11), dec_ge(eps, 11), dec_le(rho[0], 10), dec_ge(rho[0], 10), worst_rel))
    assert (norm.a <= 1) is True and (norm.b >= 1) is True
    log("check sum_b lam_b rho_b + eps in [%s, %s] (contains 1: asserted)" % (dec_le(norm, 11), dec_ge(norm, 11)))
    log("ln J in [%s, %s]; m1 <= %s" % (dec_le(lnJ, 9), dec_ge(lnJ, 9), dec_ge(m1, 8)))
    log("Psi0 <= %s   (8 decimals: Psi0 <= %s)" % (dec_ge(Psi0_up, 11), dec_ge(Psi0_up, 8)))
    log("beta_cert = -Psi0/theta >= %s   (8 decimals: >= %s)" % (dec_le(beta_cert, 11), dec_le(beta_cert, 8)))
    margins = {}
    for bt in ("2.30", "2.302", "2.3022", "2.30224"):
        v = (th * R(bt) + Psi0_up).b
        margins[bt] = dec_ge(v, 10)
        log("  beta=%s: theta*beta + Psi0 <= %s  -> %s" % (bt, margins[bt], "NEGATIVE (certified)" if (v < 0) is True else "not certified"))
    out = {"script": "certify_rigorous.py", "weights_file": name, "weights_sha256": hashlib.sha256(raw).hexdigest(),
           "delta": delta, "prec": iv.prec, "cells": n,
           "K": str(d["K"]), "B": B, "eta": str(d["eta"]), "theta": str(d["theta"]),
           "Psi0_up": hi_frac(Psi0_up), "beta_cert": lo_frac(beta_cert), "margins": margins}
    return out


def export(out, path):
    """write the hand-off file; every decimal in it is asserted against the exact endpoint."""
    p, bc = out["Psi0_up"], out["beta_cert"]
    e = dict(out)
    e["Psi0_up"] = "%d/%d" % (p.numerator, p.denominator)             # exact (binary) upper endpoint
    e["beta_cert"] = "%d/%d" % (bc.numerator, bc.denominator)         # exact lower endpoint of -Psi0_up/theta
    e["Psi0_up_dec8"] = dec_ge(p, 8); e["Psi0_up_dec11"] = dec_ge(p, 11)
    e["beta_cert_dec8"] = dec_le(bc, 8); e["beta_cert_dec11"] = dec_le(bc, 11)
    assert Fraction(e["Psi0_up_dec8"]) >= p and Fraction(e["Psi0_up_dec11"]) >= p
    assert Fraction(e["beta_cert_dec8"]) <= bc and Fraction(e["beta_cert_dec11"]) <= bc
    assert Fraction(e["Psi0_up"]) == p and Fraction(e["beta_cert"]) == bc
    json.dump(e, open(path, "w"), indent=1)
    return e


if __name__ == "__main__":
    wfile = sys.argv[1]; delta = float(sys.argv[2]) if len(sys.argv) > 2 else 0.01
    out = certify(wfile, delta, log=lambda s: print(s, flush=True))
    path = os.path.join(HERE, "psi0_export_%s" % os.path.basename(wfile))
    e = export(out, path)
    print("exported to %s: Psi0 <= %s (exact rational %s), -Psi0/theta >= %s" % (os.path.basename(path), e["Psi0_up_dec8"], e["Psi0_up"], e["beta_cert_dec8"]))
