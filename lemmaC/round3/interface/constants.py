"""Error budget for Lemma Q / (I-a*) of THEOREM_I_PA.md, in interval arithmetic (mpmath.iv).
Inputs (exact rationals): K, B, theta (tilt), the certified upper bound PSI0_UP on Psi_0 from
../firstmoment/out_certify_K5_e10.txt, and the choices (beta, c1, d1, d2, tau).
Checks, with outward rounding:
  A      = (8*c1 + 4 + 2*theta*K)/3                       (box size)
  c1    >= 1 + ln(lam_max + d1 + 2*tau)                    (hypothesis of Lemma E)
  Err    = 3*(d2 + 2*tau*(A + ln 2)) + 4*c1*(d1 + 2*tau) + (6B+2)*exp(-1-c1)
  theta*beta + PSI0_UP + Err < 0.
"""
from fractions import Fraction
from mpmath import iv
iv.prec = 80
def R(s): f = Fraction(s); return iv.mpf(f.numerator) / iv.mpf(f.denominator)
K = R("5"); B = 50; eta = R("1/10"); th = R("1.6712"); PSI0_UP = R("-3.84751328")
lam_max = eta ** 3 * (3 * B * B - 3 * B + 1) / 6
J = iv.exp(-th * K) * (K * K / (2 * th) + K / th ** 2 + 1 / th ** 3)
assert (-iv.log(J)).b < (th * K).a, "ell_max = theta*K must dominate |ln J|"
print("lam_max = %s, theta*K = %s, ln J <= %s" % (iv.nstr(lam_max, 8), iv.nstr(th * K, 8), iv.nstr(iv.mpf(iv.log(J).b), 8)))
for beta, c1, d1, d2, tau in (("2.30", "14", "1/100000", "1/5000", "1/200000"),
                              ("2.302", "17", "1/1000000", "1/50000", "1/4000000")):
    beta, c1, d1, d2, tau = map(R, (beta, c1, d1, d2, tau))
    A = (8 * c1 + 4 + 2 * th * K) / 3
    d1p = d1 + 2 * tau; d2p = d2 + 2 * tau * (A + iv.log(2))
    hyp = 1 + iv.log(lam_max + d1p)
    assert hyp.b < c1.a
    tailterm = (6 * B + 2) * iv.exp(-1 - c1)
    Err = 3 * d2p + 4 * c1 * d1p + tailterm
    margin = th * beta + PSI0_UP
    tot = margin + Err
    print("beta=%s c'=%s A<=%s delta1=%s delta2=%s tau=%s" % tuple(iv.nstr(iv.mpf(x.b), 8) for x in (beta, c1, A, d1, d2, tau)))
    print("   delta1'<=%s delta2'<=%s  (6B+2)e^(-1-c')<=%s  Err<=%s" % tuple(iv.nstr(iv.mpf(x.b), 6) for x in (d1p, d2p, tailterm, Err)))
    print("   theta*beta+Psi0 <= %s ;  theta*beta+Psi0+Err <= %s  -> %s" % (iv.nstr(iv.mpf(margin.b), 8), iv.nstr(iv.mpf(tot.b), 8), "NEGATIVE" if tot.b < 0 else "FAIL"))
