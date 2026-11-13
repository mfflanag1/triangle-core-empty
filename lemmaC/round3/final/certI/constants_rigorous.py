"""Hardened copy of ../../interface/constants.py (error budget for Lemma Q / (I-a*) of THEOREM_I_PA.md).
The mathematics is unchanged.  Differences from the original:
 1. exp and log are certify_rigorous.rexp / rlog (interval + - * / only); no iv.exp, iv.log or interval **.
 2. The constant PSI0_UP = -3.84751328 is no longer copied from printed output.  It is asserted, in exact
    rational arithmetic, to be >= the exact upper endpoint of the Psi_0 enclosure, which is taken
      (default)  from psi0_export_w_K5_e10.json written by certify_rigorous.py (weights file checked by SHA-256,
                 parameters K, B, eta, theta checked), or
      (--fresh)  from a fresh run of certify_rigorous.certify in this process (about 2 minutes).
 3. Every printed bound is an exact decimal rounded in the stated direction (dec_le / dec_ge, asserted).
Checks (unchanged), with outward rounding:
  A      = (8*c1 + 4 + 2*theta*K)/3                       (box size)
  c1    >= 1 + ln^+(lam_max + d1 + 2*tau)                  (hypothesis of Lemma E)
  Err    = 3*(d2 + 2*tau*(A + ln 2)) + 4*c1*(d1 + 2*tau) + (6B+2)*exp(-1-c1)
  theta*beta + PSI0_UP + Err < 0.
"""
import sys, os, json, hashlib
from fractions import Fraction
from mpmath import iv
import certify_rigorous as C
from certify_rigorous import R, rexp, rlog, dec_le, dec_ge

HERE = os.path.dirname(os.path.abspath(__file__))
WFILE = os.path.join(HERE, "..", "..", "firstmoment", "w_K5_e10.json"); DELTA = 0.005
K_Q, B, ETA_Q, TH_Q = Fraction(5), 50, Fraction(1, 10), Fraction("1.6712")
PSI0_UP_Q = Fraction("-3.84751328")          # the constant quoted in THEOREM_I.md and THEOREM_I_PA.md

# ---- hand-off of the certified bound on Psi_0 ----
if "--fresh" in sys.argv:
    out = C.certify(WFILE, DELTA, log=lambda s: None); psi0_exact = out["Psi0_up"]
    src = "fresh run of certify_rigorous.certify(w_K5_e10.json, delta=%g)" % DELTA
else:
    out = json.load(open(os.path.join(HERE, "psi0_export_w_K5_e10.json"))); psi0_exact = Fraction(out["Psi0_up"])
    assert out["weights_sha256"] == hashlib.sha256(open(WFILE, "rb").read()).hexdigest(), "export is for another weights file"
    src = "psi0_export_w_K5_e10.json (weights SHA-256 verified)"
assert Fraction(out["K"]) == K_Q and out["B"] == B and Fraction(out["eta"]) == ETA_Q and Fraction(out["theta"]) == TH_Q
assert psi0_exact <= PSI0_UP_Q, "PSI0_UP is not an upper bound for the certified enclosure"
print("Psi0 hand-off: source = %s" % src)
print("   exact upper endpoint = %d/%d  <=  %s   (so Psi0 <= %s; asserted in exact rational arithmetic; slack >= %s)"
      % (psi0_exact.numerator, psi0_exact.denominator, dec_ge(psi0_exact, 11), dec_ge(PSI0_UP_Q, 8), dec_le(PSI0_UP_Q - psi0_exact, 12)))

# ---- error budget ----
K = R(K_Q); eta = R(ETA_Q); th = R(TH_Q); PSI0_UP = R(PSI0_UP_Q)
lam_max = eta * eta * eta * (3 * B * B - 3 * B + 1) / 6
J = rexp(-th * K) * (K * K / (2 * th) + K / (th * th) + 1 / (th * th * th))
lnJ = rlog(J); ln2 = rlog(iv.mpf(2))
assert ((-lnJ).b < (th * K).a) is True and (lnJ.b < (th * K).a) is True, "ell_max = theta*K must dominate |ln J|"
print("lam_max in [%s, %s], theta*K in [%s, %s], ln J in [%s, %s], ln 2 in [%s, %s]"
      % (dec_le(lam_max, 9), dec_ge(lam_max, 9), dec_le(th * K, 6), dec_ge(th * K, 6), dec_le(lnJ, 8), dec_ge(lnJ, 8), dec_le(ln2, 12), dec_ge(ln2, 12)))
ok = True
for beta_s, c1_s, d1_s, d2_s, tau_s in (("2.30", "14", "1/100000", "1/5000", "1/200000"),
                                        ("2.302", "17", "1/1000000", "1/50000", "1/4000000")):
    beta, c1, d1, d2, tau = map(R, (beta_s, c1_s, d1_s, d2_s, tau_s))
    A = (8 * c1 + 4 + 2 * th * K) / 3
    d1p = d1 + 2 * tau; d2p = d2 + 2 * tau * (A + ln2)
    hyp = 1 + rlog(lam_max + d1p)
    assert (hyp.b < c1.a) is True and (c1.a >= 1) is True          # c1 >= 1 + ln^+(lam_max + d1 + 2 tau)
    tailterm = (6 * B + 2) * rexp(-1 - c1)
    Err = 3 * d2p + 4 * c1 * d1p + tailterm
    margin = th * beta + PSI0_UP
    tot = margin + Err
    neg = (tot.b < 0) is True; ok = ok and neg
    print("beta=%s c'=%s delta1=%s delta2=%s tau=%s   A <= %s" % (beta_s, c1_s, d1_s, d2_s, tau_s, dec_ge(A, 6)))
    print("   delta1' <= %s  delta2' <= %s  (6B+2)e^(-1-c') <= %s  Err <= %s" % (dec_ge(d1p, 9), dec_ge(d2p, 12), dec_ge(tailterm, 12), dec_ge(Err, 10)))
    print("   theta*beta+Psi0 <= %s (m >= %s);  theta*beta+Psi0+Err <= %s (m - Err >= %s)  -> %s"
          % (dec_ge(margin, 10), dec_le(-margin, 10), dec_ge(tot, 10), dec_le(-tot, 10), "NEGATIVE (certified)" if neg else "FAIL"))
assert ok
