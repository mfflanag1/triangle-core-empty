"""V8 referee: independent re-evaluation of the Lemma Q error budget (exact rationals + mp bounds rounded outward)."""
from fractions import Fraction as F
import mpmath as mp
mp.mp.dps = 60
def up(x):   # rational upper bound of an mp number
    return F(int(mp.ceil(x * mp.mpf(10) ** 40)), 10 ** 40)
def dn(x):
    return F(int(mp.floor(x * mp.mpf(10) ** 40)), 10 ** 40)
K = F(5); B = 50; eta = F(1, 10); th = F(16712, 10000); PSI0_UP = F(-384751328, 10 ** 8)
lam = [eta ** 3 * (3 * b * b - 3 * b + 1) / 6 for b in range(1, B + 1)]
lam_max = max(lam); assert sum(lam) == K ** 3 / 6
# ell_b = ln(kappa_b/lam_b) in [-theta*K, 0] (a conditional Laplace transform); ln J:
thf = mp.mpf(th.numerator) / th.denominator
def G(x): return -mp.exp(-thf * x) * (x * x / (2 * thf) + x / thf ** 2 + 1 / thf ** 3)
ell = []
for b in range(1, B + 1):
    kap = G(mp.mpf(b) / 10) - G(mp.mpf(b - 1) / 10)
    ell.append(mp.log(kap / (mp.mpf(lam[b - 1].numerator) / lam[b - 1].denominator)))
lnJ = mp.log(-G(mp.mpf(5)))
ell_max_true = max(max(abs(x) for x in ell), abs(lnJ))
print("lam_max =", lam_max, "=", float(lam_max), "; max|ell_b| =", mp.nstr(max(abs(x) for x in ell), 10), "; ln J =", mp.nstr(lnJ, 10), "; theta*K =", float(th * K))
assert ell_max_true < thf * 5
ln2_up = up(mp.log(2))
for beta, c1, d1, d2, tau in ((F(230, 100), 14, F(1, 10 ** 5), F(2, 10 ** 4), F(5, 10 ** 6)),
                              (F(2302, 1000), 17, F(1, 10 ** 6), F(2, 10 ** 5), F(25, 10 ** 8))):
    A = (8 * c1 + 4 + 2 * th * K) / 3
    d1p = d1 + 2 * tau
    d2p = d2 + 2 * tau * (A + ln2_up)
    hyp = 1 + mp.log(mp.mpf((lam_max + d1p).numerator) / (lam_max + d1p).denominator)
    assert hyp < c1
    tail = (6 * B + 2) * up(mp.exp(-1 - c1))
    Err = 3 * d2p + 4 * c1 * d1p + tail
    m = -(th * beta + PSI0_UP)
    tot = -m + Err
    # coefficient of |mu-mu'|_1 with the TRUE ell_max (must be <= 0)
    coef = 4 * c1 + 2 + ell_max_true - 3 * mp.mpf(A.numerator) / A.denominator / 2
    print("beta=%s c'=%d A=%.6f  d1'=%.3e d2'=%.6e tail=%.6e" % (float(beta), c1, float(A), float(d1p), float(d2p), float(tail)))
    print("   3*d2'=%.6e 4c'd1'=%.6e  Err<=%.9e  m>=%.9e  theta*beta+Psi0+Err <= %.9e  %s ; coef(|mu-mu'|)=%s" % (
        float(3 * d2p), float(4 * c1 * d1p), float(Err), float(m), float(tot), "NEGATIVE" if tot < 0 else "FAIL", mp.nstr(coef, 6)))
    # error term of Theorem T: need (tau/4) q ln q to beat 0.22 q: threshold
    print("   ln q needed for (tau/4) ln q > 0.22 + (m-Err): ln q > %.3e" % (4 * 0.23 / float(tau)))
# plain first moment constant
for bt in (2.30, 2.302, 2.31):
    print("beta=%.3f  ln(beta^3 e/27) = %.5f" % (bt, mp.log(bt ** 3 * mp.e / 27)))
