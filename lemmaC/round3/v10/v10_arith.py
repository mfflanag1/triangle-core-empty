"""V10 arithmetic attack on cert.py (the objects under test are cert.rexp, cert.law_to_lattice, cert.lo/hi,
cert.conv_nonneg, cert.evalsums; the references are my own: mp at 400 bits, Fractions)."""
import sys, random, math
sys.argv = ['cert.py', '10', '1500', '2500', '8']
sys.path.insert(0, '/Users/maxf/projects/economics-research/game-theory/coopgpt/average_case_protection_2026-11-09/lemmaC/round3/cert')
import numpy as np
import cert
from cert import iv
from mpmath import mp, mpf
from fractions import Fraction
mp.prec = 400


def fr(raw):
    s, m, e, bc = raw; v = Fraction(int(m)) * (Fraction(2) ** e); return -v if s else v


def ends(z):
    return mp.mpf(z._mpi_[0]), mp.mpf(z._mpi_[1])


# 1. rexp on point and proper intervals, special arguments
random.seed(77); bad = 0; n = 0
specials = [0.0, -0.0, 2.0 ** -10, -2.0 ** -10, 2.0 ** -11, -2.0 ** -9, -1.0, -2.0, -4.0, -32.0, -32.00002939420814, -64.0, -700.0, 1.0, 5e-324, -5e-324,
            -1e-300, -2.0 ** -1022, -1e-17, 30.0]
for i in range(6000):
    if i < len(specials):
        a = b = specials[i]
    else:
        a = -random.random() * [40, 1, 1e-3, 1e-9, 1e-30, 700][i % 6]
        b = a + (0 if i % 2 else abs(a) * random.random() * [1e-3, 1e-20, 0.5][i % 3])
    e = cert.rexp(iv.mpf([a, b])); lo_, hi_ = ends(e)
    ok = (lo_ <= mp.exp(mp.mpf(a))) and (mp.exp(mp.mpf(b)) <= hi_) and lo_ > 0
    bad += not ok; n += 1
print("1. rexp: %d failures in %d point/interval arguments (incl. %d specials)" % (bad, n, len(specials)))

# 2. iv + - * / on proper intervals of mixed sign against Fractions
bad = 0; n = 0
for i in range(20000):
    def rnd():
        a = random.uniform(-1, 1) * 10.0 ** random.randint(-30, 6); w = abs(a) * random.choice([0, 1e-25, 1e-3, 3.0])
        return a, a + w
    (a, b), (c, d) = rnd(), rnd()
    X = iv.mpf([a, b]); Y = iv.mpf([c, d])
    A, B, C, D = Fraction(a), Fraction(b), Fraction(c), Fraction(d)
    cases = [(X + Y, A + C, B + D), (X - Y, A - D, B - C), (X * Y, min(A * C, A * D, B * C, B * D), max(A * C, A * D, B * C, B * D))]
    if not (C <= 0 <= D):
        q = [A / C, A / D, B / C, B / D]; cases.append((X / Y, min(q), max(q)))
    for z, l, u in cases:
        bad += not (fr(z._mpi_[0]) <= l and u <= fr(z._mpi_[1])); n += 1
print("2. iv + - * / : %d failures in %d interval operations; iv.mpf(float) exact: %s"
      % (bad, n, all(fr(iv.mpf(v)._mpi_[0]) == Fraction(v) == fr(iv.mpf(v)._mpi_[1]) for v in [0.1, 1e-300, 5e-324, 1 / 3.0, 32.00002939420814])))

# 3. law_to_lattice outputs against my own 400-bit evaluation of the law (stage 5 of run A)
d = np.load('/Users/maxf/projects/economics-research/game-theory/coopgpt/average_case_protection_2026-11-09/lemmaC/round3/cert/stages/A_stage5.npz')
theta, shat = d['theta'], d['shat']
Kmin, p_lo, D_lo, nu_hi, EX_hi, width = cert.law_to_lattice(theta, shat)
N = len(theta) - 1
th = [mpf(float(v)) for v in theta]; sh = [mpf(float(v)) for v in shat]
xs = [None] * (N + 1); xs[N] = mpf(0)
for i in range(N, 0, -1):
    xs[i - 1] = xs[i] - (th[i] - th[i - 1]) / sh[i]
xmin = xs[0] - th[0] / sh[0]
h = mpf(2) ** -10
# cumulative integral of S from xmin, and S, at arbitrary x in [Kmin h, 0]: own closed forms
cumI = [mpf(0)]                       # int_{xmin}^{x_i} S
prevx, prevS = xmin, mpf(1)
Sx = [mp.exp(-v) for v in th]
for i in range(N + 1):
    cumI.append(cumI[-1] + (prevS - Sx[i]) / sh[i]); prevS = Sx[i]
xsf = [float(v) for v in xs]
import bisect


def S_and_I(x):
    """S(x) and int_{Kmin h}^x S for x <= 0 (left limit at 0)."""
    base = xmin - Kmin * h
    if x <= xmin:
        return mpf(1), x - Kmin * h
    i = bisect.bisect_left(xsf, float(x))          # first i with x_i >= x (float guess, corrected below)
    while i > 0 and xs[i - 1] >= x: i -= 1
    while xs[i] < x: i += 1
    Sv = mp.exp(-(th[i] - sh[i] * (xs[i] - x)))
    left = xmin if i == 0 else xs[i - 1]; Sl = mpf(1) if i == 0 else Sx[i - 1]
    return Sv, base + cumI[i] + (Sl - Sv) / sh[i]


ncell = -Kmin
SI = [S_and_I((Kmin + c) * h) for c in range(ncell + 1)]
worst = [mpf(0), mpf(0), mpf(0)]; relgap = [0.0, 0.0, 0.0]
for c in range(ncell):
    p = SI[c][0] - SI[c + 1][0]; J = (SI[c + 1][1] - SI[c][1]) / h; Dk = h * (J - SI[c + 1][0])
    Jm = mpf(1) if c == 0 else (SI[c][1] - SI[c - 1][1]) / h
    nu = Jm - J
    worst[0] = max(worst[0], mpf(float(p_lo[c])) - p); worst[1] = max(worst[1], mpf(float(D_lo[c])) - Dk); worst[2] = max(worst[2], nu - mpf(float(nu_hi[c])))
    relgap[0] = max(relgap[0], float((p - mpf(float(p_lo[c]))) / p)); relgap[1] = max(relgap[1], float((Dk - mpf(float(D_lo[c]))) / Dk))
    relgap[2] = max(relgap[2], float((mpf(float(nu_hi[c])) - nu) / nu))
atom = Sx[N]; Jlast = (SI[ncell][1] - SI[ncell - 1][1]) / h
EXref = Kmin * h + SI[ncell][1]
print("3. stage-5 lattice (%d cells): max(p_lo - p) = %s, max(D_lo - D) = %s, max(nu - nu_hi) = %s (all must be <= 0); largest relative slack %.1e %.1e %.1e"
      % (ncell, mp.nstr(worst[0], 3), mp.nstr(worst[1], 3), mp.nstr(worst[2], 3), relgap[0], relgap[1], relgap[2]))
print("   atom: p_lo - atom = %s (<= 0), nu_hi[0] - J_{-1} = %s (>= 0); E X: certified upper %.15f, mine %s, diff %s (>= 0)"
      % (mp.nstr(mpf(float(p_lo[ncell])) - atom, 3), mp.nstr(mpf(float(nu_hi[ncell])) - Jlast, 3), EX_hi, mp.nstr(EXref, 18), mp.nstr(mpf(EX_hi) - EXref, 3)))

# 4. float64 sums against exact rational arithmetic (same law): convolution entries and the three sums
conv1 = cert.conv_nonneg(p_lo, p_lo); conv2 = cert.conv_nonneg(D_lo, p_lo); conv3 = cert.conv_nonneg(nu_hi, nu_hi)
P2lo = cert.lo(conv1); CDlo = cert.lo(conv2); NU2hi = cert.hi(conv3)
F = lambda arr: [Fraction(float(v)) for v in arr]
Fp, FD, Fn = F(p_lo), F(D_lo), F(nu_hi)
nn = len(p_lo); relerr = 0.0; bad = 0
ms = [0, 1, 2, nn // 3, nn // 2, nn - 1, nn, nn + 7, 3 * nn // 2, 2 * nn - 3, 2 * nn - 2] + [random.randrange(2 * nn - 1) for _ in range(20)]
for m in ms:
    ks = range(max(0, m - nn + 1), min(nn - 1, m) + 1)
    e1 = sum(Fp[k] * Fp[m - k] for k in ks); e2 = sum(FD[k] * Fp[m - k] for k in ks); e3 = sum(Fn[k] * Fn[m - k] for k in ks)
    bad += not (Fraction(float(P2lo[m])) <= e1 and Fraction(float(CDlo[m])) <= e2 and Fraction(float(NU2hi[m])) >= e3)
    for c, e in ((conv1[m], e1), (conv2[m], e2), (conv3[m], e3)):
        if e > 0: relerr = max(relerr, abs(float((Fraction(float(c)) - e) / e)))
print("4a. convolutions: %d wrong-side results in %d entries x 3; largest observed relative float error %.2e (bound used: 3e-11, deflation 1e-9)" % (bad, len(ms), relerr))
y = (2 * Kmin + np.arange(len(P2lo))) * cert.hf
d6 = np.load('/Users/maxf/projects/economics-research/game-theory/coopgpt/average_case_protection_2026-11-09/lemmaC/round3/cert/stages/A_stage6.npz')
t = d6['t']
assert all(Fraction(float(v)) * 2 ** 20 % 1 == 0 for v in t) and all(Fraction(float(v)) * 2 ** 10 % 1 == 0 for v in y)
glo, gplo, gphi = cert.evalsums(t, y, P2lo, CDlo, NU2hi)
print("    stage-6 parameters reproduced from the stage-5 law: theta %s, shat %s" % (np.array_equal(glo, d6['theta']), np.array_equal(np.maximum.accumulate(gphi), d6['shat'])))
FP, FC, FN, Fy = F(P2lo), F(CDlo), F(NU2hi), F(y)
bad = 0; marg = [1.0, 1.0, 1.0]
for i in [0, 1, 2, len(t) // 4, len(t) // 2, 3 * len(t) // 4, len(t) - 2, len(t) - 1] + [random.randrange(len(t)) for _ in range(12)]:
    ti = Fraction(float(t[i])); a3 = a2 = b2 = b1 = c2 = Fraction(0)
    for m in range(len(Fy)):
        u = ti + Fy[m]
        if u <= 0: continue
        u2 = u * u
        a3 += FP[m] * u2 * u; a2 += FP[m] * u2; b2 += FC[m] * u2; b1 += FC[m] * u; c2 += FN[m] * u2
    eg, egp, egh = a3 / 6 + b2, a2 / 2 + 2 * b1, c2 / 2
    bad += not (Fraction(float(glo[i])) <= eg and Fraction(float(gplo[i])) <= egp and Fraction(float(gphi[i])) >= egh)
    marg = [min(marg[0], float(1 - Fraction(float(glo[i])) / eg)), min(marg[1], float(1 - Fraction(float(gplo[i])) / egp)), min(marg[2], float(Fraction(float(gphi[i])) / egh - 1))]
print("4b. sums at 20 mesh points (exact rationals): %d wrong-side results; smallest relative safety margins: g_lo %.3e, g'_lo %.3e, g'_hi %.3e (designed: about 2e-9, 2e-9, 2e-9)" % (bad, marg[0], marg[1], marg[2]))
