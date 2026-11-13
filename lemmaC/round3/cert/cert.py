"""Certified upper bound on c* via one-sided (dispersive) rounding.  See CERTIFICATE.md for the proof.

Usage: cert.py HEXP N1 N2 STAGES [outprefix]
  lattice spacing h = 2^-HEXP; mesh of levels: N1 log-spaced targets on [TH_MIN, TH_MAX] plus N2 on [0.01, 8].

Rigorous parts
  * mpmath.iv (prec 120; only its +, -, *, / and comparisons; exp is rexp below): everything about a rounded law (breakpoints x_i, survival values, cell
    integrals, mean), and the final integral of exp(-g).
  * float64, round-to-nearest: sums of products of NONNEGATIVE numbers only (convolutions and the
    sums  sum_m W_m (t+y_m)_+^k).  Arguments t+y_m are exact (dyadic).  Every such sum has at most
    2^17 terms and at most 4 rounded operations per term, so its relative error is below
    2^18 * 2^-53 < 3e-11; results are deflated/inflated by EPS = 1e-9 and an absolute TINY = 1e-290
    (underflow).  iv -> float conversions go through nextafter and are checked exactly.
Plain floating point (no effect on validity): the choice of mesh points t_i.
"""
import sys, time, math
import numpy as np
from mpmath import iv

HEXP = int(sys.argv[1]); N1 = int(sys.argv[2]); N2 = int(sys.argv[3]); STAGES = int(sys.argv[4])
PREFIX = sys.argv[5] if len(sys.argv) > 5 else None
TH_MIN, TH_MAX = 1e-6, 32.0
EPS, TINY = 1e-9, 1e-290
TGRID = 2.0 ** 20
iv.prec = 120
hf = 2.0 ** -HEXP


_C = [iv.mpf(1) / math.factorial(j) for j in range(12)]            # enclosures of 1/j!
_REM = iv.mpf([-1, 1]) * (iv.mpf(2.0 ** -120) / math.factorial(12) * 2)   # |r|^12/12! * e^|r| <= this for |r| <= 2^-10


def rexp(x):
    """Rigorous enclosure of exp(x) for an iv number x, built from interval +, *, / only:
    exp(x) = exp(r)^(2^k), r = x/2^k, |r| <= 2^-10 (asserted), exp(r) in sum_{j<12} r^j/j! + [-rho, rho]
    with rho >= |r|^12/12! * e^|r| (Lagrange remainder).
    mpmath's own iv.exp is NOT used: in mpmath 1.3.0 at prec 120 its upper endpoint was observed to be
    below the true value (by a relative 1e-40, in 5 of 20000 random arguments)."""
    xm = max(abs(float(x.a)), abs(float(x.b)))
    k = 0 if xm == 0.0 else max(0, math.frexp(xm)[1] + 11)
    r = x / iv.mpf(2.0 ** k)
    assert (abs(r) <= iv.mpf(2.0 ** -10)) is True
    P = _C[11]
    for j in range(10, -1, -1):
        P = _C[j] + r * P
    E = P + _REM
    for _ in range(k):
        E = E * E
    return E


def f_hi(x):
    """float >= upper endpoint of the iv number x (checked exactly)."""
    f = float(np.nextafter(float(x.b), np.inf))
    assert (iv.mpf(f) >= x.b) is True
    return f


def f_lo(x):
    """float in [0, lower endpoint of x] (checked exactly); 0 if the lower endpoint is <= 0."""
    f = float(np.nextafter(float(x.a), -np.inf))
    if f <= 0.0:
        return 0.0
    assert (iv.mpf(f) <= x.a) is True
    return f


def lo(a):   # rigorous lower bound for a float64 nonnegative sum-of-products
    return np.maximum(a * (1 - EPS) - TINY, 0.0)


def hi(a):
    return a * (1 + EPS) + TINY


def conv_nonneg(a, b):
    """out[m] = sum_k a[k] b[m-k], plain sequential accumulation (elementwise IEEE mult/add)."""
    n, m = len(a), len(b)
    out = np.zeros(n + m - 1)
    for k in range(n):
        if a[k] != 0.0:
            out[k:k + m] += a[k] * b
    return out


def law_to_lattice(theta, shat):
    """Exact law: G(x) piecewise linear, G(x_i) = theta_i, slope shat_i on (x_{i-1}, x_i], x_N = 0,
    G = 0 at x_min = x_0 - theta_0/shat_0, atom exp(-theta_N) at 0.
    Returns Kmin, p_lo[k], D_lo[k], nu_hi[k] for lattice index k = Kmin..0 (array index k-Kmin), and
    an upper bound on E X.   S = survival, J_k = (1/h) int_{kh}^{(k+1)h} S:
       p_k = S_k - S_{k+1} (mass of cell k),  D_k = h (J_k - S_{k+1}) = int_cell (x - kh) dmu,
       nu_j = J_{j-1} - J_j  (J_{Kmin-1} = 1, J_0 = 0)."""
    N = len(theta) - 1
    h = iv.mpf(hf)
    th = [iv.mpf(float(v)) for v in theta]; sh = [iv.mpf(float(v)) for v in shat]
    x = [None] * (N + 1); x[N] = iv.mpf(0)
    for i in range(N, 0, -1):
        x[i - 1] = x[i] - (th[i] - th[i - 1]) / sh[i]
    xmin = x[0] - th[0] / sh[0]
    bps = [xmin] + x                       # breakpoints, increasing
    lev = [iv.mpf(0)] + th                 # G at breakpoints
    Sbp = [iv.mpf(1)] + [rexp(-v) for v in th]
    slope = [None] + sh                    # slope[j] on (bps[j-1], bps[j])
    Kmin = int(math.floor(float(xmin.a) / hf))
    assert (iv.mpf(Kmin) * h <= xmin.a) is True and (iv.mpf(Kmin + 1) * h > xmin.b) is True
    ncell = -Kmin
    S = [None] * (ncell + 1); J = [None] * ncell
    pos = iv.mpf(Kmin) * h; Spos = iv.mpf(1); j = 0; tot = iv.mpf(0)
    for c in range(ncell):
        k = Kmin + c
        right = iv.mpf(k + 1) * h
        S[c] = Spos; a, Sa = pos, Spos; integ = iv.mpf(0)
        while True:
            cmp = bps[j] < right
            if cmp is None:
                raise RuntimeError("undecidable breakpoint/lattice comparison")
            if not cmp:
                break
            integ += (bps[j] - a) if j == 0 else (Sa - Sbp[j]) / slope[j]
            a, Sa = bps[j], Sbp[j]; j += 1
        # now a <= right <= bps[j], right lies on piece j (j >= 1 because bps[0] = xmin < (Kmin+1) h)
        Sr = rexp(-(lev[j] - slope[j] * (bps[j] - right)))
        integ += (Sa - Sr) / slope[j]
        J[c] = integ / h; tot += integ
        pos, Spos = right, Sr
    S[ncell] = Spos                         # = exp(-theta_N) up to enclosure (left limit at 0)
    assert j == N + 1
    A = Sbp[N + 1]
    # E X = -int_{Kmin h}^0 (1 - S) dx = Kmin*h + tot   (S = 1 on [Kmin h, xmin] is included in tot)
    p_lo = np.zeros(ncell + 1); D_lo = np.zeros(ncell + 1); nu_hi = np.zeros(ncell + 1)
    for c in range(ncell):
        p_lo[c] = f_lo(S[c] - S[c + 1])
        D_lo[c] = f_lo(h * (J[c] - S[c + 1]))
        Jm = iv.mpf(1) if c == 0 else J[c - 1]
        nu_hi[c] = f_hi(Jm - J[c])
    p_lo[ncell] = f_lo(A)
    nu_hi[ncell] = f_hi(J[ncell - 1])
    EX_hi = f_hi(iv.mpf(Kmin) * h + tot)
    return Kmin, p_lo, D_lo, nu_hi, EX_hi, f_hi(-xmin)


def evalsums(t, y, P2lo, CDlo, NU2hi, chunk=64):
    """g_lo, gp_lo, gp_hi at the points t (t + y exact).  g = E(t+Y)_+^3/6, g' = E(t+Y)_+^2/2.
       g  >= sum P2 u^3/6 + sum CD u^2,   g' >= sum P2 u^2/2 + 2 sum CD u,   g' <= sum NU2 u^2/2,
       u = (t + y_m)_+,  P2 = p*p,  CD = D*p,  NU2 = nu*nu."""
    assert np.all(t * TGRID == np.round(t * TGRID)) and np.max(np.abs(t)) < 32 and np.max(np.abs(y)) < 32
    glo = np.empty(len(t)); gplo = np.empty(len(t)); gphi = np.empty(len(t))
    for s in range(0, len(t), chunk):
        tt = t[s:s + chunk]
        m0 = np.searchsorted(y, -tt.max(), side='right')        # y[m] <= -t for all m < m0: u = 0
        U = np.maximum(tt[:, None] + y[None, m0:], 0.0)
        U2 = U * U; U3 = U2 * U
        a, b, c = P2lo[m0:], CDlo[m0:], NU2hi[m0:]
        glo[s:s + chunk] = lo(lo(U3 @ a) / 6 + lo(U2 @ b))
        gplo[s:s + chunk] = lo(lo(U2 @ a) / 2 + 2 * lo(U @ b))
        gphi[s:s + chunk] = hi(hi(U2 @ c) / 2)
    return glo, gplo, gphi


def choose_mesh(y, P2lo, CDlo, NU2hi):
    """Plain floating point: mesh points t_i (multiples of 2^-20) with g(t_i) roughly at the target levels."""
    tc = np.arange(1, 32 * 128) / 128.0
    g = evalsums(tc, y, P2lo, CDlo, NU2hi)[0]
    ok = g > 1e-200
    targets = np.concatenate([np.linspace(math.log(TH_MIN), math.log(TH_MAX), N1),
                              np.linspace(math.log(0.01), math.log(8.0), N2)])
    t = np.interp(np.sort(targets), np.log(g[ok]), tc[ok])
    t = np.unique(np.round(t * TGRID) / TGRID)
    return t[t > 0]


def mean_T_hi(t, glo, gplo):
    """Upper bound on int_0^inf exp(-g):  t_0 + sum_i exp(-glo_i) (1 - exp(-gplo_i d_i))/gplo_i + exp(-glo_N)/gplo_N."""
    tot = iv.mpf(float(t[0]))
    n = len(t)
    for i in range(n):
        a = iv.mpf(float(gplo[i])); e = rexp(-iv.mpf(float(glo[i])))
        assert gplo[i] > 0
        if i < n - 1:
            d = iv.mpf(float(t[i + 1])) - iv.mpf(float(t[i]))
            tot += e * (1 - rexp(-a * d)) / a
        else:
            tot += e / a
    return f_hi(tot)


def spot_check(p_lo, D_lo, nu_hi, P2lo, CDlo, NU2hi, y, t, glo, gplo, gphi, rng, ncheck=2):
    """Exact (interval) re-evaluation of a few entries of the float64 sums: direct evidence for the
    rounding model.  Not needed for the proof; raises if a float bound is on the wrong side."""
    n = len(p_lo)
    for m in rng.integers(0, 2 * n - 1, size=ncheck):
        s1 = iv.mpf(0); s2 = iv.mpf(0); s3 = iv.mpf(0)
        for k in range(max(0, m - n + 1), min(n - 1, m) + 1):
            a = iv.mpf(float(p_lo[m - k]))
            s1 += iv.mpf(float(p_lo[k])) * a; s2 += iv.mpf(float(D_lo[k])) * a
            s3 += iv.mpf(float(nu_hi[k])) * iv.mpf(float(nu_hi[m - k]))
        assert (iv.mpf(float(P2lo[m])) <= s1.a) is True and (iv.mpf(float(CDlo[m])) <= s2.a) is True
        assert (iv.mpf(float(NU2hi[m])) >= s3.b) is True
    for i in rng.integers(0, len(t), size=ncheck):
        a3 = iv.mpf(0); a2 = iv.mpf(0); b2 = iv.mpf(0); b1 = iv.mpf(0); c2 = iv.mpf(0)
        for m in range(len(y)):
            u = float(t[i]) + float(y[m])
            if u <= 0:
                continue
            u = iv.mpf(u); u2 = u * u
            P = iv.mpf(float(P2lo[m])); C = iv.mpf(float(CDlo[m]))
            a3 += P * u2 * u; a2 += P * u2; b2 += C * u2; b1 += C * u; c2 += iv.mpf(float(NU2hi[m])) * u2
        assert (iv.mpf(float(glo[i])) <= (a3 / 6 + b2).a) is True
        assert (iv.mpf(float(gplo[i])) <= (a2 / 2 + 2 * b1).a) is True
        assert (iv.mpf(float(gphi[i])) >= (c2 / 2).b) is True


def main():
    t00 = time.time()
    print("cert.py HEXP=%d h=2^-%d N1=%d N2=%d levels [%g,%g] EPS=%g prec=%d" % (HEXP, HEXP, N1, N2, TH_MIN, TH_MAX, EPS, iv.prec), flush=True)
    Kmin, p_lo, D_lo, nu_hi, EX_hi, width = 0, np.array([1.0]), np.array([0.0]), np.array([1.0]), 0.0, 0.0   # delta_0
    rng = np.random.default_rng(1)
    for n in range(STAGES + 1):
        tic = time.time()
        assert len(p_lo) <= 2 ** 16
        P2lo = lo(conv_nonneg(p_lo, p_lo)); CDlo = lo(conv_nonneg(D_lo, p_lo)); NU2hi = hi(conv_nonneg(nu_hi, nu_hi))
        y = (2 * Kmin + np.arange(len(P2lo))) * hf
        t = choose_mesh(y, P2lo, CDlo, NU2hi)
        glo, gplo, gphi = evalsums(t, y, P2lo, CDlo, NU2hi)
        keep = np.concatenate([[True], np.diff(glo) > 0]) & (glo > 0)
        t, glo, gplo, gphi = t[keep], glo[keep], gplo[keep], gphi[keep]
        assert np.all(np.diff(glo) > 0) and glo[0] > 0
        spot_check(p_lo, D_lo, nu_hi, P2lo, CDlo, NU2hi, y, t, glo, gplo, gphi, rng)
        ET_hi = mean_T_hi(t, glo, gplo)
        c_hi = f_hi(3 * (iv.mpf(ET_hi) + 2 * iv.mpf(EX_hi)))
        print("stage %d: lattice cells %d (support width <= %.4f), mesh %d points, levels [%.3g, %.4g], "
              "sum p_lo = %.12f, sum nu_hi = %.12f, E Q~ <= %.9f, E T(Q~) <= %.9f, CERTIFIED c* <= 3*Gamma(Q~_%d) <= %.6f   [%.0fs]"
              % (n, len(p_lo) - 1, width, len(t), glo[0], glo[-1], p_lo.sum(), nu_hi.sum(), EX_hi, ET_hi, n, c_hi, time.time() - tic), flush=True)
        if n == STAGES:
            break
        theta = glo; shat = np.maximum.accumulate(gphi)
        assert np.all(shat > 0) and np.all(np.diff(shat) >= 0)
        if PREFIX:
            np.savez_compressed("%s_stage%d.npz" % (PREFIX, n + 1), theta=theta, shat=shat, t=t)
        Kmin, p_lo, D_lo, nu_hi, EX_hi, width = law_to_lattice(theta, shat)
    print("total time %.0fs" % (time.time() - t00))


if __name__ == '__main__':
    main()
