"""Validation of the validator (floating point, evidence only): on a deliberately COARSE configuration
(lattice h = 2^-HEXP with HEXP = 3, 4, 5; about 20 pieces) compare cert.py's one-sided bounds with a
reference computed by a different method: tensor Gauss-Legendre quadrature (2 x 16 nodes per piece) of the
exact piecewise-exponential law.  With a coarse lattice the bounds are visibly not tight, so a wrong
direction anywhere (tangent/chord bounds, cell masses, iv->float rounding) would show.
Usage: test_small.py HEXP"""
import sys
import numpy as np
HEXP = int(sys.argv[1])
sys.argv = ['cert.py', str(HEXP), '12', '8', '1']
import cert
from cert import iv

GLx, GLw = np.polynomial.legendre.leggauss(16)


def ref_nodes(theta, shat):
    """quadrature nodes/weights of the exact law (float): pieces + atom at 0."""
    dxs = np.diff(theta) / shat[1:]
    x = np.concatenate([[0.0], -np.cumsum(dxs[::-1])])[::-1]
    xb = np.concatenate([[x[0] - theta[0] / shat[0]], x]); lv = np.concatenate([[0.0], theta]); sl = shat
    X = [np.array([0.0])]; W = [np.array([np.exp(-theta[-1])])]
    for j in range(len(sl)):
        a, b = xb[j], xb[j + 1]
        # split each piece further so that kinks of (t+x+x')_+ matter little
        for (aa, bb) in zip(np.linspace(a, b, 3)[:-1], np.linspace(a, b, 3)[1:]):
            xs = 0.5 * (aa + bb) + 0.5 * (bb - aa) * GLx
            X.append(xs); W.append(0.5 * (bb - aa) * GLw * sl[j] * np.exp(-(lv[j + 1] - sl[j] * (xb[j + 1] - xs))))
    return np.concatenate(X), np.concatenate(W)


def ref_g(t, X, W):
    Y = (X[:, None] + X[None, :]).ravel(); P = (W[:, None] * W[None, :]).ravel()
    g = np.array([(P * np.maximum(tt + Y, 0) ** 3).sum() / 6 for tt in t])
    gp = np.array([(P * np.maximum(tt + Y, 0) ** 2).sum() / 2 for tt in t])
    return g, gp


# stage 1 parameters on a coarse mesh (from delta_0), then the law Q~_1 and its lattice bounds
y0 = np.array([0.0]); one = np.array([1.0]); zero = np.array([0.0])
t1 = cert.choose_mesh(y0, one, zero, one)
glo, gplo, gphi = cert.evalsums(t1, y0, cert.lo(one), zero, cert.hi(one))
theta, shat = glo, np.maximum.accumulate(gphi)
print("coarse law Q~_1: %d pieces, levels [%.2g, %.3g], h = 2^-%d" % (len(theta), theta[0], theta[-1], HEXP))
Kmin, p_lo, D_lo, nu_hi, EX_hi, width = cert.law_to_lattice(theta, shat)
X, W = ref_nodes(theta, shat)
print("reference: total mass %.12f, E X = %.9f ; certified E X <= %.9f ; sum p_lo = %.12f, sum nu_hi = %.12f"
      % (W.sum(), (X * W).sum(), EX_hi, p_lo.sum(), nu_hi.sum()))
assert EX_hi >= (X * W).sum() - 1e-10
# cell masses and first moments against a reference: survival function sampled on a fine uniform grid
# (4096 trapezoid panels per lattice cell), independent of the closed forms used in cert.py
dxs = np.diff(theta) / shat[1:]
xx = np.concatenate([[0.0], -np.cumsum(dxs[::-1])])[::-1]
xb = np.concatenate([[xx[0] - theta[0] / shat[0]], xx]); lv = np.concatenate([[0.0], theta])
ncell = -Kmin; sub = 4096
grid = Kmin * cert.hf + np.arange(ncell * sub + 1) * (cert.hf / sub)
Sg = np.exp(-np.interp(grid, xb, lv, left=0.0)).reshape(-1)       # S(0-) = exp(-theta_N) at the last point
Sc = Sg[::sub]
Jc = np.array([(Sg[c * sub:(c + 1) * sub + 1].sum() - 0.5 * (Sg[c * sub] + Sg[(c + 1) * sub])) / sub for c in range(ncell)])
pm = np.concatenate([Sc[:-1] - Sc[1:], [Sc[-1]]]); Dm = np.concatenate([cert.hf * (Jc - Sc[1:]), [0.0]])
nu_ref = np.concatenate([[1.0], Jc]) - np.concatenate([Jc, [0.0]])
print("max (p_lo - p_ref) = %.2e (want <= 0), max (D_lo - D_ref) = %.2e (want <= 0), min (nu_hi - nu_ref) = %.2e (want >= 0); "
      "max |p_lo - p_ref| = %.1e, max |D_lo - D_ref| = %.1e, max |nu_hi - nu_ref| = %.1e  (reference accuracy about 1e-9)"
      % ((p_lo - pm).max(), (D_lo - Dm).max(), (nu_hi - nu_ref).min(), np.abs(p_lo - pm).max(), np.abs(D_lo - Dm).max(), np.abs(nu_hi - nu_ref).max()))
assert (p_lo - pm).max() < 1e-8 and (D_lo - Dm).max() < 1e-8 and (nu_hi - nu_ref).min() > -1e-8
# g, g' bounds against the reference
P2lo = cert.lo(cert.conv_nonneg(p_lo, p_lo)); CDlo = cert.lo(cert.conv_nonneg(D_lo, p_lo)); NU2hi = cert.hi(cert.conv_nonneg(nu_hi, nu_hi))
y = (2 * Kmin + np.arange(len(P2lo))) * cert.hf
t = cert.choose_mesh(y, P2lo, CDlo, NU2hi)
t = np.unique(np.concatenate([t, np.round(np.linspace(0.02, t[-1], 40) * 2 ** 20) / 2 ** 20]))
glo, gplo, gphi = cert.evalsums(t, y, P2lo, CDlo, NU2hi)
g, gp = ref_g(t, X, W)
ok = g > 1e-12
r1 = glo[ok] / g[ok]; r2 = gplo[ok] / gp[ok]; r3 = gphi[ok] / gp[ok]
print("%d test points t in [%.3f, %.3f]:  g_lo/g in [%.6f, %.6f] (want <= 1);  g'_lo/g' in [%.6f, %.6f] (want <= 1);  "
      "g'_hi/g' in [%.6f, %.6f] (want >= 1)" % (ok.sum(), t[ok][0], t[ok][-1], r1.min(), r1.max(), r2.min(), r2.max(), r3.min(), r3.max()))
assert r1.max() <= 1 + 1e-9 and r2.max() <= 1 + 1e-9 and r3.min() >= 1 - 1e-9
# int exp(-g): certified upper bound on the certificate mesh vs reference trapezoid on a fine uniform grid
tm = cert.choose_mesh(y, P2lo, CDlo, NU2hi)
gl, gpl, _ = cert.evalsums(tm, y, P2lo, CDlo, NU2hi)
tf = np.linspace(0, tm[-1] + 2, 801); gf, _ = ref_g(tf, X, W); e = np.exp(-gf)
ETref = (0.5 * (e[1:] + e[:-1]) * np.diff(tf)).sum()
print("E T(Q~_1): reference %.6f, certified upper bound %.6f" % (ETref, cert.mean_T_hi(tm, gl, gpl)))
assert cert.mean_T_hi(tm, gl, gpl) >= ETref - 1e-6
print("all direction checks passed")
