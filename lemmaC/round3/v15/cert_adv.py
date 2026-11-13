"""v15: adversarial checks of cert.py pieces against independent high-precision references."""
import sys, numpy as np, mpmath as mp
from mpmath import iv
sys.argv = ['cert.py', '6', '10', '10', '1']
sys.path.insert(0, '../cert')
import cert
# 1. comparison semantics
a = iv.mpf([1, 3]); b = iv.mpf([2, 4]); c = iv.mpf([5, 6])
print("iv compare overlap  a<b :", a < b, " a<=b:", a <= b, " a>=b:", a >= b, " | disjoint a<c:", a < c, " c<a:", c < a, " | point 0<0:", iv.mpf(0) < iv.mpf(0), " touching [1,2]<[2,3]:", iv.mpf([1,2]) < iv.mpf([2,3]), " [1,2]<=[2,3]:", iv.mpf([1,2]) <= iv.mpf([2,3]))
# 2. law_to_lattice vs exact (mp, 60 digits)
mp.mp.dps = 50
rng = np.random.default_rng(7)
worst = dict(p=0, D=0, nu=0)
for trial in range(6):
    N = int(rng.integers(3, 12))
    theta = np.sort(rng.uniform(0.01, 6, size=N + 1)); theta[0] = rng.uniform(1e-4, 0.05)
    theta = np.sort(theta)
    shat = np.sort(rng.uniform(0.05, 6, size=N + 1))
    Kmin, p_lo, D_lo, nu_hi, EX_hi, width = cert.law_to_lattice(theta, shat)
    h = mp.mpf(cert.hf)
    th = [mp.mpf(float(v)) for v in theta]; sh = [mp.mpf(float(v)) for v in shat]
    x = [None] * (N + 1); x[N] = mp.mpf(0)
    for i in range(N, 0, -1): x[i - 1] = x[i] - (th[i] - th[i - 1]) / sh[i]
    xmin = x[0] - th[0] / sh[0]
    def Gt(z):
        if z <= xmin: return mp.mpf(0)
        if z <= x[0]: return sh[0] * (z - xmin)
        for i in range(1, N + 1):
            if z <= x[i]: return th[i] - sh[i] * (x[i] - z)
        return mp.inf
    S = lambda z: mp.e ** (-Gt(z)) if z < 0 else mp.mpf(0)
    Sleft0 = mp.e ** (-th[N])
    bp = [xmin] + x
    ncell = -Kmin
    Sk = [S(mp.mpf(Kmin + c_) * h) for c_ in range(ncell)] + [Sleft0]
    J = []
    for c_ in range(ncell):
        l = mp.mpf(Kmin + c_) * h; r = l + h
        pts = [l] + [z for z in bp if l < z < r] + [r]
        val = mp.mpf(0)
        for u, w in zip(pts[:-1], pts[1:]):
            # closed form on a linear piece
            mid = (u + w) / 2
            if mid <= xmin: val += w - u
            else:
                # slope at mid
                if mid <= x[0]: s = sh[0]
                else:
                    for i in range(1, N + 1):
                        if mid <= x[i]: s = sh[i]; break
                val += (mp.e ** (-Gt(u)) - (mp.e ** (-Gt(w)) if w < 0 else Sleft0)) / s
        J.append(val / h)
    EX = mp.mpf(Kmin) * h + h * sum(J)
    ok = True
    for c_ in range(ncell):
        p = Sk[c_] - Sk[c_ + 1]; D = h * (J[c_] - Sk[c_ + 1]); nu = (mp.mpf(1) if c_ == 0 else J[c_ - 1]) - J[c_]
        ok &= (mp.mpf(float(p_lo[c_])) <= p) and (mp.mpf(float(D_lo[c_])) <= D) and (mp.mpf(float(nu_hi[c_])) >= nu)
        if p > 0: worst['p'] = max(worst['p'], float((p - p_lo[c_]) / p)) if p > 1e-250 else worst['p']
        if nu > 1e-250: worst['nu'] = max(worst['nu'], float((nu_hi[c_] - nu) / nu))
        if D > 1e-250: worst['D'] = max(worst['D'], float((D - D_lo[c_]) / D))
    ok &= (mp.mpf(float(p_lo[ncell])) <= Sleft0) and (mp.mpf(float(nu_hi[ncell])) >= J[-1]) and (mp.mpf(float(EX_hi)) >= EX) and D_lo[ncell] == 0
    print("trial %d N=%d cells=%d one-sided: %s   EX_hi-EX=%.2e  sum p_lo=%.15f sum nu_hi=%.15f" % (trial, N, ncell, ok, float(EX_hi - EX), p_lo.sum(), nu_hi.sum()))
    # 3. evalsums vs independent tensor Gauss-Legendre of the exact law
    nodes, wts = np.polynomial.legendre.leggauss(60)
    X = [0.0]; W = [float(Sleft0)]
    bpf = [float(z) for z in bp]
    for j in range(1, len(bpf)):
        u, w = bpf[j - 1], bpf[j]
        s = float(sh[j - 1])
        xs = (u + w) / 2 + (w - u) / 2 * nodes
        Gs = np.array([float(Gt(mp.mpf(float(z)))) for z in xs])
        X += list(xs); W += list((w - u) / 2 * wts * s * np.exp(-Gs))
    X = np.array(X); W = np.array(W)
    Y = (X[:, None] + X[None, :]).ravel(); WW = (W[:, None] * W[None, :]).ravel()
    P2lo = cert.lo(cert.conv_nonneg(p_lo, p_lo)); CDlo = cert.lo(cert.conv_nonneg(D_lo, p_lo)); NU2hi = cert.hi(cert.conv_nonneg(nu_hi, nu_hi))
    y = (2 * Kmin + np.arange(len(P2lo))) * cert.hf
    t = np.round(np.linspace(0.02, min(25.0, 3 - 2 * float(xmin)), 60) * cert.TGRID) / cert.TGRID
    glo, gplo, gphi = cert.evalsums(t, y, P2lo, CDlo, NU2hi)
    U = np.maximum(t[:, None] + Y[None, :], 0)
    g = (U ** 3 @ WW) / 6; gp = (U ** 2 @ WW) / 2
    print("   mass(quad)=%.12f  min (g-glo)/g = %.3e  min (gp-gplo)/gp = %.3e  min (gphi-gp)/gp = %.3e   max: %.2e %.2e %.2e" % (W.sum(), ((g - glo) / g).min(), ((gp - gplo) / gp).min(), ((gphi - gp) / gp).min(), ((g - glo) / g).max(), ((gp - gplo) / gp).max(), ((gphi - gp) / gp).max()))
    # R5 bound vs numeric integral of exp(-g)
    tt = np.linspace(0, 60, 200001)
    # crude exact: g on fine grid by interpolation of exact tensor evaluation at 400 points is too slow; use bound check only at mesh: mean_T_hi >= trapezoid of exp(-g) on dense grid
    td = np.linspace(1e-6, max(30.0, t[-1] + 10), 3000)
    Ud = np.maximum(td[:, None] + Y[None, :], 0); gd = (Ud ** 3 @ WW) / 6
    ET = td[0] + np.trapz(np.exp(-gd), td)
    keep = np.concatenate([[True], np.diff(glo) > 0]) & (glo > 0)
    print("   E T (trapz) = %.8f   mean_T_hi(mesh of 60) = %.8f   (must be >=)" % (ET, cert.mean_T_hi(t[keep], glo[keep], gplo[keep])))
print("worst relative slack: ", worst)
# rexp vs mp
mp.mp.dps = 60; bad = 0
for _ in range(3000):
    xv = float(rng.uniform(-700, 5)) * float(rng.choice([1, 1e-3, 1e-8]))
    e = cert.rexp(iv.mpf(xv)); ref = mp.e ** mp.mpf(xv)
    if not (mp.mpf(e.a.a if hasattr(e.a,'a') else e.a) <= ref): bad += 1
    if not (ref <= mp.mpf(e.b.b if hasattr(e.b,'b') else e.b)): bad += 1
print("rexp enclosure failures (3000 args):", bad)
