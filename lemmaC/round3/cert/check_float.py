"""Independent plain-float re-computation of the rounded scheme of cert.py (sanity check, NOT a proof).
Different code path: cells of width DX = 0.0007 (not dyadic), mass placed at cell MIDPOINTS, FFT
convolution, trapezoid rule for int exp(-g).  No directed rounding anywhere.

Usage: check_float.py PREFIX STAGES      (reads PREFIX_stage{n}.npz written by cert.py)
 (B) for each stage n: take the certified law (theta^(n), shat^(n)) and compute in floats
     c(Q~_n) = 3 (E T(Q~_n) + 2 E Q~_n), and g, g' at the next mesh t^(n+1); compare with the certified
     theta^(n+1) (must be <= g) and shat^(n+1) (must be >= g').
 (A) from scratch: iterate the scheme in floats (theta = g(t_i), s = g'(t_i) on the same meshes t^(n)),
     never reading the certified theta, shat."""
import sys, numpy as np
from scipy.signal import fftconvolve
DX = 0.0007
PREFIX, STAGES = sys.argv[1], int(sys.argv[2])


def law(theta, shat):
    """midpoint masses on the half-lattice: a[r] = mass at position -r*DX/2; returns a, E X."""
    if theta is None:
        return np.array([1.0]), 0.0
    dxs = np.diff(theta) / shat[1:]
    x = np.concatenate([[0.0], -np.cumsum(dxs[::-1])])[::-1]          # x_0..x_N, x_N = 0
    xb = np.concatenate([[x[0] - theta[0] / shat[0]], x]); lv = np.concatenate([[0.0], theta])
    M = int(np.ceil(-xb[0] / DX)) + 1
    edges = -DX * np.arange(M + 1)                                     # 0, -DX, ...
    S = np.exp(-np.interp(edges, xb, lv, left=0.0))                    # S(0-) = exp(-theta_N)
    mass = S[1:] - S[:-1]                                              # cell [-(j+1)DX, -jDX]
    a = np.zeros(2 * M + 1); a[0] = S[0]; a[1::2] = mass
    return a, -(a * np.arange(2 * M + 1)).sum() * DX / 2


def g_at(t, a, chunk=64):
    P = np.clip(fftconvolve(a, a), 0, None); y = -np.arange(len(P)) * DX / 2
    g = np.empty(len(t)); gp = np.empty(len(t))
    for s in range(0, len(t), chunk):
        tt = t[s:s + chunk]; m1 = np.searchsorted(-y, tt.max()) + 1
        U = np.maximum(tt[:, None] + y[None, :m1], 0.0); U2 = U * U
        g[s:s + chunk] = (U2 * U) @ P[:m1] / 6; gp[s:s + chunk] = U2 @ P[:m1] / 2
    return g, gp


def mean_T(t, g):
    e = np.exp(-g)
    return t[0] + (0.5 * (e[1:] + e[:-1]) * np.diff(t)).sum()


d = {n: np.load("%s_stage%d.npz" % (PREFIX, n)) for n in range(1, STAGES + 1)}
print("(B) certified laws re-evaluated in floats")
for n in range(0, STAGES):
    th, sh = (None, None) if n == 0 else (d[n]['theta'], d[n]['shat'])
    a, EX = law(th, sh)
    t = d[n + 1]['t']; g, gp = g_at(t, a)
    c = 3 * (mean_T(t, g) + 2 * EX)
    r1 = d[n + 1]['theta'] / g; r2 = d[n + 1]['shat'] / gp
    print(" stage %d: float c(Q~_%d) = %.6f | next mesh: theta_cert/g_float in [%.6f, %.6f] (want <= 1), "
          "shat_cert/g'_float in [%.6f, %.6f] (want >= 1)" % (n, n, c, r1.min(), r1.max(), r2.min(), r2.max()), flush=True)
a, EX = law(d[STAGES]['theta'], d[STAGES]['shat']); t = d[STAGES]['t']; g, gp = g_at(t, a)   # last law: reuse the last mesh for the trapezoid rule
print(" stage %d: float c(Q~_%d) = %.6f  (g on the stage-%d mesh spans [%.2g, %.3g])" % (STAGES, STAGES, 3 * (mean_T(t, g) + 2 * EX), STAGES, g[0], g[-1]))
print("(A) scheme iterated from scratch in floats on the same meshes")
th = sh = None
for n in range(0, STAGES):
    a, EX = law(th, sh)
    t = d[n + 1]['t']; g, gp = g_at(t, a)
    print(" stage %d: float c = %.6f" % (n, 3 * (mean_T(t, g) + 2 * EX)), flush=True)
    th, sh = g, np.maximum.accumulate(gp)
a, EX = law(th, sh); g, gp = g_at(t, a)
print(" stage %d: float c = %.6f" % (STAGES, 3 * (mean_T(t, g) + 2 * EX)))
