"""V10: accuracy of v10_exact.py's float64 quadrature, against the same integrals in mpmath (30 digits,
12-point Gauss-Legendre per subinterval, all piece tables and breakpoints in mp).  Run B, law 9, mesh 10."""
import sys, bisect
import numpy as np
from mpmath import mp, mpf, exp as mpexp
sys.path.insert(0, '/Users/maxf/projects/economics-research/game-theory/coopgpt/average_case_protection_2026-11-09/lemmaC/round3/v10')
import v10_exact
mp.dps = 30
pre = sys.argv[1]; n = int(sys.argv[2])
d0 = np.load("%s_stage%d.npz" % (pre, n)); d1 = np.load("%s_stage%d.npz" % (pre, n + 1))
theta, shat = d0['theta'], d0['shat']
L = v10_exact.Law(theta, shat)
N = len(theta) - 1
th = [mpf(float(v)) for v in theta]; sh = [mpf(float(v)) for v in shat]
b = [mpf(0)]; lev = []; sig = []
for k in range(N):
    i = N - k; lev.append(th[i]); sig.append(sh[i]); b.append(b[-1] + (th[i] - th[i - 1]) / sh[i])
lev.append(th[0]); sig.append(sh[0]); b.append(b[-1] + th[0] / sh[0]); lev.append(mpf(0)); sig.append(mpf(0))
m = [mpexp(-v) for v in lev]; M1 = [mpf(0)]; M2 = [mpf(0)]
for k in range(N + 1):
    dd = b[k + 1] - b[k]; I1 = (m[k + 1] - m[k]) / sig[k]
    M2.append(M2[k] + M1[k] * dd + (I1 - m[k] * dd) / sig[k]); M1.append(M1[k] + I1)
bf = [float(v) for v in b]
from mpmath import calculus
xs_w = mp.quad.__self__.get_nodes if False else None
import mpmath
nodes = mpmath.calculus.quadrature.GaussLegendre(mp).calc_nodes(3, mp.prec)   # degree 3: 3*2^2 = 12 nodes on [-1,1]


def piece(a):
    k = min(max(bisect.bisect_right(bf, float(a)) - 1, 0), N + 1)
    while k > 0 and b[k] > a: k -= 1
    while k < N + 1 and b[k + 1] <= a: k += 1
    return k


def M0(a):
    k = piece(a); return m[k] * mpexp(sig[k] * (a - b[k]))


def M12(u):
    k = piece(u); dd = u - b[k]
    if sig[k] == 0:
        return M1[k] + dd, M2[k] + M1[k] * dd + dd * dd / 2
    e = mpexp(sig[k] * dd); I1 = m[k] * (e - 1) / sig[k]
    return M1[k] + I1, M2[k] + M1[k] * dd + (I1 - m[k] * dd) / sig[k]


for i in [0, len(d1['t']) // 3, len(d1['t']) - 1]:
    t = mpf(float(d1['t'][i]))
    br = sorted(set([mpf(0), t] + [v for v in b if v < t] + [t - v for v in b if v < t]))
    g = mpf(0); gp = mpf(0)
    for lo, hi in zip(br[:-1], br[1:]):
        half = (hi - lo) / 2; mid = (hi + lo) / 2
        for xx, ww in nodes:
            a = mid + half * xx; m1, m2 = M12(t - a); m0 = M0(a) * ww * half
            g += m0 * m2; gp += m0 * m1
    gf, gpf = L.g_gp(float(t))
    print("point %d, t = %.6f: g mp %s  float64 rel.diff %.2e | g' mp %s  float64 rel.diff %.2e | theta/g - 1 = %s, shat/g' - 1 = %s"
          % (i, float(t), mp.nstr(g, 18), float(gf / g - 1), mp.nstr(gp, 18), float(gpf / gp - 1),
             mp.nstr(mpf(float(d1['theta'][i])) / g - 1, 6), mp.nstr(mpf(float(d1['shat'][i])) / gp - 1, 6)), flush=True)
