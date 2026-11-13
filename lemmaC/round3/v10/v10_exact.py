"""V10 independent recomputation (no code shared with cert.py).

Exact g, g' of the image of a rounded law by a different method:
  M0(a) = P(X > -a) (a > 0), M1 = int_0^u M0, M2 = int_0^u M1 (closed forms per piece),
  g(t)  = E(t+X+X')_+^3/6 = int_0^t M0(a) M2(t-a) da,   g'(t) = int_0^t M0(a) M1(t-a) da,
obtained from two layer-cake integrations by parts.  The t-integrals are done by 10-point
Gauss-Legendre on every subinterval between consecutive breakpoints of a -> M0(a) M2(t-a)
(so the integrand is analytic on each subinterval).  Piece tables are built with mpmath (40 digits).

Usage: v10_exact.py PREFIX LASTSTAGE STRIDE [final]
"""
import sys, time
import numpy as np
from mpmath import mp, mpf, exp as mpexp

mp.dps = 40
GLx, GLw = np.polynomial.legendre.leggauss(10)


class Law:
    def __init__(self, theta, shat):
        N = len(theta) - 1
        th = [mpf(float(v)) for v in theta]; sh = [mpf(float(v)) for v in shat]
        # pieces in a = -x, increasing a: piece k = [b_k, b_{k+1}], level at b_k = lev_k, exponent rate sig_k
        b = [mpf(0)]; lev = []; sig = []
        for k in range(N):                      # x in [x_{N-k-1}, x_{N-k}], slope shat_{N-k}
            i = N - k
            lev.append(th[i]); sig.append(sh[i]); b.append(b[-1] + (th[i] - th[i - 1]) / sh[i])
        lev.append(th[0]); sig.append(sh[0]); b.append(b[-1] + th[0] / sh[0])      # piece N: down to level 0
        lev.append(mpf(0)); sig.append(mpf(0))                                     # piece N+1: M0 = 1
        m = [mpexp(-v) for v in lev]
        M1 = [mpf(0)]; M2 = [mpf(0)]
        for k in range(N + 1):
            d = b[k + 1] - b[k]
            I1 = (m[k + 1] - m[k]) / sig[k]
            M2.append(M2[k] + M1[k] * d + (I1 - m[k] * d) / sig[k])
            M1.append(M1[k] + I1)
        self.amin = b[-1]
        self.EX = -(b[-1] - M1[-1])            # E X = -int_0^amin (1 - M0)
        f = lambda v: np.array([float(u) for u in v])
        self.b, self.sig, self.m, self.M1, self.M2 = f(b), f(sig), f(m), f(M1), f(M2)
        self.atom = float(m[0])

    def _idx(self, a):
        return np.clip(np.searchsorted(self.b, a, side='right') - 1, 0, len(self.b) - 1)

    def M0(self, a):
        k = self._idx(a)
        return self.m[k] * np.exp(self.sig[k] * (a - self.b[k]))

    def M12(self, u):
        k = self._idx(u); d = u - self.b[k]; z = self.sig[k] * d
        small = z < 1.0
        zs = np.where(small, z, 0.0)
        p1 = np.ones_like(z); p2 = np.full_like(z, 0.5); term1 = np.ones_like(z); term2 = np.full_like(z, 0.5)
        for j in range(1, 24):                  # phi1 = sum z^j/(j+1)!, phi2 = sum z^j/(j+2)!
            term1 = term1 * zs / (j + 1); term2 = term2 * zs / (j + 2); p1 += term1; p2 += term2
        zb = np.where(small, 1.0, z)
        p1 = np.where(small, p1, np.expm1(zb) / zb); p2 = np.where(small, p2, (np.expm1(zb) - zb) / (zb * zb))
        mk = self.m[k]
        return self.M1[k] + mk * d * p1, self.M2[k] + self.M1[k] * d + mk * d * d * p2

    def g_gp(self, t):
        """g(t), g'(t) for a scalar t > 0."""
        bb = self.b[self.b < t]
        br = np.unique(np.concatenate([bb, t - bb, [0.0, t]]))
        br = br[(br >= 0) & (br <= t)]
        lo, hi = br[:-1], br[1:]
        half = 0.5 * (hi - lo); mid = 0.5 * (hi + lo)
        a = mid[:, None] + half[:, None] * GLx[None, :]
        w = half[:, None] * GLw[None, :]
        m0 = self.M0(a); m1, m2 = self.M12(t - a)
        return float((w * m0 * m2).sum()), float((w * m0 * m1).sum())


def main():
    prefix, last, stride = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    final = len(sys.argv) > 4
    d = {n: np.load("%s_stage%d.npz" % (prefix, n)) for n in range(1, last + 1)}
    print("checks of theta_i <= g(t_i), shat_i >= g'(t_i): prefix %s" % prefix, flush=True)
    # stage 1 comes from delta_0: g = t^3/6, g' = t^2/2 exactly (rational)
    from fractions import Fraction
    t = d[1]['t']; th = d[1]['theta']; sh = d[1]['shat']
    r1 = max(Fraction(float(a)) / (Fraction(float(b)) ** 3 / 6) for a, b in zip(th, t))
    r2 = min(Fraction(float(a)) / (Fraction(float(b)) ** 2 / 2) for a, b in zip(sh, t))
    print(" stage 1 (exact rationals, all %d points): max theta/g - 1 = %.3e, min shat/g' - 1 = %.3e" % (len(t), float(r1 - 1), float(r2 - 1)), flush=True)
    for n in range(1, last):
        tic = time.time()
        L = Law(d[n]['theta'], d[n]['shat'])
        t = d[n + 1]['t']; th = d[n + 1]['theta']; sh = d[n + 1]['shat']
        assert np.all(np.diff(th) > 0) and np.all(np.diff(sh) >= 0) and th[0] > 0 and sh[0] > 0 and len(t) == len(th) == len(sh)
        idx = np.unique(np.concatenate([np.arange(0, len(t), stride), np.arange(0, min(15, len(t))), np.arange(max(0, len(t) - 15), len(t))]))
        G = np.array([L.g_gp(float(t[i])) for i in idx])
        r1 = th[idx] / G[:, 0] - 1; r2 = sh[idx] / G[:, 1] - 1
        print(" stage %d -> %d: %d of %d points; theta/g - 1 in [%.3e, %.3e] (need <= 0); shat/g' - 1 in [%.3e, %.3e] (need >= 0); "
              "law: atom %.3e, support %.6f, E X = %.9f  [%.0fs]"
              % (n, n + 1, len(idx), len(t), r1.min(), r1.max(), r2.min(), r2.max(), L.atom, float(L.amin), float(L.EX), time.time() - tic), flush=True)
    if final:
        for n in sorted(set([last] + [int(v) for v in sys.argv[5:]])):
            tic = time.time()
            L = Law(d[n]['theta'], d[n]['shat'])
            out = []
            for dt in (1 / 32.0, 1 / 64.0):
                tt = np.arange(0, 20 + dt / 2, dt)
                e = np.array([1.0] + [np.exp(-L.g_gp(float(v))[0]) for v in tt[1:]])
                wS = np.ones(len(tt)); wS[1:-1:2] = 4; wS[2:-1:2] = 2
                out.append(float((wS * e).sum() * dt / 3))
            ET = out[1] + (out[1] - out[0]) / 15
            print(" final law n=%d: int exp(-g) = %.10f (Simpson dt=1/32: %.10f, 1/64: %.10f); E Q~ = %.10f; c(Q~_%d) = %.8f  [%.0fs]"
                  % (n, ET, out[0], out[1], float(L.EX), n, 3 * (ET + 2 * float(L.EX)), time.time() - tic), flush=True)


if __name__ == '__main__':
    main()
