"""Route B test: monotone iteration in the dispersive order, started from delta_0.
Q_n = T^n(delta_0) modulo shift.  c_n = 3*Gamma(Q_n), Gamma(Q) = mean(T(Q)) + 2*mean(Q).
Claim under test: (i) slope-at-level s_n(theta) = G_n'(G_n^{-1}(theta)) is nonincreasing in n
for every theta (Q_n increasing in dispersive order); (ii) c_n decreases to 2.24954561.
Floating point: evidence only."""
import numpy as np, sys
from scipy.signal import fftconvolve

def T(f, x, h):
    fY = fftconvolve(f, f) * h
    y = 2 * x[0] + h * np.arange(len(fY))
    def rc(a): return np.concatenate([np.cumsum(a[::-1])[::-1], [0.0]])
    c0, c1, c2, c3 = rc(fY), rc(fY * y), rc(fY * y * y), rc(fY * y ** 3)
    j0 = np.searchsorted(y, -x - 1e-12, side='left')
    G = h * (x ** 3 * c0[j0] + 3 * x ** 2 * c1[j0] + 3 * x * c2[j0] + c3[j0]) / 6
    G1 = h * (x ** 2 * c0[j0] + 2 * x * c1[j0] + c2[j0]) / 2
    return np.maximum(G, 0), np.maximum(G1, 0)

def slope_at_level(G, G1, th):
    """s(theta) = G'(G^{-1}(theta)) by interpolation in log G."""
    ok = (G > 1e-300) & (G1 > 0)
    lg = np.log(G[ok]); return np.exp(np.interp(np.log(th), lg, np.log(G1[ok])))

if __name__ == '__main__':
    L, h = 16.0, 0.01
    n = int(round(2 * L / h)) + 1
    x = -L + h * np.arange(n)
    th = np.exp(np.linspace(np.log(1e-12), np.log(60.0), 400))
    # rho_1 = T(delta_0): G = x_+^3/6
    G = np.clip(x, 0, None) ** 3 / 6; G1 = np.clip(x, 0, None) ** 2 / 2
    s_prev = slope_at_level(G, G1, th)
    f = G1 * np.exp(-G); m_in = 0.0   # mean of the input delta_0
    from math import gamma
    print('n=0  c_0 = 3*E W = %.6f' % (3 * 6 ** (1 / 3) * gamma(4 / 3)))
    for it in range(1, 16):
        m = (x * f).sum() * h; mass = f.sum() * h
        sd = np.sqrt((x * x * f).sum() * h - m * m)
        # recentre f to mean 0 by interpolation shift (keeps it on the grid)
        fc = np.interp(x + m, x, f, left=0, right=0)
        G, G1 = T(fc, x, h)
        fn = G1 * np.exp(-G)
        mo = (x * fn).sum() * h
        s = slope_at_level(G, G1, th)
        viol = np.max((s - s_prev) / s_prev)        # should be <= 0 (slopes decrease)
        print('n=%d  sd(Q_n)=%.6f  c_n=3*Gamma(Q_n)=%.7f  max rel increase of s: %+.2e  mass %.10f'
              % (it, sd, 3 * mo, viol, mass), flush=True)
        s_prev = s; f = fn
