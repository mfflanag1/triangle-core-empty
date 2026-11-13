"""c_1 by a third method: x-grid, no recentring/interpolation (W>=0 has density w^2/2 exp(-w^3/6))."""
import numpy as np
from math import gamma
for h in [0.004, 0.002, 0.001]:
    w = h * np.arange(int(8 / h) + 1); f = w * w / 2 * np.exp(-w ** 3 / 6)
    fY = np.convolve(f, f) * h; y = h * np.arange(len(fY))            # density of Y on [0,16]
    t = -16 + h * np.arange(int(26 / h) + 1)
    suf = lambda a: np.concatenate([np.cumsum(a[::-1])[::-1], [0.0]])
    S0, S1, S2, S3 = suf(fY) * h, suf(fY * y) * h, suf(fY * y * y) * h, suf(fY * y ** 3) * h
    j = np.searchsorted(y, -t, side='right')
    g = np.maximum((t ** 3 * S0[j] + 3 * t * t * S1[j] + 3 * t * S2[j] + S3[j]) / 6, 0)
    sv = np.exp(-g)
    ET = sv[t >= 0].sum() * h - (1 - sv[t < 0]).sum() * h - 0.5 * h * sv[t >= 0][0] + 0.5*h*(1-sv[t<0][-1])*0
    EQ1 = 6 ** (1 / 3) * gamma(4 / 3)
    print('h=%.3f  E T(Q1) = %.6f  E Q1 = %.6f  c_1 = %.6f' % (h, ET, EQ1, 3 * (ET + 2 * EQ1)))
