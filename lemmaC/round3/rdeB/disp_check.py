"""Locate apparent violations of slope monotonicity; test Theorem M on random IFR pairs."""
import numpy as np
from disp_iter import T, slope_at_level
L, h = 16.0, 0.01
n = int(round(2 * L / h)) + 1
x = -L + h * np.arange(n)
th = np.exp(np.linspace(np.log(1e-9), np.log(40.0), 300))
G = np.clip(x, 0, None) ** 3 / 6; G1 = np.clip(x, 0, None) ** 2 / 2
S = [slope_at_level(G, G1, th)]
f = G1 * np.exp(-G)
for it in range(1, 9):
    m = (x * f).sum() * h
    fc = np.interp(x + m, x, f, left=0, right=0)
    G, G1 = T(fc, x, h); f = G1 * np.exp(-G)
    S.append(slope_at_level(G, G1, th))
S = np.array(S)
for k in range(1, len(S)):
    r = (S[k] - S[k - 1]) / S[k - 1]
    i = np.argmax(r)
    print('step %d->%d: max rel change %+.3e at theta=%.3e ; min %+.3e ; #viol(>1e-6)=%d, theta range of viol: %s'
          % (k, k + 1, r[i], th[i], r.min(), (r > 1e-6).sum(),
             (th[r > 1e-6].min(), th[r > 1e-6].max()) if (r > 1e-6).any() else None))
