"""Deterministic solver for G(x) = E[(x+X'+X'')_+^3]/6, P(X>x)=exp(-G(x)).
Unknown: density f on a uniform grid (trapezoid on the line: spectrally accurate for the
convolution, O(h^4) for the kinked kernels). Iteration: mixture damping 1/3, which kills
the shift mode (multiplier 1-3*lam = 0).  Floating point: evidence only."""
import numpy as np, sys
from scipy.signal import fftconvolve

def T(f, x, h):
    """return G, G' on grid x for input density f on grid x (x symmetric, x[0] = -L)."""
    fY = fftconvolve(f, f) * h                      # density of Y on y = 2*x[0] + h*j
    y = 2 * x[0] + h * np.arange(len(fY))
    # tail sums over y_j > -x_i ; kernels vanish at y=-x so side does not matter
    def rc(a): return np.concatenate([np.cumsum(a[::-1])[::-1], [0.0]])
    c0, c1, c2, c3 = rc(fY), rc(fY * y), rc(fY * y * y), rc(fY * y ** 3)
    j0 = np.searchsorted(y, -x - 1e-12, side='left')
    G = h * (x ** 3 * c0[j0] + 3 * x ** 2 * c1[j0] + 3 * x * c2[j0] + c3[j0]) / 6
    G1 = h * (x ** 2 * c0[j0] + 2 * x * c1[j0] + c2[j0]) / 2
    return np.maximum(G, 0), np.maximum(G1, 0), fY, y

def solve(L=14.0, h=0.02, iters=4000, lam=1 / 3, tol=1e-14):
    n = int(round(2 * L / h)) + 1
    x = -L + h * np.arange(n)
    f = np.exp(-(x - 0.25) ** 2 / (2 * 0.75)); f /= f.sum() * h
    for it in range(iters):
        G, G1, fY, y = T(f, x, h)
        fn = G1 * np.exp(-G)
        d = np.abs(fn - f).sum() * h
        f = (1 - lam) * f + lam * fn
        if d < tol: break
    G, G1, fY, y = T(f, x, h)
    fS = fftconvolve(fY, f) * h; s = 3 * x[0] + h * np.arange(len(fS)); sp = np.clip(s, 0, None)
    EX = (x * f).sum() * h
    out = dict(h=h, L=L, iters=it + 1, resid=d, mass=f.sum() * h, EX=EX,
               sd=np.sqrt((x * x * f).sum() * h - EX ** 2),
               norm=(fS * sp ** 3).sum() * h / 6, c_direct=(fS * sp ** 4).sum() * h / 8,
               c_9EX=9 * EX, EW=(fS * sp ** 4).sum() * h / 18)
    return out, x, f, G, G1

if __name__ == '__main__':
    res = []
    for h in [0.04, 0.02, 0.01, 0.005]:
        o, x, f, G, G1 = solve(h=h)
        res.append(o)
        print({k: (float('%.12g' % v) if isinstance(v, float) else v) for k, v in o.items()}, flush=True)
    for k in ['EX', 'c_direct', 'c_9EX', 'norm', 'sd']:
        v = [r[k] for r in res]
        # Richardson assuming h^4 and observed order
        d1, d2 = v[-2] - v[-3], v[-1] - v[-2]
        order = np.log2(abs(d1 / d2)) if d2 != 0 else float('nan')
        print(k, 'h=.005: %.12f' % v[-1], 'diffs %.3e %.3e' % (d1, d2), 'obs order %.2f' % order,
              'rich(h^4) %.12f' % (v[-1] + d2 / 15))
    np.save(__file__.replace('rde_solve.py', 'fixed_point_h005.npy'), np.vstack([x, f, G, G1]))
