"""Feasibility of a certified c_n: iterate with ONE-SIDED ROUNDING at every stage.
After each T, replace Q by the minorant Qt (in the dispersive order, still concave):
slope of Qt on the cell (th_j, th_{j+1}] = 1/s(th_{j+1}) (s is nondecreasing in theta), slope 1/s(th_0) below th_0,
slope 0 above th_N (atom at the top).  By Theorem M + Lemma 3 the resulting costs are >= the true c_n >= c(fixed point).
Floating point here: evidence that a coarse mesh already gives < 2.30; not a certificate."""
import numpy as np, sys
from disp_iter import T, slope_at_level
L, h = 16.0, 0.005
n = int(round(2 * L / h)) + 1
x = -L + h * np.arange(n)
tg = np.exp(np.linspace(np.log(1e-15), np.log(200.0), 60001))
def run(N, tmin=1e-6, tmax=30.0, steps=8):
    nodes = np.exp(np.linspace(np.log(tmin), np.log(tmax), N + 1))
    f = None; out = []
    G = np.clip(x, 0, None) ** 3 / 6; Gp = np.clip(x, 0, None) ** 2 / 2      # T(delta_0)
    for it in range(1, steps + 1):
        s = slope_at_level(G, Gp, nodes)
        cell = np.searchsorted(nodes, tg, side='left')                     # tg in (nodes[cell-1], nodes[cell]]
        sl = np.where(cell == 0, 1 / s[0], np.where(cell > N, 0.0, 1 / s[np.minimum(cell, N)]))
        Q = np.concatenate([[0], np.cumsum(0.5 * (sl[1:] + sl[:-1]) * np.diff(tg))])
        Q = Q + 1e-9 * np.arange(len(Q)) / len(Q)                          # strictly increasing for interp (atom smeared by 1e-9)
        lG = np.interp(x + Q[len(Q) // 2], Q, np.log(tg), left=-800, right=800)
        F = -np.expm1(-np.exp(np.minimum(lG, 700))); f = np.clip(np.gradient(F, h), 0, None); f /= f.sum() * h
        m = (x * f).sum() * h; fc = np.interp(x + m, x, f, left=0, right=0)
        G, Gp = T(fc, x, h); fo = Gp * np.exp(-G)
        out.append(3 * (x * fo).sum() * h)
    return out
for N in [50, 100, 200, 400]:
    print('N=%4d cells: rounded c_1..c_8 =' % N, ' '.join('%.4f' % c for c in run(N)), flush=True)
