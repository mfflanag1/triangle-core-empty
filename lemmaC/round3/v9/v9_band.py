"""Investigate the largest apparent violation in the 'band' family by grid refinement; and a non-IFR control with power."""
import numpy as np, sys
from v9_core import *
def F(b, s): return np.log(s) if abs(b - 1) < 1e-12 else (s ** (1 - b) - 1) / (1 - b)
def run(G, terms, kap, c1, c2, levels):
    Q = lambda extra_c: sum(a * F(b, np.minimum(G.th, c)) for a, b, c in terms) + kap * np.minimum(G.th, extra_c)
    out = []
    for q in (Q(c1), Q(c2)):
        t, s, P = T(q, G, levels=levels); out.append((s, P.last_cond))
    return (out[1][0] - out[0][0]) / out[0][0], np.maximum(out[0][1], out[1][1])
rng = np.random.default_rng(1)
G0 = Grid(panel=0.25, k=6); e0 = np.exp(G0.edges)
levels = np.exp(np.linspace(np.log(1e-9), np.log(30), 600))
worst = []
for trial in range(60):
    terms = []
    for _ in range(rng.integers(1, 4)):
        a = np.exp(rng.uniform(-3, 2)); b = rng.choice([0.0, 1.0, rng.uniform(0, 1), 2 / 3])
        c = np.inf if rng.random() < 0.5 else e0[rng.integers(60, len(e0) - 4)]
        terms.append((a, b, c))
    i1 = rng.integers(90, len(e0) - 10); c1 = e0[i1]; c2 = e0[i1 + rng.integers(1, 6)]; kap = np.exp(rng.uniform(-3, 2))
    r, cond = run(G0, terms, kap, c1, c2, levels); r = np.where(cond < 1e7, r, -9); i = np.argmax(r)
    worst.append((r[i], levels[i], terms, kap, c1, c2))
worst.sort(key=lambda z: -z[0])
print('coarse grid, worst 3 of 60 band pairs:')
for z in worst[:3]:
    print('  max rel %+.3e at theta=%.2e ; kap=%.3g c1=%.3g c2=%.3g terms=%s' % (z[0], z[1], z[3], z[4], z[5], [(round(a, 3), round(b, 3), c) for a, b, c in z[2]]))
    for panel, k in [(0.25, 6), (0.25, 12), (0.125, 12)]:
        G = Grid(panel=panel, k=k); r, cond = run(G, z[2], z[3], z[4], z[5], levels); r = np.where(cond < 1e7, r, -9)
        j = np.argmin(abs(levels - z[1])); print('     panel %.3f k=%2d: max rel %+.3e (at theta=%.2e); at the coarse-worst level %+.3e' % (panel, k, r.max(), levels[np.argmax(r)], r[j]))
# control with power: two-atom-like non-IFR input (upper atom mass e^-c < 1/2), dispersive increase below the gap
G = Grid(panel=0.25, k=8); e = np.exp(G.edges); viol = 0; w = -9
for trial in range(30):
    c = e[rng.integers(143, 152)]; J = np.exp(rng.uniform(0, 2))
    q0 = J * (G.th > c) + 1e-3 * np.minimum(G.th, c); q1 = q0 + np.exp(rng.uniform(-3, 0)) * np.minimum(G.th, c)
    o = []
    for q in (q0, q1):
        t, s, P = T(q, G, levels=levels); o.append(s)
    r = np.where(P.last_cond < 1e7, (o[1] - o[0]) / o[0], -9); w = max(w, r.max()); viol += r.max() > 1e-6
print('non-IFR control (gap in support, top mass < 1/2): %d/30 violate, max (s1-s0)/s0 = %+.3e' % (viol, w))
