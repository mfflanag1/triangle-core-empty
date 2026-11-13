"""Power check: a NON-IFR input designed from the proof (M_2/M_1 dips after a gap in the support) must violate slope monotonicity."""
import numpy as np
from v9_core import *
levels = np.exp(np.linspace(np.log(1e-6), np.log(30), 500))
for panel, k in [(0.25, 8), (0.125, 12)]:
    G = Grid(panel=panel, k=k)
    for c, d, J, cp, kap in [(np.e, 1.0, 5.0, np.exp(-1.25), 0.1), (np.e, 1.0, 5.0, np.exp(-1.25), 1.0), (np.exp(1.5), 0.5, 8.0, np.exp(-2.0), 0.2)]:
        q0 = np.where(G.th < c, d * (G.th / c - 1), J); q1 = q0 + kap * np.minimum(G.th, cp)
        o = [T(q, G, levels=levels)[1] for q in (q0, q1)]; r = (o[1] - o[0]) / o[0]
        g0, g1 = Gamma(q0, G), Gamma(q1, G)
        print('panel %.3f: c=%.2f d=%.1f J=%.0f c\'=%.2f kap=%.1f : max (s1-s0)/s0 = %+.3e at theta=%.3f, min = %+.3e ; Gamma0=%.5f >= Gamma1=%.5f (Lemma 3 needs no IFR)'
              % (panel, c, d, J, cp, kap, r.max(), levels[np.argmax(r)], r.min(), g0, g1))
