"""Recompute c_n = 3*Gamma(T^n delta_0), two resolutions; check RDE 1.5 identities at the converged point."""
import numpy as np, sys, time
from v9_core import *
ref = {1: 3.4152, 2: 2.6601, 3: 2.3980, 4: 2.3016, 5: 2.2678, 6: 2.2558}
for panel, k in [(0.5, 6)]:
    G = Grid(panel=panel, k=k); t0 = time.time()
    q = np.zeros_like(G.th); out = []
    for n in range(0, 41):
        t, s, P = T(q, G)
        c = 3 * (mean(t, G) + 2 * mean(q, G)); out.append(c)
        kshift = mean(t, G) - mean(q, G)          # T(q) ~ q + kshift at convergence (q centred)
        qprev = q; q = t - mean(t, G)
    print('panel %.3f k=%d nodes=%d (%.1fs)' % (panel, k, len(G.th), time.time() - t0))
    print('  c_0..c_10:', ' '.join('%.5f' % c for c in out[:11]))
    print('  claimed    ', ' '.join('%d:%.4f(diff %+.1e)' % (n, ref[n], out[n] - ref[n]) for n in ref))
    print('  c_20=%.8f c_40=%.8f ; ratios of successive differences n=6..10:' % (out[20], out[40]),
          ' '.join('%.4f' % ((out[n + 1] - out[n + 2]) / (out[n] - out[n + 1])) for n in range(6, 11)))
    print('  monotone decreasing c_n:', all(out[i] >= out[i + 1] - 1e-12 for i in range(40)))
    # true fixed point: qstar = qprev(centred) + kshift/3, check T(qstar)=qstar and RDE identities
    qs = qprev + kshift / 3; ts, ss, Ps = T(qs, G)
    EX = mean(qs, G); u = qs + 2 * Ps.shift; ES4 = 24 * float((G.p * Ps.M(4, u)).sum()); ES3 = 6 * float((G.p * Ps.M(3, u)).sum())
    print('  fixed point: max|T(q*)-q*| on theta in [1e-8,30] = %.2e ; E X = %.9f ; 9EX = %.8f ; E S_+^4/8 = %.8f ; E S_+^3 = %.8f'
          % (np.abs(ts - qs)[(G.th > 1e-8) & (G.th < 30)].max(), EX, 9 * EX, ES4 / 8, ES3))
    # direct Mecke-type checks: E[(X2+X3) 1{xi<S}] integral = E S_+^4/9 ; via E[X G(X)] = E S_+^4/18
    print('  E[X*G(X)] = %.8f  vs  E S_+^4/18 = %.8f' % (float((G.p * qs * Ps.M(3, u)).sum()), ES4 / 18))
