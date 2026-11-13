"""(1) one-sided rounding principle; (2) Lemma 4 super-solution aR and the decreasing iteration."""
import numpy as np, sys
from v9_core import *
G = Grid(panel=0.25, k=8); th = G.th; E = np.exp(G.edges); rng = np.random.default_rng(5)
def step(q):
    """exact T at nodes and at panel edges; returns centred T(q) at nodes, same shift at edges, slopes at edges, c = 3*Gamma(q)"""
    lv = np.concatenate([th, E]); t, s, P = T(q, G, levels=lv)
    tn, te, se = t[:len(th)], t[len(th):], s[len(th):]
    c = 3 * (mean(tn, G) + 2 * mean(q, G)); m = mean(tn, G)
    return tn - m, te - m, se, c
# exact iterates
q = np.zeros_like(th); exact = []; Qn = [q]; Qe = [np.zeros_like(E)]
for n in range(0, 9):
    tn, te, se, c = step(q); exact.append(c); q = tn; Qn.append(tn); Qe.append(te)
print('exact  c_1..c_8:', ' '.join('%.5f' % c for c in exact[1:9]))
def rounded(se, i0, i1, stride, up=False):
    """minorant: slope 1/s(theta_j) on (theta_{j-1},theta_j], 1/s(theta_0) below, 0 above. (up=True: the WRONG direction, slope 1/s(theta_{j-1}))"""
    idx = np.arange(i0, i1 + 1, stride); nodes = E[idx]; s = se[idx]
    cell = np.searchsorted(nodes, th, side='left')              # th in (nodes[cell-1], nodes[cell]]
    N = len(nodes) - 1
    if not up:
        sl = np.concatenate([[1 / s[0]], 1 / s[1:], [0.0]])       # slope in cell 0 (below theta_0), cells 1..N, cell N+1 (above)
    else:
        sl = np.concatenate([[1 / s[0]], 1 / s[:-1], [0.0]])
    left = np.concatenate([[0.0], nodes])                       # left end of each cell
    Qleft = np.concatenate([[0.0], np.cumsum(sl[:-1] * np.diff(left))])
    return Qleft[cell] + sl[cell] * (th - left[cell])
for (wlo, whi, stride) in [(-14, 3.5, 4), (-14, 3.5, 2), (-14, 3.5, 1), (-9, 3.0, 1)]:
    i0 = int(round((wlo - G.edges[0]) / 0.25)); i1 = int(round((whi - G.edges[0]) / 0.25))
    for up in (False, True):
        q = np.zeros_like(th); out = []; okdisp = True
        for n in range(0, 7):
            tn, te, se, c = step(q); out.append(c)
            q = rounded(se, i0, i1, stride, up)
            if not up: d = np.diff(Qn[n + 1] - q); okdisp &= bool((d >= -1e-9).all())   # rounded iterate <=_disp exact iterate
        print('%s mesh w in [%g,%g], %3d cells: c~_1..c~_6 = %s | all >= exact: %s%s' % ('WRONG-way' if up else 'minorant ', wlo, whi, (i1 - i0) // stride,
              ' '.join('%.5f' % c for c in out[1:7]), all(out[i] >= exact[i] - 1e-9 for i in range(1, 7)), '' if up else ' | Q~_n <=_disp Q_n: %s' % okdisp))
# random less-dispersed replacements: Q~' = Q_n' * phi, phi nonincreasing step function in [0,1] (constant on panels)
bad = 0; tot = 0; mn = 9
for trial in range(12):
    n = rng.integers(1, 5); phi = np.sort(rng.uniform(0, 1, len(E) - 1))[::-1]
    if rng.random() < 0.5: phi[rng.integers(100, 159):] = 0.0            # atom at the top
    pan = np.repeat(np.arange(len(E) - 1), len(th) // (len(E) - 1))
    inc = np.concatenate([[0.0], np.cumsum(phi * np.diff(Qe[n]))])       # Q~ at edges
    q = inc[pan] + phi[pan] * (Qn[n] - Qe[n][pan])
    for m in range(3):
        tn, te, se, c = step(q); q = tn; tot += 1; bad += c < exact[n + m] - 1e-9; mn = min(mn, c - exact[n + m])
print('random minorants Q~ <=_disp Q_n (12 trials x 3 further steps): 3Gamma(T^m Q~) >= c_{n+m} fails in %d/%d, min excess %.3e' % (bad, tot, mn))
# Lemma 4
R = np.where(th <= 1, np.log(th), 3 * (th ** (1 / 3) - 1))
lev = np.exp(np.linspace(np.log(1e-12), np.log(40), 800))
t, s, P = T(R, G, levels=lev); ell = lev / s
i1 = np.argmax(ell > 1); print('Lemma 4: l_TR(theta) at 1e-12,1e-9,1e-6,1e-3,1: ', ' '.join('%.4f' % ell[np.argmin(abs(lev - v))] for v in [1e-12, 1e-9, 1e-6, 1e-3, 1]),
      '; first theta with l_TR>1: %.3f ; A = sup l_TR/theta^(1/3) beyond = %.4f' % (lev[i1], (ell / lev ** (1 / 3))[i1:].max()))
a = max(1.0, (ell / lev ** (1 / 3))[i1:].max())
for aa in [a, 1.0, 0.8, 0.7, 0.6]:
    t, s, P = T(aa * R, G, levels=lev); ratio = (aa * np.where(lev <= 1, 1 / lev, lev ** (-2 / 3))) * s      # (aR)'/T(aR)' must be >= 1
    q = aa * R; cs = []
    for m in range(8): tn, te, se, c = step(q); cs.append(c); q = tn
    print('  a=%.4f: min (aR)\'/T(aR)\' = %.4f at theta=%.2e ; 3Gamma(T^k aR), k=0..7: %s' % (aa, ratio.min(), lev[np.argmin(ratio)], ' '.join('%.5f' % c for c in cs)))
