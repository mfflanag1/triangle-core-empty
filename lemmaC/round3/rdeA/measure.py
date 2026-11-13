"""FLOAT64 measurement of the constants of EXISTENCE.md (Theorem A).  Evidence, not proof:
the cell-wise worst-case structure is the one a rigorous run needs, but arithmetic is float64,
the within-cell defect is sampled at 3 points per cell, and the tail term uses the reference law
only (factor TAILSAFE) instead of the envelope T^+ of the theorem.
Reference law mu_0: piecewise-constant density, cells [Jh,(J+1)h), J = o..o+n-1, masses p_J."""
import numpy as np, sys
h = 0.01; A0, B0 = -6.0, 7.0; AL, BE = float(sys.argv[1]) if len(sys.argv) > 1 else -4.5, float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
TAILSAFE = 4.0
d = np.load('../rde/fixed_point_h005.npy'); xf, ff = d[0], d[1]; hf = xf[1] - xf[0]
o = int(round(A0 / h)); n = int(round((B0 - A0) / h)); J = o + np.arange(n)
edges = (o + np.arange(n + 1)) * h
k0 = int(round((A0 - xf[0]) / hf)); st = int(round(h / hf)); assert st == 2
p = hf * (ff[k0:k0 + st * n:st] / 2 + ff[k0 + 1:k0 + st * n:st] + ff[k0 + 2:k0 + st * n + 1:st] / 2); p /= p.sum()   # cell-wise trapezoid (no cancellation)
ja, jb = int(round(AL / h)) - o, int(round(BE / h)) - o          # bulk cells ja..jb-1
pB = p.copy(); pB[:ja] = 0; pB[jb:] = 0
cl = lambda s: np.clip(s, 0, None)
c3 = lambda s: (cl(s + 1) ** 3 - cl(s) ** 3) / 6                # E (s+V)_+^2/2, V~U[0,1]
c4 = lambda s: (cl(s + 2) ** 4 - 2 * cl(s + 1) ** 4 + cl(s) ** 4) / 24
c5 = lambda s: (cl(s + 2) ** 5 - 2 * cl(s + 1) ** 5 + cl(s) ** 5) / 120
pY = np.convolve(p, p); pYB = np.convolve(pB, pB); K = 2 * o + np.arange(2 * n - 1)
def Gs(I, py):                                                  # G0, G0' at x = I*h (I real array)
    S = I[:, None] + K[None, :]
    return h ** 3 * (c5(S) @ py), h ** 2 * (c4(S) @ py)
def g1(M, pp): return h ** 2 * (c3(M[:, None] + J[None, :]) @ pp)   # g1 at t = M*h
# ---- defect: b(u) - a(u) at cell edges and midpoints of bulk cells
F = np.concatenate([[0], np.cumsum(p)]); Fb = np.concatenate([np.cumsum(p[::-1])[::-1], [0]])
Ie = (o + np.arange(ja, jb + 1)).astype(float); Im = Ie[:-1] + 0.5
def disp(I, fbar):
    G, G1 = Gs(I, pY); return (-np.log(fbar) - G) / G1          # one Newton step from a(u)
de = disp(Ie, Fb[ja:jb + 1]); dm = disp(Im, Fb[ja:jb] - p[ja:jb] / 2)
EX0 = (p * (J + 0.5) * h).sum()
print('reference: h=%g support [%g,%g] bulk [%g,%g] cells %d  E X0 = %.10f  9EX0 = %.8f' % (h, A0, B0, AL, BE, jb - ja, EX0, 9 * EX0))
print('tail masses of mu_0 outside bulk: lower %.3e upper %.3e' % (F[ja], Fb[jb]))
print('defect displacement b-a: edges max|.| %.3e, midpoints max|.| %.3e' % (abs(de).max(), abs(dm).max()))
# ---- kernel envelopes on bulk x bulk
I = o + np.arange(ja, jb); JB = J[ja:jb]; pb = p[ja:jb]
Mall = np.arange(2 * I[0] - 2, 2 * I[-1] + 6); g = g1(Mall.astype(float), pB); off = Mall[0]
gmid = g1((I[:, None] + JB[None, :] + 1.0).ravel(), pB).reshape(len(I), len(JB))
Pm = pb[None, :] * gmid; Pm /= Pm.sum(1, keepdims=True)
Nlo = pb[None, :] * g[(I[:, None] + JB[None, :] - 1) - off]
Nhi = pb[None, :] * g[(I[:, None] + JB[None, :] + 3) - off]
Plo = Nlo / Nhi.sum(1, keepdims=True); Phi = Nhi / Nlo.sum(1, keepdims=True)
ev, V = np.linalg.eig(Pm.T); q = np.real(V[:, np.argmax(ev.real)]); q /= q.sum()
Amid = 2 * np.abs(Pm - q[None, :])
Aenv = 2 * np.maximum(Phi - q[None, :], q[None, :] - Plo)
THETA = float(sys.argv[3]) if len(sys.argv) > 3 else 0.7
print('specrad midpoint 2|P-q| = %.4f ; specrad envelope matrix = %.4f ; theta = %g' % (max(abs(np.linalg.eigvals(Amid))), max(abs(np.linalg.eigvals(Aenv))), THETA))
W = np.linalg.solve(np.eye(len(I)) - Aenv / THETA, np.ones(len(I))); assert W.min() >= 1
print('resolvent weight W = (1 - Aenv/theta)^-1 1: range [%.3f, %.1f]; at bulk edges %.2f / %.2f' % (W.min(), W.max(), W[0], W[-1]))
Zsum = Aenv @ W
Zi = Zsum / W; print('Z (cell-wise worst-case envelope, shift <= h) = max_i %.4f  (argmax x=%.2f); median %.4f' % (Zi.max(), I[np.argmax(Zi)] * h, np.median(Zi)))
# ---- tail term (reference only): (G0 - G0^BB)/G0'^BB at output points
Gall, _ = Gs(I.astype(float), pY); GBB, G1BB = Gs(I.astype(float), pYB)
tail0 = (Gall - GBB) / G1BB; tail = TAILSAFE * tail0
print('tail term tau0/E kappa (reference): max %.3e at x=%.2f; max of tail0/W %.3e' % (tail0.max(), I[np.argmax(tail0)] * h, (tail0 / W).max()))
# ---- radius
m0 = 0.0
Y = np.maximum(np.maximum(abs(de[:-1] - m0), abs(de[1:] - m0)), abs(dm - m0))
r = ((Y + tail) / (W - Zsum)).max(); i = np.argmax((Y + tail) / (W - Zsum))
print('max_i Y_i/W_i = %.3e ; radius r = max_i (Y_i + tail_i)/(W_i - Zsum_i) = %.3e (at x=%.2f)' % ((Y / W).max(), r, I[i] * h))
shift = 6 * r * W.max() + 2 * tail.max() + abs(de).max()
print('envelope shift needed 6 r Wmax + 2 tail + defect = %.3e  (must be <= h = %g): %s' % (shift, h, shift <= h))
costerr = r * ((pb * W).sum() + (2 / 3) * (q * W).sum())
print('sum p W = %.3f, sum q W = %.3f ; |E X* - E X0 - m0/3| <= %.3e + tail integrals ; 9x = %.3e' % ((pb * W).sum(), (q * W).sum(), costerr, 9 * costerr))
print('=> (float, conditional on tail envelopes) c in [%.5f, %.5f]' % (9 * (EX0 - costerr), 9 * (EX0 + costerr)))
