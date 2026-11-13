"""Float exploration: weighted Dobrushin (sup-norm) and Hilbert-Schmidt constants of the
linearised shape map at the numerical fixed point.  Evidence only."""
import numpy as np
d = np.load('../rde/fixed_point_h005.npy'); x0, f0, G0, G10 = d
st = 8; x = x0[::st]; f = f0[::st]; G = G0[::st]; G1 = G10[::st]; h = x[1] - x[0]
n = len(x); print('h', h, 'n', n, 'mass', f.sum() * h)
Fbar = np.exp(-G)
for xx in [-6, -5, -4, -3.5, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6]:
    i = np.argmin(abs(x - xx)); print('x %5.1f  F %.3e  Fbar %.3e  f %.3e G1 %.3e' % (x[i], 1 - Fbar[i] if xx > -2 else G[i], Fbar[i], f[i], G1[i]))
t = 2 * x[0] + h * np.arange(2 * n - 1)
g1 = np.array([(np.clip(tt + x, 0, None) ** 2 * f).sum() * h / 2 for tt in t])
def consts(lo, hi):
    m = (x >= lo) & (x <= hi); idx = np.where(m)[0]
    K = g1[np.add.outer(idx, idx)]
    P = h * f[None, idx] * K; P /= P.sum(1, keepdims=True)       # rows: output x_i, cols: input a_j
    wu = h * f[idx]; wu = wu / wu.sum()                           # du measure on output (fixed pt: same law)
    res = {}
    for name, q in [('f', wu), ('pbar', wu @ P), ('pi', None)]:
        if q is None:
            ev, V = np.linalg.eig(P.T); q = np.real(V[:, np.argmax(ev.real)]); q = q / q.sum()
        A = 2 * np.abs(P - q[None, :])
        rho = max(abs(np.linalg.eigvals(A)))
        # HS in L2(du): kernel density p(u,v) = P_ij / wu_j
        C = (P - q[None, :]) / wu[None, :]
        hs = np.sqrt((wu[:, None] * C ** 2 * wu[None, :]).sum())
        op = np.linalg.norm(np.sqrt(wu)[:, None] * C * np.sqrt(wu)[None, :], 2)
        res[name] = (rho, 2 * A.sum(1).max() / 2, 2 * hs, 2 * op)
    return res
for lo, hi in [(-3, 4), (-3.5, 4.5), (-4, 5), (-5, 6)]:
    r = consts(lo, hi)
    print('bulk', lo, hi, {k: tuple(round(float(v), 4) for v in vv) for k, vv in r.items()}, '(specrad 2|P-q|, unweighted sup, 2*HS_L2(du), 2*opnorm)')
