"""V9 referee: independent discretisation of T in QUANTILE space (no x-grid, no FFT).
Nodes: composite Gauss-Legendre in w = ln(theta); law of Q(E) = sum_i p_i delta_{q_i}.
g(t) = sum_ij p_i p_j (t+q_i+q_j)_+^3/6 evaluated exactly for the discrete law by sorted suffix sums.
Floating point: evidence only."""
import numpy as np
from numpy.polynomial.legendre import leggauss

class Grid:
    def __init__(self, wmin=-36.0, wmax=4.0, panel=0.25, k=8):
        npan = int(round((wmax - wmin) / panel))
        self.edges = wmin + panel * np.arange(npan + 1)
        gx, gw = leggauss(k)
        w = (self.edges[:-1, None] + panel * (gx[None, :] + 1) / 2).ravel()
        wt = np.tile(panel / 2 * gw, npan)
        self.w = w; self.th = np.exp(w)
        self.p = wt * np.exp(w - self.th)
        self.p /= self.p.sum()

class Pair:
    """sorted pair sums of a discrete law (q,p) with suffix moment sums."""
    def __init__(self, q, p):
        Y = (q[:, None] + q[None, :]).ravel(); W = (p[:, None] * p[None, :]).ravel()
        i = np.argsort(Y, kind='stable'); self.Y = Y[i]; W = W[i]
        suf = lambda a: np.concatenate([np.cumsum(a[::-1])[::-1], [0.0]])
        Yc = self.Y
        self.S = [suf(W), suf(W * Yc), suf(W * Yc ** 2), suf(W * Yc ** 3), suf(W * Yc ** 4)]
        self.A = [suf(W), suf(W * abs(Yc)), suf(W * Yc ** 2), suf(W * abs(Yc) ** 3)]
        self.tot = W.sum()
    def M(self, k, t):
        """E (t+Y)_+^k / k!  for k = 0..4"""
        t = np.asarray(t, float); j = np.searchsorted(self.Y, -t, side='right'); S = self.S
        if k == 0: return S[0][j]
        if k == 1: return t * S[0][j] + S[1][j]
        if k == 2: return (t * t * S[0][j] + 2 * t * S[1][j] + S[2][j]) / 2
        if k == 3: return (t ** 3 * S[0][j] + 3 * t * t * S[1][j] + 3 * t * S[2][j] + S[3][j]) / 6
        if k == 4: return (t ** 4 * S[0][j] + 4 * t ** 3 * S[1][j] + 6 * t * t * S[2][j] + 4 * t * S[3][j] + S[4][j]) / 24
    def cond(self, t):
        """cancellation amplification factor of the suffix-sum evaluation of M_3 at t"""
        t = np.asarray(t, float); j = np.searchsorted(self.Y, -t, side='right'); A = self.A; a = abs(t)
        return (a ** 3 * A[0][j] + 3 * a * a * A[1][j] + 3 * a * A[2][j] + A[3][j]) / 6 / np.maximum(self.M(3, t), 1e-300)
    def ginv(self, theta):
        theta = np.asarray(theta, float)
        lo = np.full(theta.shape, -self.Y[-1]); hi = (6 * theta / self.tot) ** (1 / 3) - self.Y[0] + 1.0
        for _ in range(90):
            mid = 0.5 * (lo + hi); up = self.M(3, mid) < theta
            lo = np.where(up, mid, lo); hi = np.where(up, hi, mid)
        return 0.5 * (lo + hi)

def T(q, G, levels=None):
    """returns T(Q) at the nodes (or at `levels`) and the slope-at-level s there."""
    m = float((q * G.p).sum())                      # centre first (T(Q+a) = T(Q) - 2a): avoids cancellation
    P = Pair(q - m, G.p); th = G.th if levels is None else levels
    t = P.ginv(th); P.last_cond = P.cond(t); P.shift = m
    return t - 2 * m, P.M(2, t), P

def mean(q, G): return float((q * G.p).sum())

def Gamma(q, G):
    t, _, _ = T(q, G); return mean(t, G) + 2 * mean(q, G)
