"""Endogeny test for the replica-symmetric recursive tree process X = min_k(xi_k - X_j - X_l) on the K-truncated
Poisson hypertree (team costs at an agent: intensity xi^2/2 on [0,K]).  Fixed point G_K by damped grid iteration;
then sample depth-g trees and feed B independent iid-G_K boundary conditions to the SAME tree.
Reports corr of root message between replicas and P(same root team).  Endogenous => both -> 1 as g grows."""
import numpy as np, sys, time
from scipy.signal import fftconvolve
K = float(sys.argv[1]); gmax = int(sys.argv[2]); LBIG = 30.0
def fixed_point(K, L=12.0, h=0.01, iters=300):
    x = np.arange(-L, L + h / 2, h); n = len(x)
    f = np.exp(-(x - 0.6) ** 2 / 0.5); f /= f.sum() * h
    s = np.arange(2 * n - 1) * h - 2 * L
    for it in range(iters):
        g = fftconvolve(f, f) * h
        Phi = np.empty(n)
        for a in range(0, n, 256):
            xx = np.clip(np.minimum(x[a:a + 256, None] + s[None, :], K), 0, None)
            Phi[a:a + 256] = (xx ** 3 * g[None, :]).sum(1) * h / 6
        G = np.exp(-Phi); fn = np.clip(-np.gradient(G, h), 0, None); fn /= fn.sum() * h
        diff = np.abs(fn - f).sum() * h; f = 0.5 * fn + 0.5 * f
        if diff < 1e-9: break
    # energy: E[xi; root team] = int_0^K xi * xi^2/2 * P(xi < X1+X2+X3) dxi
    g2 = fftconvolve(f, f) * h; g3 = fftconvolve(g2, f) * h; s3 = np.arange(len(g3)) * h - 3 * L
    tail3 = 1 - np.cumsum(g3) * h
    t = np.arange(0, K, h); rho = np.interp(t, s3, tail3)
    return x, f, diff, (t ** 3 / 2 * rho).sum() * h, (t ** 2 / 2 * rho).sum() * h
x, f, diff, en, norm = fixed_point(K)
cdf = np.cumsum(f) * (x[1] - x[0]); cdf /= cdf[-1]
print('K=%g fixed point: diff=%.1e energy=%.4f norm=%.4f meanX=%.4f sdX=%.4f' % (K, diff, en, norm, (x * f).sum() * (x[1] - x[0]), np.sqrt(((x - (x * f).sum() * (x[1] - x[0])) ** 2 * f).sum() * (x[1] - x[0]))), flush=True)
rng = np.random.default_rng(5)
def draw(shape): return np.interp(rng.random(shape), cdf, x)
lam = K ** 3 / 6
def one_tree(g, B):
    levels = []; n = 1
    for l in range(g):
        deg = rng.poisson(lam, n); nt = int(deg.sum())
        par = np.repeat(np.arange(n), deg); xi = K * rng.random(nt) ** (1 / 3)
        levels.append((n, deg, par, xi)); n = 2 * nt
    X = draw((n, B))
    for l in range(g - 1, -1, -1):
        n, deg, par, xi = levels[l]
        val = xi[:, None] - X[0::2] - X[1::2]
        if l == 0:
            if len(xi) == 0: return None
            return val.min(0), xi[val.argmin(0)], val.argmin(0)
        Xn = np.full((n, B), LBIG)
        nz = deg > 0; starts = (np.cumsum(deg) - deg)[nz]
        Xn[nz] = np.minimum.reduceat(val, starts, axis=0)
        X = Xn
for g in range(1, gmax + 1):
    ntree = {1: 2000, 2: 2000, 3: 800, 4: 200, 5: 24, 6: 4}[g]; B = 8 if g <= 4 else 4
    t0 = time.time(); xs = []; same = []; cost = []
    for _ in range(ntree):
        r = one_tree(g, B)
        if r is None: continue
        xo, xi_o, am = r; xs.append(xo); cost.append(xi_o.mean())
        same.append(np.mean([am[i] == am[j] for i in range(B) for j in range(i + 1, B)]))
    xs = np.array(xs); within = xs.var(1, ddof=1).mean(); total = xs.var()
    print('g=%d trees=%d B=%d  E[X_o]=%.3f  var_total=%.3f  var_within_tree=%.3f  corr=%.3f  P(same root team)=%.3f +- %.3f  E[xi root team]=%.3f  (%.0fs)'
          % (g, len(xs), B, xs.mean(), total, within, 1 - within / total, np.mean(same), np.std(same) / np.sqrt(len(same)), np.mean(cost), time.time() - t0), flush=True)
