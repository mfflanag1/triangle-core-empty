"""v15 independent numerics: iterate T from delta_0 on an x-lattice, report c_n = 3*(E T(Q_n) + 2 E Q_n).
Law n is kept as a density on the lattice x_j = j h; G_{n+1}(x) = E(x+Y)_+^3/6 is computed from the lattice
convolution f*f by four nested cumulative sums (exact discrete identity for n^3, n^2); laws are re-centred
by an integer lattice shift at each step (Gamma is shift invariant)."""
import sys, numpy as np
from math import gamma

def step(f, h, J):
    # f: density at lattice j=-J..J.  returns G, G' at lattice j=-J..J for the image law (unshifted)
    fY = h * np.convolve(f, f)                # lattice index m=-2J..2J  (density of Y)
    v = fY[::-1]                              # v_j = fY(-j h), j=-2J..2J
    c1 = np.cumsum(v); c2 = np.cumsum(c1); c3 = np.cumsum(c2); c4 = np.cumsum(c3)
    S3 = 6*c4 - 12*c3 + 7*c2 - c1             # sum_j v_j (i-j)_+^3
    S2 = 2*c3 - 3*c2 + c1                     # sum_j v_j (i-j)_+^2
    S3 = np.maximum(S3, 0); S2 = np.maximum(S2, 0)
    # want x_i for i=-3J..3J ; cumsum index = i + 2J  -> pad beyond 4J using polynomial continuation: just take i in [-2J,2J]
    G = h * h**3 / 6 * S3                     # at i = -2J..2J
    Gp = h * h**2 / 2 * S2
    return G, Gp

def run(hexp, L, nmax, save=None):
    h = 2.0 ** -hexp; J = int(round(L / h)); x = np.arange(-J, J + 1) * h
    xp = np.maximum(x, 0)
    f = xp**2 / 2 * np.exp(-xp**3 / 6)       # Q_1 = T(delta_0)
    mean_prev = h * np.sum(x * f)             # E Q_1
    out = [3 * mean_prev]                     # c_0 = 3*(E T(0) + 0)
    # shift law to be centred
    for n in range(1, nmax + 1):
        m = int(round(mean_prev / h))
        fs = np.zeros_like(f)
        if m >= 0: fs[:len(f) - m] = f[m:]
        else: fs[-m:] = f[:len(f) + m]
        lost = h * (f.sum() - fs.sum())
        EQ = h * np.sum(x * fs); mass = h * fs.sum()
        G, Gp = step(fs, h, J)                # at lattice i=-2J..2J
        x2 = np.arange(-2 * J, 2 * J + 1) * h
        fn = Gp * np.exp(-G)
        ET = h * np.sum(x2 * fn); massn = h * fn.sum()
        out.append(3 * (ET + 2 * EQ))
        if n <= 12 or n == nmax:
            print("h=2^-%d n=%2d c_n=%.8f  (mass in %.3e, out %.3e, lost %.1e, E Q_n=%.6f, E T=%.6f)" % (hexp, n, out[-1], mass - 1, massn - 1, lost, EQ, ET), flush=True)
        # next: need law on lattice -J..J: recentre fn (index over -2J..2J)
        mn = int(round(ET / h))
        idx = np.arange(-J, J + 1) + mn + 2 * J
        ok = (idx >= 0) & (idx <= 4 * J)
        f = np.zeros(2 * J + 1); f[ok] = fn[idx[ok]]
        mean_prev = h * np.sum(x * f)         # residual mean (|.|<h)
        if save and n == nmax:
            Gs = np.full(2 * J + 1, np.nan); Gps = np.full(2 * J + 1, np.nan)
            Gs[ok] = G[idx[ok]]; Gps[ok] = Gp[idx[ok]]
            np.savez(save, x=x, G=Gs, Gp=Gps, f=f, c=out[-1])
    return out

if __name__ == '__main__':
    print("c_0 exact = %.8f" % (3 * 6 ** (1 / 3) * gamma(4 / 3)))
    res = {}
    for hexp in (6, 7, 8):
        res[hexp] = run(hexp, 16.0, 40, save=(sys.argv[1] if (len(sys.argv) > 1 and hexp == 8) else None))
    print("n   h=2^-6      h=2^-7      h=2^-8      Richardson(7,8; h^2)")
    for n in range(0, 13):
        a, b, c = res[6][n], res[7][n], res[8][n]
        print("%2d  %.8f  %.8f  %.8f  %.8f   ratio(6-7)/(7-8)=%.2f" % (n, a, b, c, c + (c - b) / 3, (a - b) / (b - c) if b != c else float('nan')))
    print("limit (n=40): %.8f %.8f %.8f" % (res[6][40], res[7][40], res[8][40]))
