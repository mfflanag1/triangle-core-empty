"""V10: exact (unrounded) iteration c_n = 3 Gamma(T^n delta_0), own code, plain float64.
Grid x_i = -L + i dx.  Law of X_n by G_n (survival exp(-G_n)) and G_n' on the grid; density f = G' exp(-G).
f_Y = f*f (trapezoid = plain discrete convolution, f vanishes at both ends),
G_{n+1}(x) = dx sum_j f_Y(y_j) (x+y_j)_+^3/6,  G_{n+1}'(x) = dx sum_j f_Y(y_j) (x+y_j)_+^2/2  (kink on a grid node: O(dx^4)).
Each iterate is re-centred by a whole number of cells (shift law T(Q+a) = T(Q) - 2a; Gamma is shift invariant):
c_n = 3 (E T(X_n) + 2 E X_n).
Usage: v10_unrounded.py NMAX INVDX [L]"""
import sys
import numpy as np
NMAX = int(sys.argv[1]); INV = int(sys.argv[2]); L = float(sys.argv[3]) if len(sys.argv) > 3 else 16.0
dx = 1.0 / INV
x = np.arange(-int(L * INV), int(L * INV) + 1) * dx
n = len(x)
y = np.arange(-2 * int(L * INV), 2 * int(L * INV) + 1) * dx


def mean(f):
    return float((x * f).sum() * dx)


def T(f):
    fY = np.convolve(f, f) * dx
    G = np.empty(n); Gp = np.empty(n)
    for s in range(0, n, 256):
        U = np.maximum(x[s:s + 256, None] + y[None, :], 0.0); U2 = U * U
        G[s:s + 256] = (U2 * U) @ fY * (dx / 6); Gp[s:s + 256] = U2 @ fY * (dx / 2)
    return Gp * np.exp(-G)


# X_0 = delta_0: T gives G = x_+^3/6
xp = np.maximum(x, 0.0)
f = xp ** 2 / 2 * np.exp(-xp ** 3 / 6)
print("dx = 1/%d, L = %g" % (INV, L))
print(" c_0 = %.9f   (mass of X_1: %.12f)" % (3 * mean(f), f.sum() * dx), flush=True)
for k in range(1, NMAX + 1):
    mu = mean(f); sh = int(round(mu / dx))
    fc = np.zeros(n)                      # law of X_k - sh*dx
    if sh >= 0:
        fc[:n - sh] = f[sh:]
    else:
        fc[-sh:] = f[:n + sh]
    muc = mean(fc)
    fn = T(fc)
    print(" c_%d = %.9f   (mass %.12f, edge densities %.1e %.1e)" % (k, 3 * (mean(fn) + 2 * muc), fn.sum() * dx, fn[0], fn[-1]), flush=True)
    f = fn
