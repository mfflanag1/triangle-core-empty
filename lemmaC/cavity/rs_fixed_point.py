"""Deterministic solution of the replica-symmetric (RS) cavity fixed point for the
local limit of random three-index assignment with pair-additive values.

Local limit: each agent's incident teams have rescaled costs xi = q^(2/3) * team deficit
forming a Poisson process of intensity xi^2/2 on [0, inf). Integer cavity recursion
    X = min_k (xi_k - X'_k - X''_k)
gives  P(X > x) = exp( - E[(x + X' + X'')_+^3] / 6 ).
Energy per team (= per agent of one side), two formulas that must agree:
    a_direct = E[S_+^4] / 8,                 S = X1 + X2 + X3
    a_bethe  = 3 E[X] + E[S_+^4] / 12
Normalisation (each agent in exactly one team): E[S_+^3] / 6 = 1.
"""
import numpy as np

def solve(L=12.0, h=0.002, iters=400, damp=0.5, verbose=False):
    x = np.arange(-L, L + h / 2, h)
    n = len(x)
    # initial guess: X ~ N(1, 1)
    Fbar = 0.5 * (1 - np.tanh((x - 1.0)))
    for it in range(iters):
        p = -np.diff(Fbar, prepend=1.0)          # mass of X in each cell, X ~ x
        p = np.clip(p, 0, None); p[-1] += max(0.0, 1 - p.sum())
        pY = np.convolve(p, p)                   # Y = X'+X'' at grid -2L .. 2L
        yv = 2 * x[0] + h * np.arange(len(pY))
        # g(x) = E[(x+Y)_+^3] for each x on grid, via convolution of pY with t_+^3
        # (x+Y) values: index sum
        t = h * np.arange(-(len(pY) - 1), n) # offsets
        # direct: g(x_i) = sum_j pY[j] * (x_i + yv[j])_+^3
        # compute with cumulative moments: sum over j with yv[j] > -x_i
        c0 = np.cumsum(pY[::-1])[::-1]
        c1 = np.cumsum((pY * yv)[::-1])[::-1]
        c2 = np.cumsum((pY * yv**2)[::-1])[::-1]
        c3 = np.cumsum((pY * yv**3)[::-1])[::-1]
        # first index j with yv[j] > -x_i
        j0 = np.searchsorted(yv, -x, side='right')
        ok = j0 < len(yv)
        jj = np.minimum(j0, len(yv) - 1)
        g = np.where(ok, x**3 * c0[jj] + 3 * x**2 * c1[jj] + 3 * x * c2[jj] + c3[jj], 0.0)
        new = np.exp(-g / 6.0)
        diff = np.max(np.abs(new - Fbar))
        Fbar = damp * Fbar + (1 - damp) * new
        if verbose and it % 50 == 0:
            print(it, diff)
        if diff < 1e-12:
            break
    p = -np.diff(Fbar, prepend=1.0); p = np.clip(p, 0, None); p /= p.sum()
    xm = x - h / 2   # cell midpoints for the mass p (X in (x_{i-1}, x_i])
    pS = np.convolve(np.convolve(p, p), p)
    sv = 3 * xm[0] + h * np.arange(len(pS))
    Sp = np.clip(sv, 0, None)
    EX = (p * xm).sum()
    m3 = (pS * Sp**3).sum(); m4 = (pS * Sp**4).sum()
    return dict(x=xm, p=p, EX=EX, sdX=np.sqrt((p * xm**2).sum() - EX**2),
                norm=m3 / 6, a_direct=m4 / 8, a_bethe=3 * EX + m4 / 12, iters=it, diff=diff)

if __name__ == '__main__':
    for h in [0.01, 0.005, 0.0025]:
        r = solve(h=h)
        print(f"h={h}: iters={r['iters']} diff={r['diff']:.1e} E[X]={r['EX']:.6f} sd={r['sdX']:.6f} "
              f"norm={r['norm']:.6f} a_direct={r['a_direct']:.6f} a_bethe={r['a_bethe']:.6f}")
    np.save(__file__.replace('rs_fixed_point.py', 'out_rs_density.npy'), np.vstack([r['x'], r['p']]))
