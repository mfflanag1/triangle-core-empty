"""Annealed rates, variants. Floating point; evidence only.
mode 'cont' : phi = rate of q^{2q} * E[rescaled (2q)-volume of the certificate set of the identity]
mode 'grid2': (s,t) on (eps Z)^2, third price x continuous (fine mesh eta = eps/m), diagonal weight (s+t+x+2eps)^2/2
              [certificate rounded down in two sectors, third sector eliminated]
Value = H + E log w - T3/6 evaluated directly at the fixed point of rho ~ w exp(-(F1+F2+F3)/6)."""
import sys, numpy as np

def solve(mode, eps, m=4, lo=-5.0, hi=7.0, iters=3000, damp=0.3):
    eta = eps if mode == 'cont' else eps / m
    k = np.arange(round(lo / eta), round(hi / eta) + 1); g = k * eta; n = len(g)
    S = g[:, None, None] + g[None, :, None] + g[None, None, :]
    if mode == 'cont':
        logw = np.where(S > 0, np.log(np.maximum(S, 1e-300) ** 2 / 2), -np.inf) + 3 * np.log(eta)
    else:
        coarse = (k % m == 0)
        logw = np.where(S + 2 * eps > 0, np.log(np.maximum(S + 2 * eps, 1e-300) ** 2 / 2), -np.inf) + np.log(eta)
        logw = np.where(coarse[:, None, None] & coarse[None, :, None], logw, -np.inf)
    F = [0.05 * np.maximum(g + 2, 0) ** 3 for _ in range(3)]
    two = (np.arange(2 * n - 1) + 2 * k[0]) * eta
    cube = np.maximum(g[:, None] + two[None, :], 0) ** 3
    val = None
    for it in range(iters):
        L = logw - (F[0][:, None, None] + F[1][None, :, None] + F[2][None, None, :]) / 6
        mx = L.max(); R = np.exp(L - mx); rho = R / R.sum()
        mm = [rho.sum((1, 2)), rho.sum((0, 2)), rho.sum((0, 1))]
        Fn = [cube @ np.convolve(mm[1], mm[2]), cube @ np.convolve(mm[0], mm[2]), cube @ np.convolve(mm[0], mm[1])]
        T3 = float(mm[0] @ Fn[0])
        pos = rho > 0
        newval = float(-(rho[pos] * np.log(rho[pos])).sum() + (rho[pos] * logw[pos]).sum() - T3 / 6)
        F = [(1 - damp) * F[i] + damp * Fn[i] for i in range(3)]
        if val is not None and abs(newval - val) < 1e-11 and it > 50: break
        val = newval
    return dict(mode=mode, eps=eps, value=round(newval, 4), rate=round(newval - 2, 4), T3=round(T3, 3), mean_sum=round(float((rho * S).sum()), 3),
                edge=float(max(mm[0][0], mm[0][-1], mm[2][0], mm[2][-1])), it=it)

if __name__ == '__main__':
    mode = sys.argv[1]
    for eps in [float(x) for x in sys.argv[2:]]:
        print(solve(mode, eps), flush=True)
