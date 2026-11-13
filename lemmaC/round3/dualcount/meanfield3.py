"""Rate of the simplex-relaxed dual first moment (Prop. 1/2 of REPORT.md). Floating point; evidence only.
r(tau) = sup_rho { h(rho) + E log w_tau(s+t+u) - (1/6) E_{rho1 x rho2 x rho3}(s+t'+u'')_+^3 } - log(tau^3/6) - 2,
w_tau(x) = ((x+tau)_+^3 - x_+^3)/6, rho a density on R^3 (mesh eta << tau). Optional tilt theta*sigma to pin the mean."""
import sys, numpy as np

def solve(tau, eta=0.15, lo=-5.0, hi=7.0, iters=3000, damp=0.3, theta=0.0):
    k = np.arange(round(lo / eta), round(hi / eta) + 1); g = k * eta; n = len(g)
    S = g[:, None, None] + g[None, :, None] + g[None, None, :]
    w = (np.maximum(S + tau, 0) ** 3 - np.maximum(S, 0) ** 3) / 6
    logw = np.where(w > 0, np.log(np.maximum(w, 1e-300)), -np.inf) + 3 * np.log(eta)
    F = [0.05 * np.maximum(g + 2, 0) ** 3 for _ in range(3)]
    two = (np.arange(2 * n - 1) + 2 * k[0]) * eta
    cube = np.maximum(g[:, None] + two[None, :], 0) ** 3
    val = None
    for it in range(iters):
        L = logw + theta * S - (F[0][:, None, None] + F[1][None, :, None] + F[2][None, None, :]) / 6
        mx = L.max(); R = np.exp(L - mx); rho = R / R.sum()
        mm = [rho.sum((1, 2)), rho.sum((0, 2)), rho.sum((0, 1))]
        Fn = [cube @ np.convolve(mm[1], mm[2]), cube @ np.convolve(mm[0], mm[2]), cube @ np.convolve(mm[0], mm[1])]
        T3 = float(mm[0] @ Fn[0]); pos = rho > 0
        newval = float(-(rho[pos] * np.log(rho[pos])).sum() + (rho[pos] * logw[pos]).sum() - T3 / 6)
        F = [(1 - damp) * F[i] + damp * Fn[i] for i in range(3)]
        if val is not None and abs(newval - val) < 1e-11 and it > 50: break
        val = newval
    r = newval - np.log(tau ** 3 / 6) - 2
    return dict(tau=tau, eta=eta, theta=theta, sup=round(newval, 4), rate=round(float(r), 4), mu_over_q=round(T3 / 6, 3), mean_sigma=round(float((rho * S).sum()), 3),
                edge=float(max(mm[0][0], mm[0][-1])), it=it)

if __name__ == '__main__':
    eta = float(sys.argv[1])
    for tau in [float(x) for x in sys.argv[2:]]:
        print(solve(tau, eta=eta), flush=True)
