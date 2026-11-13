"""Annealed rate of the grid first moment. Floating point; evidence only.
Phi_eps = sup_rho { H(rho) + E_rho log w_eps(s+t+u) - (1/6) E_{rho1 x rho2 x rho3} (s+t'+u'')_+^3 },  rho on (eps Z)^3,
w_eps(x) = ((x+3eps)_+^3 - x_+^3)/6.   (q!)^2 E N_eps = exp(q (Phi_eps - 2) + o(q)).
Fixed point: rho ~ w(s+t+u) exp(-(F1(s)+F2(t)+F3(u))/6); value = log Z + T3/3."""
import sys, numpy as np

def solve(eps, lo=-5.0, hi=7.0, iters=4000, damp=0.3, seed=None, cont=False):
    g = np.arange(round(lo / eps), round(hi / eps) + 1) * eps
    n = len(g)
    S = g[:, None, None] + g[None, :, None] + g[None, None, :]
    if cont:   # continuous-volume normalisation: weight x^2/2, counting measure eps^3 per point
        logw = np.where(S > 0, np.log(np.maximum(S, 1e-300) ** 2 / 2), -np.inf) + 3 * np.log(eps)
    else:
        w = (np.maximum(S + 3 * eps, 0) ** 3 - np.maximum(S, 0) ** 3) / 6
        logw = np.where(w > 0, np.log(np.maximum(w, 1e-300)), -np.inf)
    rng = np.random.default_rng(seed)
    F = [0.05 * np.maximum(g + 2, 0) ** 3 + (0 if seed is None else rng.normal(0, 0.5, n)) for _ in range(3)]
    # sums index: value g_i + g_j = eps*(k_i + k_j); cube table for s + (t+u)
    k = np.round(g / eps).astype(int)
    two = (np.arange(2 * n - 1) + 2 * k[0]) * eps                 # values of t'+u''
    cube = np.maximum(g[:, None] + two[None, :], 0) ** 3          # (s, t'+u'')
    val = None
    for it in range(iters):
        L = logw - (F[0][:, None, None] + F[1][None, :, None] + F[2][None, None, :]) / 6
        mx = L.max(); R = np.exp(L - mx); Z = R.sum(); rho = R / Z; logZ = np.log(Z) + mx
        m = [rho.sum((1, 2)), rho.sum((0, 2)), rho.sum((0, 1))]
        Fn = [cube @ np.convolve(m[1], m[2]), cube @ np.convolve(m[0], m[2]), cube @ np.convolve(m[0], m[1])]
        T3 = float(m[0] @ Fn[0])
        newval = logZ + T3 / 3
        # recentre gauge: nothing needed (shift-equivariant); damped update
        F = [(1 - damp) * F[i] + damp * Fn[i] for i in range(3)]
        if val is not None and abs(newval - val) < 1e-12 and it > 50: break
        val = newval
    # direct evaluation of the functional at rho (not using the fixed-point identity)
    H = -(rho[rho > 0] * np.log(rho[rho > 0])).sum()
    direct = H + (rho[rho > 0] * logw[rho > 0]).sum() - T3 / 6
    msig = float((rho * S).sum()); edge = float(max(m[0][0], m[0][-1], m[1][0], m[1][-1], m[2][0], m[2][-1]))
    return dict(eps=eps, Phi=float(direct), fp=float(newval), it=it, T3=T3, mean_sigma=msig, mean_s=float(m[0] @ g), edge=edge,
                sd_s=float(np.sqrt(m[0] @ g ** 2 - (m[0] @ g) ** 2)))

if __name__ == '__main__':
    for eps in [float(x) for x in sys.argv[1:]]:
        r = solve(eps)
        print('grid ', {k_: (round(v, 5) if isinstance(v, float) else v) for k_, v in r.items()}, 'rate Phi-2 =', round(r['Phi'] - 2, 4), flush=True)
