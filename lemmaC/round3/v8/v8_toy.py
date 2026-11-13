"""V8 referee: independent toy test of Lemma L (box Legendre identity), Lemma Q (adversarial search for a violation),
and of Prop 4.5 used at mu' (3H - 2S0 + <ell,.> <= Psi0(u0)).  Evidence only.  Own implementation."""
import itertools, math, sys
import numpy as np
from scipy.optimize import minimize, differential_evolution
from scipy.stats import poisson
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
def phi(x): return -x * math.log(x) if x > 0 else 0.0
def build(K, B, th, dmax):
    eta = K / B; b = np.arange(1, B + 1)
    lam = ((eta * b) ** 3 - (eta * (b - 1)) ** 3) / 6
    Gf = lambda x: -np.exp(-th * x) * (x * x / (2 * th) + x / th ** 2 + 1 / th ** 3)
    kap = Gf(eta * b) - Gf(eta * (b - 1)); J = float(-Gf(K))
    Dm = np.array(list(itertools.product(range(dmax + 1), repeat=B)), dtype=float)
    pr = np.prod(poisson.pmf(Dm, lam[None, :]), axis=1)
    assert 1 - pr.sum() < 1e-12, 1 - pr.sum()
    return lam, kap, J, Dm, pr
class M:
    def __init__(s, K, B, th, dmax):
        s.K, s.B, s.th = K, B, th
        s.lam, s.kap, s.J, s.Dm, s.pr = build(K, B, th, dmax)
        s.ell = np.concatenate([np.log(s.kap / s.lam), [math.log(s.J)]]); s.ellmax = np.abs(s.ell).max()
    def Phi(s, u):
        w = np.exp(u); W = w[-1] + s.Dm @ w[:-1]
        g = np.concatenate([(s.pr / W) @ (s.Dm * w[None, :-1]), [s.pr @ (w[-1] / W)]])
        return s.pr @ np.log(W), g
    def Psi0(s, u):
        w = np.exp(u); W = w[-1] + s.Dm @ w[:-1]
        rho = np.array([s.pr @ (w[b] / (W + w[b])) for b in range(s.B)]); eps = s.pr @ (w[-1] / W)
        C = np.exp(-3 * u[:-1]) * s.kap / s.lam; Cinf = math.exp(-3 * u[-1]) * s.J
        return 3 * (s.pr @ np.log(W)) + s.lam @ np.log((1 - rho) ** 3 + C * rho ** 3) - 3 * eps + Cinf * eps ** 3
    def HA(s, mu, A):
        f = lambda u: (lambda v, g: (v - mu @ u, g - mu))(*s.Phi(u))
        best = None
        for st in (np.zeros(s.B + 1), rng.uniform(-A / 2, A / 2, s.B + 1), -A / 2 * np.sign(mu - 1.0 / (s.B + 1))):
            r = minimize(f, st, jac=True, method="L-BFGS-B", bounds=[(-A / 2, A / 2)] * (s.B + 1), options={"ftol": 1e-16, "gtol": 1e-12, "maxiter": 10000})
            if best is None or r.fun < best.fun: best = r
        return best.fun, best.x
def Stil(mu, nu):
    return sum(phi(m) + phi(max(n - m, 0.0)) - phi(n) for m, n in zip(mu[:-1], nu)) + phi(mu[-1]) + mu[-1]
def S0(mu, lam): return Stil(mu, lam)
tot_viol = -1e9; worst_id = 0.0
for (K, B, th, dmax) in ((1.5, 1, 2.0, 60), (2.2, 1, 1.0, 60), (2.0, 2, 1.0, 34), (2.6, 2, 1.5, 40)):
    m = M(K, B, th, dmax)
    # min of Psi0 over u0 (tightest right-hand side)
    best = min((minimize(m.Psi0, rng.normal(0, 1, B + 1), method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-13, "maxiter": 40000}) for _ in range(4)), key=lambda r: r.fun)
    psimin = best.fun
    # (a) consistency of Prop 4.5 as used: sup over mu' in D of 3H - 2S0 + <ell,mu'>, param by u (mu'=grad Phi(u))
    def Lfun(u):
        v, g = m.Phi(u); H = v - g @ u
        return -(3 * H - 2 * S0(g, m.lam) + m.ell @ g)
    bestL = min((minimize(Lfun, rng.normal(0, 1, B + 1), method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-13, "maxiter": 40000}) for _ in range(4)), key=lambda r: r.fun)
    print("K=%g B=%d theta=%g: lam=%s  min_u Psi0 = %.8f ; sup_{mu' in D} [3H-2S0+<ell,mu'>] = %.8f (must be <= , tight if equal)" % (K, B, th, np.round(m.lam, 4), psimin, -bestL.fun))
    for A in (6.0, 10.0, 16.0):
        c1 = (3 * A - 2 * m.ellmax - 4) / 8
        if c1 < 1: continue
        def gap(x, ret_id=False):
            # x in [0,1]^(2B): fr (mu_b = fr*nu_b), nu = lam*(1+pert), pert in [-.6,.6]
            fr = x[:B]; nu = m.lam * (1 + 1.2 * (x[B:] - 0.5))
            if c1 < 1 + max(0.0, math.log(nu.max())): return 1e3 if not ret_id else (-1e3, 0.0)  # hypothesis of Lemma E fails: point skipped
            mub = fr * nu
            if mub.sum() > 1: mub = mub / mub.sum()          # eps = 0 boundary
            mu = np.concatenate([mub, [max(1 - mub.sum(), 0.0)]])
            HA, us = m.HA(mu, A)
            lhs = 3 * HA - 2 * Stil(mu, nu) + m.ell @ mu
            rhs = psimin + 4 * c1 * np.abs(nu - m.lam).sum() + (6 * B + 2) * math.exp(-1 - c1)
            if ret_id:
                v, mup = m.Phi(us); return lhs - rhs, HA - ((v - mup @ us) - A / 2 * np.abs(mu - mup).sum())
            return -(lhs - rhs)
        # random + boundary samples
        wg = -1e9
        for trial in range(120):
            x = rng.uniform(0, 1, 2 * B)
            if trial % 4 == 0: x[:B] = rng.choice([0.0, 1.0, 1e-9, 1 - 1e-9], B)
            if trial % 5 == 0: x[B:] = rng.choice([0.0, 1.0, 0.5], B)
            g_, id_ = gap(x, True); wg = max(wg, g_); worst_id = max(worst_id, abs(id_))
        r = differential_evolution(gap, [(0, 1)] * (2 * B), seed=int(rng.integers(1 << 30)), maxiter=25, popsize=10, tol=1e-9, polish=False)
        wg = max(wg, -r.fun); tot_viol = max(tot_viol, wg)
        print("   A=%g c'=%.3f: max(lhs-rhs) over random+DE = %.5f at x=%s ; max |identity residual| so far %.2e" % (A, c1, wg, np.round(r.x, 3), worst_id))
print("OVERALL max(lhs - rhs) = %.5f (Lemma Q violated iff > 0); max identity residual %.2e" % (tot_viol, worst_id))
