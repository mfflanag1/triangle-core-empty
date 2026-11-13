"""Floating-point sanity check of Lemma Q of THEOREM_I_PA.md in toy models (evidence only).
Model: B bins on [0,K], Poisson reference star, high option.  For random (mu, nu):
 (1) identity  H_A(mu) = H(mu') - (A/2)|mu - mu'|_1  with mu' = grad Phi(u*), u* the box minimiser;
 (2) bound  3 H_A(mu) - 2 S~(mu,nu) + <ell,mu> <= Psi_0(u_0) + 4 c' d1 + (6B+2) exp(-1-c'),
     c' = (3A - 2 ell_max - 4)/8, d1 = |nu - lam|_1, for a random u_0;
 (3) Lemma E: phi(a) - phi(b) <= c'|a-b| + exp(-1-c') for a,b in [0,M], c' >= 1 + ln^+ M.
"""
import itertools, math, sys
import numpy as np
from scipy.optimize import minimize
from scipy.stats import poisson
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
def model(K, B, th, dmax=28):
    eta = K / B
    lam = np.array([eta ** 3 * (3 * b * b - 3 * b + 1) / 6 for b in range(1, B + 1)])
    G = lambda x: -math.exp(-th * x) * (x * x / (2 * th) + x / th ** 2 + 1 / th ** 3)
    kap = np.array([G(eta * b) - G(eta * (b - 1)) for b in range(1, B + 1)])
    J = -G(K)
    grids = list(itertools.product(range(dmax + 1), repeat=B))
    Dm = np.array(grids, dtype=float)
    pr = np.prod(poisson.pmf(Dm, lam[None, :]), axis=1)
    return lam, kap, J, Dm, pr
def Phi_grad(u, Dm, pr):
    w = np.exp(u); W = w[-1] + Dm @ w[:-1]
    val = pr @ np.log(W)
    g = np.concatenate([(pr / W) @ (Dm * w[None, :-1]), [pr @ (w[-1] / W)]])
    return val, g
def Psi0(u, lam, kap, J, th, Dm, pr):
    w = np.exp(u); W = w[-1] + Dm @ w[:-1]
    Phi = pr @ np.log(W)
    rho = np.array([pr @ (w[b] / (W + w[b])) for b in range(len(lam))]); eps = pr @ (w[-1] / W)
    C = np.exp(-3 * u[:-1]) * kap / lam; Cinf = math.exp(-3 * u[-1]) * J
    return 3 * Phi + lam @ np.log((1 - rho) ** 3 + C * rho ** 3) - 3 * eps + Cinf * eps ** 3
def phi(x): return -x * math.log(x) if x > 0 else 0.0
def Stil(mu, nu):
    s = sum(phi(m) + phi(n - m) - phi(n) for m, n in zip(mu[:-1], nu))
    return s + phi(mu[-1]) + mu[-1]
worst_id = 0; worst_gap = -1e9; ntest = 0
for (K, B, th, A) in ((2.0, 2, 1.0, 8.0), (2.0, 2, 1.0, 14.0), (2.4, 2, 1.5, 12.0), (1.5, 1, 2.0, 10.0), (2.5, 3, 1.2, 14.0)):
    lam, kap, J, Dm, pr = model(K, B, th, dmax=28 if B < 3 else 16)
    ell = np.concatenate([np.log(kap / lam), [math.log(J)]]); ellmax = max(th * K, abs(math.log(J)))
    c1 = (3 * A - 2 * ellmax - 4) / 8
    psis = [Psi0(rng.normal(0, 1.5, B + 1), lam, kap, J, th, Dm, pr) for _ in range(6)]
    r = minimize(lambda u: Psi0(u, lam, kap, J, th, Dm, pr), np.zeros(B + 1), method="Nelder-Mead", options={"xatol": 1e-9, "fatol": 1e-12, "maxiter": 20000})
    psi_min = min(min(psis), r.fun)
    for trial in range(40):
        d = rng.uniform(-1, 1, B) * rng.choice([0.0, 0.01, 0.1, 0.3]) * lam
        nu = lam + d; d1 = np.abs(d).sum()
        if c1 < 1 + max(0.0, math.log(max(1.0, nu.max()))): continue
        # random profile: 0<=mu_b<=nu_b, sum<=1; include boundary cases
        fr = rng.uniform(0, 1, B) ** rng.choice([0.2, 1.0, 5.0])
        mub = fr * nu
        if rng.random() < 0.2: mub[rng.integers(B)] = nu[rng.integers(B)] if B == 1 else mub[rng.integers(B)]
        if rng.random() < 0.15: mub[rng.integers(B)] = 0.0
        if mub.sum() > 1: mub *= rng.uniform(0.3, 1.0) / mub.sum()
        mu = np.concatenate([mub, [1 - mub.sum()]])
        f = lambda u: (lambda v, g: (v - mu @ u, g - mu))(*Phi_grad(u, Dm, pr))
        best = None
        for start in (np.zeros(B + 1), rng.uniform(-A / 2, A / 2, B + 1)):
            res = minimize(f, start, jac=True, method="L-BFGS-B", bounds=[(-A / 2, A / 2)] * (B + 1), options={"ftol": 1e-15, "gtol": 1e-11, "maxiter": 5000})
            if best is None or res.fun < best.fun: best = res
        HA = best.fun; us = best.x
        v, mup = Phi_grad(us, Dm, pr); Hmup = v - mup @ us
        ident = HA - (Hmup - A / 2 * np.abs(mu - mup).sum())
        worst_id = max(worst_id, abs(ident))
        lhs = 3 * HA - 2 * Stil(mu, nu) + ell @ mu
        rhs = psi_min + 4 * c1 * d1 + (6 * B + 2) * math.exp(-1 - c1)
        worst_gap = max(worst_gap, lhs - rhs); ntest += 1
    print("K=%g B=%d theta=%g A=%g c'=%.3f: min Psi0 found %.5f; tests so far %d; max|identity residual| %.2e; max(lhs-rhs) %.4f" % (K, B, th, A, c1, psi_min, ntest, worst_id, worst_gap))
wE = -1e9
for _ in range(200000):
    M = rng.choice([0.5, 1.0, 1.3, 3.0]); a, b = rng.uniform(0, M, 2) * rng.choice([1, 1e-3, 1e-6], 2)
    c = 1 + max(0, math.log(M)) + rng.exponential(3)
    wE = max(wE, phi(a) - phi(b) - c * abs(a - b) - math.exp(-1 - c))
print("Lemma E: max of phi(a)-phi(b)-c|a-b|-exp(-1-c) over 200000 samples: %.3e (must be <= 0)" % wE)
