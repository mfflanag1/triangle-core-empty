"""V12 adversarial test of Lemma 4.2 (M <= Rt) and of the tangent inequality (4.1) on small models
(B = 1, 2, 3 bins, arbitrary lambda_b, high option included).  Poisson expectations by truncated exact sums.
Searches (differential evolution + polish) maximise
  E1: lambda_max( Rt^-1/2 M Rt^-1/2 )                      (claim: <= 1)
  E2: [LHS - RHS of (4.1)] / |mu1 - mu0|^2, mu1 = mu(u1)   (claim: <= 0)
  E3: LHS - RHS with H(mu) from a direct minimisation over u at mu on a grid of P (not parametrised by u1)."""
import sys, time, itertools
import numpy as np
from scipy.optimize import differential_evolution, minimize
from scipy.special import gammaln
rng = np.random.default_rng(12)

class Model:
    def __init__(self, lam):
        self.lam = np.asarray(lam, float); self.B = len(lam)
        axes = []; pm = []
        for l in self.lam:
            n = int(l + 10 * np.sqrt(l) + 22); k = np.arange(n + 1)
            axes.append(k); pm.append(np.exp(-l + k * np.log(l) - gammaln(k + 1)))
        G = np.meshgrid(*axes, indexing="ij"); self.D = np.stack([g.ravel() for g in G]).astype(float)
        P = pm[0]
        for q in pm[1:]: P = np.multiply.outer(P, q)
        self.P = P.ravel()
    def stats(self, u, hess=False):
        B = self.B; w = np.exp(u[:B]); winf = np.exp(u[B])
        W = winf + w @ self.D
        Phi = self.P @ np.log(W)
        rho = np.array([self.P @ (w[i] / (W + w[i])) for i in range(B)])
        omr = np.array([self.P @ (W / (W + w[i])) for i in range(B)])          # 1 - rho, computed directly
        eps = self.P @ (winf / W)
        out = dict(Phi=Phi, rho=rho, omr=omr, eps=eps, mu=self.lam * rho)
        if hess:
            pi = np.vstack([self.D * w[:, None] / W, winf / W])
            mu = pi @ self.P
            M = np.diag(mu) - (pi * self.P) @ pi.T
            out["M"] = M; out["mu_direct"] = mu
        return out

def eig_ratio(mdl, u):
    s = mdl.stats(u, hess=True)
    Rt = np.concatenate([mdl.lam * s["rho"] * s["omr"], [s["eps"]]])
    if np.any(Rt <= 1e-280): return 0.0
    A = s["M"] / np.sqrt(np.outer(Rt, Rt))
    return float(np.linalg.eigvalsh((A + A.T) / 2)[-1])

def h(r, omr): return -(r * np.log(r) + omr * np.log(omr))

def gap(mdl, u0, u1):
    B = mdl.B; s0 = mdl.stats(u0); s1 = mdl.stats(u1)
    mu1, e1 = s1["mu"], s1["eps"]
    S1 = (mdl.lam * h(s1["rho"], s1["omr"])).sum() + e1 - e1 * np.log(e1)
    lhs = s1["Phi"] - mu1 @ u1[:B] - e1 * u1[B] - S1
    rhs = s0["Phi"] + (mdl.lam * np.log(s0["omr"])).sum() - s0["eps"] - mu1 @ (u0[:B] + np.log(s0["omr"] / s0["rho"])) - e1 * (u0[B] - np.log(s0["eps"]))
    dist2 = ((mu1 - s0["mu"]) ** 2).sum() + (e1 - s0["eps"]) ** 2
    return lhs - rhs, dist2

TB = float(sys.argv[1]) if len(sys.argv) > 1 else 60.0     # seconds per search
def run(name, f, bounds, seed):
    t = time.time(); best = [-np.inf, None]
    def obj(x):
        v = f(x)
        if not np.isfinite(v): return 1e9
        if v > best[0]: best[0], best[1] = v, x.copy()
        return -v
    def cb(xk, convergence=None): return time.time() - t > TB
    differential_evolution(obj, bounds, seed=seed, popsize=14, maxiter=100000, tol=0, mutation=(0.4, 1.2), recombination=0.8, polish=False, callback=cb)
    r = minimize(obj, best[1], method="Nelder-Mead", options=dict(maxfev=3000, xatol=1e-10, fatol=1e-14))
    print("%-34s max = %+.12e   at %s   [%.0fs]" % (name, best[0], np.array2string(best[1], precision=3), time.time() - t), flush=True)
    return best

LL, LH, UM = np.log(1e-3), np.log(6.0), 9.0
for B in (1, 2, 3):
    lh = LH if B < 3 else np.log(3.0)
    bl = [(LL, lh)] * B; bu = [(-UM, UM)] * (B + 1)
    run("E1 eig ratio  B=%d" % B, lambda x: eig_ratio(Model(np.exp(x[:B])), x[B:]), bl + bu, 100 + B)
    def e2(x):
        mdl = Model(np.exp(x[:B])); g, d2 = gap(mdl, x[B:2 * B + 1], x[2 * B + 1:])
        return g / d2 if d2 > 1e-8 else -1e3
    run("E2 gap/|dmu|^2 B=%d" % B, e2, bl + bu + bu, 200 + B)
    def e2b(x):
        mdl = Model(np.exp(x[:B])); g, d2 = gap(mdl, x[B:2 * B + 1], x[2 * B + 1:])
        return g if d2 > 1e-4 else -1e3
    run("E2b raw gap (|dmu|>0.01) B=%d" % B, e2b, bl + bu + bu, 300 + B)

# E3: H(mu) by direct minimisation, mu on a grid of P; many tangent points u0
print("E3: H by direct minimisation over u (shift fixed by u_inf = 0 ... all coordinates free, BFGS)")
worst = -np.inf; nfin = 0; ntot = 0
for B, lam in ((1, [0.7]), (1, [3.0]), (2, [0.5, 1.5]), (2, [2.0, 0.2])):
    mdl = Model(lam); lam = np.array(lam)
    u0s = [rng.uniform(-6, 6, B + 1) for _ in range(40)]
    s0s = [mdl.stats(v) for v in u0s]
    grid = np.linspace(0.0, 1.0, 13)
    for fr in itertools.product(grid, repeat=B):
        mu = lam * np.array(fr); e = 1 - mu.sum()
        if e < 0: continue
        ntot += 1
        f = lambda v: mdl.stats(v)["Phi"] - mu @ v[:B] - e * v[B]
        r = minimize(f, np.zeros(B + 1), method="BFGS", options=dict(gtol=1e-10, maxiter=400))
        Hval = r.fun                                   # >= true H (an overestimate of the LHS)
        if Hval < -30: continue                        # H = -inf (or hugely negative): nothing to check
        nfin += 1
        S0 = sum(l * (0 if x in (0.0, 1.0) else -(x * np.log(x) + (1 - x) * np.log(1 - x))) for l, x in zip(lam, fr)) + (e - e * np.log(e) if e > 0 else 0.0)
        for v, s in zip(u0s, s0s):
            rhs = s["Phi"] + (lam * np.log(s["omr"])).sum() - s["eps"] - mu @ (v[:B] + np.log(s["omr"] / s["rho"])) - e * (v[B] - np.log(s["eps"]))
            worst = max(worst, Hval - S0 - rhs)
print("   grid points of P: %d, with finite H: %d; max of LHS - RHS over 40 tangent points each: %+.3e" % (ntot, nfin, worst))
