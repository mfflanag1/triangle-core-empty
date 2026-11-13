"""V12, second pass (noise-aware).  Uses Model / gap / eig_ratio from tangent_adv.py.
 E2'  max of [LHS-RHS of (4.1)]/|dmu|^2 subject to |dmu|^2 > 1e-4 (roundoff then < 1e-9), bounded search.
 E4   Proposition 4.5 directly: max over (lambda, ell, u0, u1) of  L(mu(u1)) - theta*beta - Psi_0(u0)   (claim: <= 0),
      ell arbitrary reals (the proof uses nothing about ell).
 E3   H(mu) by direct minimisation on a grid of P."""
import sys, os, time, itertools
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "tangent_adv.py")).read()
exec(src[:src.index("TB = float")])
TB = float(sys.argv[1]) if len(sys.argv) > 1 else 25.0

def run(name, f, bounds, seed):
    t = time.time(); best = [-np.inf, None]; lo = np.array([b[0] for b in bounds]); hi = np.array([b[1] for b in bounds])
    def obj(x):
        x = np.clip(x, lo, hi); v = f(x)
        if not np.isfinite(v): return 1e9
        if v > best[0]: best[0], best[1] = v, x.copy()
        return -v
    differential_evolution(obj, bounds, seed=seed, popsize=12, maxiter=100000, tol=0, mutation=(0.4, 1.2), recombination=0.8, polish=False,
                           callback=lambda xk, convergence=None: time.time() - t > TB)
    minimize(obj, best[1], method="Nelder-Mead", options=dict(maxfev=2500, xatol=1e-10, fatol=1e-15))
    print("%-30s max = %+.6e   at %s   [%.0fs]" % (name, best[0], np.array2string(best[1], precision=3), time.time() - t), flush=True)

def hh(r, omr): return -(r * np.log(r) + omr * np.log(omr))
LL, LH, UM = np.log(1e-2), np.log(6.0), 7.0
for B in (1, 2):
    bl = [(LL, LH)] * B; bu = [(-UM, UM)] * (B + 1)
    def e2(x):
        g, d2 = gap(Model(np.exp(x[:B])), x[B:2 * B + 1], x[2 * B + 1:])
        return g / d2 if d2 > 1e-4 else -1e3
    run("E2' gap/|dmu|^2 (|dmu|>.01) B=%d" % B, e2, bl + bu + bu, 400 + B)
    def e4(x):
        mdl = Model(np.exp(x[:B])); lam = mdl.lam; ell = x[B:2 * B + 1]; u0 = x[2 * B + 1:3 * B + 2]; u1 = x[3 * B + 2:]
        s0 = mdl.stats(u0); s1 = mdl.stats(u1); mu1, e1 = s1["mu"], s1["eps"]
        S1 = (lam * hh(s1["rho"], s1["omr"])).sum() + e1 - e1 * np.log(e1)
        Lval = 3 * (s1["Phi"] - mu1 @ u1[:B] - e1 * u1[B]) - 2 * S1 + ell[:B] @ mu1 + ell[B] * e1
        Psi = 3 * s0["Phi"] + (lam * np.log(s0["omr"] ** 3 + np.exp(ell[:B] - 3 * u0[:B]) * s0["rho"] ** 3)).sum() - 3 * s0["eps"] + np.exp(ell[B] - 3 * u0[B]) * s0["eps"] ** 3
        return Lval - Psi
    run("E4 L(mu1)-Psi0(u0) B=%d" % B, e4, bl + [(-6, 3)] * (B + 1) + bu + bu, 500 + B)

worst = -np.inf; nfin = 0; ntot = 0
for B, lam in ((1, [0.7]), (1, [3.0]), (2, [0.5, 1.5]), (2, [2.0, 0.2])):
    mdl = Model(lam); lam = np.array(lam)
    u0s = [rng.uniform(-6, 6, B + 1) for _ in range(40)]; s0s = [mdl.stats(v) for v in u0s]
    for fr in itertools.product(np.linspace(0.0, 1.0, 13), repeat=B):
        mu = lam * np.array(fr); e = 1 - mu.sum()
        if e < 0: continue
        ntot += 1
        r = minimize(lambda v: mdl.stats(v)["Phi"] - mu @ v[:B] - e * v[B], np.zeros(B + 1), method="BFGS", options=dict(gtol=1e-10, maxiter=400))
        if r.fun < -30 or not np.isfinite(r.fun): continue      # H = -inf: nothing to check
        nfin += 1
        S0 = sum(l * (0 if x in (0.0, 1.0) else -(x * np.log(x) + (1 - x) * np.log(1 - x))) for l, x in zip(lam, fr)) + (e - e * np.log(e) if e > 0 else 0.0)
        for v, s in zip(u0s, s0s):
            rhs = s["Phi"] + (lam * np.log(s["omr"])).sum() - s["eps"] - mu @ (v[:B] + np.log(s["omr"] / s["rho"])) - e * (v[B] - np.log(s["eps"]))
            worst = max(worst, r.fun - S0 - rhs)
print("E3 grid points of P: %d, finite H (BFGS value, an overestimate of H): %d; max LHS-RHS over 40 tangent points each: %+.3e" % (ntot, nfin, worst))
