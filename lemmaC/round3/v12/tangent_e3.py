"""V12: E3 redone with a convergence diagnosis.  The first E3 run reported LHS-RHS = +0.86 using BFGS values that
were NOT converged (overflow warnings).  Here every grid point mu of P is classified:
  interior  : BFGS reaches |grad Phi(u) - mu| < 1e-7  -> H(mu) = value, inequality (4.1) tested;
  otherwise : the objective is followed along the ray c * u_last (c = 1, 2, 4, ..., clipped) -- if it keeps
              decreasing, H(mu) = -inf or a boundary limit; the smallest value seen is still an UPPER bound for H,
              and (4.1) is tested with it."""
import os, itertools, warnings
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "tangent_adv.py")).read()
exec(src[:src.index("TB = float")])
warnings.filterwarnings("ignore")
for B, lam in ((1, [0.7]), (1, [3.0]), (2, [0.5, 1.5]), (2, [2.0, 0.2])):
    mdl = Model(lam); lam = np.array(lam)
    u0s = [rng.uniform(-6, 6, B + 1) for _ in range(40)]; s0s = [mdl.stats(v) for v in u0s]
    worst_in = -np.inf; worst_out = -np.inf; nin = nout = 0; worst_out_at = None
    for fr in itertools.product(np.linspace(0.0, 1.0, 13), repeat=B):
        mu = lam * np.array(fr); e = 1 - mu.sum()
        if e < 0: continue
        f = lambda v: (lambda s: s["Phi"] - mu @ v[:B] - e * v[B])(mdl.stats(np.clip(v, -300, 300)))
        r = minimize(f, np.zeros(B + 1), method="BFGS", options=dict(gtol=1e-11, maxiter=2000))
        s = mdl.stats(np.clip(r.x, -300, 300)); res = np.abs(np.concatenate([s["mu"] - mu, [s["eps"] - e]])).max()
        Hup = r.fun
        if res >= 1e-7:
            v = r.x - r.x.mean()
            for c in (1, 2, 4, 8, 16, 32, 64, 128, 256):
                val = f(c * v)
                if np.isfinite(val): Hup = min(Hup, val)
        S0 = sum(l * (0 if x in (0.0, 1.0) else -(x * np.log(x) + (1 - x) * np.log(1 - x))) for l, x in zip(lam, fr)) + (e - e * np.log(e) if e > 0 else 0.0)
        g = max(Hup - S0 - (s0["Phi"] + (lam * np.log(s0["omr"])).sum() - s0["eps"] - mu @ (v0[:B] + np.log(s0["omr"] / s0["rho"])) - e * (v0[B] - np.log(s0["eps"])))
                for v0, s0 in zip(u0s, s0s))
        if res < 1e-7:
            nin += 1
            if g > worst_in: worst_in, win_at = g, (fr, res)
        else:
            nout += 1
            if g > worst_out: worst_out, worst_out_at = g, (fr, Hup, res)
    print("lam=%s  P(D=0)=%.4f: interior points %d, max LHS-RHS = %+.3e | non-converged (boundary or H=-inf) %d, max (upper bound for H) - S0 - RHS = %+.3e at %s"
          % (lam, np.exp(-lam.sum()), nin, worst_in, nout, worst_out, worst_out_at), "| worst interior at (fractions mu_b/lam_b, residual):", win_at, flush=True)
