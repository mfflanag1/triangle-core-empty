"""Per-sweep (a_direct, norm) pairs; the exact fixed point has norm = 1, so regress a on
(norm - 1) and report the intercept (drift-corrected energy) with its standard error.
Sweeps are autocorrelated, so the SE is inflated by sqrt(2 tau) with tau from lag-1 corr."""
import json, sys, numpy as np
for f in sys.argv[1:]:
    r = [json.loads(l) for l in open(f) if l.startswith('{')]
    a = np.array([x['a'] for x in r]); m = np.array([x['norm'] for x in r]) - 1
    A = np.vstack([np.ones_like(m), m]).T; c, res, *_ = np.linalg.lstsq(A, a, rcond=None)
    e = a - A @ c; s2 = e.var(ddof=2); cov = s2 * np.linalg.inv(A.T @ A)
    rho = np.corrcoef(e[:-1], e[1:])[0, 1] if len(e) > 3 else 0
    infl = np.sqrt(max(1, (1 + rho) / (1 - rho)))
    extra = ''
    if 'supp_mass' in r[0]:
        extra = f" support teams/agent={np.mean([x['supp_mass'] for x in r]):.3f} full(x=1) teams/agent={np.mean([x['full_mass'] for x in r]):.3f} kink={np.mean([x['kink'] for x in r]):.3f}"
    print(f"{f}: sweeps={len(a)} raw mean a={a.mean():.4f} intercept a(norm=1)={c[0]:.4f} +- {np.sqrt(cov[0,0])*infl:.4f} slope={c[1]:.2f} resid sd={np.sqrt(s2):.4f}{extra}")
