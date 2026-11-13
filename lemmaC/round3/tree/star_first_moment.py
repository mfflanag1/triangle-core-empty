"""Star-constrained (depth-1, configuration-model) first moment for integral perfect matchings
of the Poisson 3-hypertree environment.  Evidence only (floating point quadrature).

Sigma[rho] = 3*H - 2*int h(rho) dLambda,  Lambda(dxi) = xi^2/2 dxi,
H = E ln W - int ln w * rho dLambda,  W = sum_{xi in Pi} w(xi),  rho(xi) = E[w/(w+W)].
Stationarity under a cost multiplier: ln w = -beta*xi + (2/3) ln(rho/(1-rho)).
Unconstrained (annealed) exponent for the same rho: -2 - int rho ln rho dLambda.
"""
import numpy as np, sys
Xi, nx = 14.0, 1400
xi = (np.arange(nx) + 0.5) * Xi / nx
lam = xi**2 / 2 * (Xi / nx)
u = np.linspace(-25.0, 45.0, 2800); du = u[1] - u[0]; t = np.exp(u)

def stats(lnw):
    w = np.exp(lnw)
    tw = np.outer(t, w)
    phi = np.exp(-((-np.expm1(-tw)) @ lam))
    rho = ((tw * np.exp(-tw)) * phi[:, None]).sum(0) * du
    ElnW = (np.exp(-t) - phi).sum() * du
    return rho, ElnW

def solve(beta, iters=400, damp=0.5):
    lnw = -3 * beta * xi
    for it in range(iters):
        rho, _ = stats(lnw)
        rho = np.clip(rho, 1e-300, 1 - 1e-12)
        new = -beta * xi + (2.0 / 3.0) * (np.log(rho) - np.log1p(-rho))
        new -= new.max()
        diff = np.abs(new - lnw)[xi < 8].max()
        lnw = damp * new + (1 - damp) * lnw
        if diff < 1e-9: break
    rho, ElnW = stats(lnw)
    rho = np.clip(rho, 1e-300, 1 - 1e-12)
    norm = (rho * lam).sum()
    H = ElnW - (lnw * rho * lam).sum()
    h = -(rho * np.log(rho) + (1 - rho) * np.log1p(-rho))
    Sig = 3 * H - 2 * (h * lam).sum()
    e = (xi * rho * lam).sum()
    ann = -2 - (rho * np.log(rho) * lam).sum()
    return dict(beta=beta, e=e, Sigma=Sig, ann_same_rho=ann, ann_at_e=np.log(e**3 * np.e / 27), norm=norm, it=it, diff=diff, rho0=rho[0])

if __name__ == "__main__":
    for beta in [float(b) for b in sys.argv[1:]] or [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.70]:
        r = solve(beta)
        print(" ".join(f"{k}={v:.6g}" for k, v in r.items()), flush=True)
