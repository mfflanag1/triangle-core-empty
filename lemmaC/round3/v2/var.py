"""V2 referee, front 3: independent evaluation of Sigma[rho], e[rho] (torch, autograd),
Monte Carlo cross-check of the quadrature, own solve of e_1*, Hessian test, and a
global search for profiles with Sigma >= 0 and e < 2.3025."""
import sys, math, numpy as np, torch
from scipy.optimize import minimize
torch.set_default_dtype(torch.float64)

class Model:
    def __init__(s, Xmax=14.0, nx=280, nt=1000, slo=-30.0, shi=55.0):
        s.xi = (torch.arange(nx) + 0.5) * Xmax / nx
        s.lam = s.xi ** 2 / 2 * (Xmax / nx)
        sg = torch.linspace(slo, shi, nt); s.t = torch.exp(sg); s.ds = (sg[1] - sg[0]).item()
    def stats(s, u):
        u = u - u.max().detach()
        w = torch.exp(u); tw = s.t[:, None] * w[None, :]; E = torch.exp(-tw)
        phi = torch.exp(-((-torch.expm1(-tw)) @ s.lam))
        rho = (tw * E * phi[:, None]).sum(0) * s.ds
        ElnW = (torch.exp(-s.t) - phi).sum() * s.ds
        return u, rho, ElnW
    def SE(s, u):
        u, rho, ElnW = s.stats(u); rho = rho.clamp(1e-300, 1 - 1e-13)
        H = ElnW - (u * rho * s.lam).sum()
        h = -(rho * torch.log(rho) + (1 - rho) * torch.log1p(-rho))
        return 3 * H - 2 * (h * s.lam).sum(), (s.xi * rho * s.lam).sum(), rho
    def np_fun(s, f):
        def g(x):
            u = torch.tensor(x, requires_grad=True); v = f(u); v.backward()
            return v.item(), u.grad.numpy().copy()
        return g

def lag_max(M, beta, u0, maxiter=3000):
    f = M.np_fun(lambda u: -(M.SE(u)[0] - 3 * beta * M.SE(u)[1]))
    r = minimize(f, u0, jac=True, method="L-BFGS-B", options=dict(maxiter=maxiter, maxfun=maxiter * 2, ftol=1e-15, gtol=1e-10))
    return -r.fun, r.x

def mc_check(M, u, N=200000, seed=3):
    rng = np.random.default_rng(seed); un = (u - u.max()); w = np.exp(un); lam = M.lam.numpy(); xi = M.xi.numpy()
    keep = xi < 8.0
    W = np.zeros(N)
    for i in np.where(keep)[0]: W += rng.poisson(lam[i], size=N) * w[i]
    ElnW = np.log(W[W > 0]).mean(); se = np.log(W[W > 0]).std() / math.sqrt(N)
    rho = np.array([(w[i] / (w[i] + W)).mean() for i in range(len(w))])
    return ElnW, se, rho, (W == 0).mean()

if __name__ == "__main__":
    mode = sys.argv[1]
    M = Model()
    xi = M.xi.numpy()
    if mode == "solve":
        # 1. own stationary solve: maximise L_beta, secant on beta for Sigma = 0
        u = -1.5 * xi
        def at(beta, u):
            L, u = lag_max(M, beta, u); S, e, rho = M.SE(torch.tensor(u)); return L, u, S.item(), e.item(), rho.numpy()
        b0, b1 = 0.54, 0.58
        L0, u, S0, e0, _ = at(b0, u); L1, u, S1, e1, _ = at(b1, u)
        for _ in range(6):
            b2 = b1 - S1 * (b1 - b0) / (S1 - S0); b0, S0 = b1, S1; b1 = b2; L1, u, S1, e1, rho = at(b1, u)
        print("grid nx=%d: beta*=%.5f e*=%.5f Sigma=%.2e L*=%.6f rho(0)=%.4f norm=%.6f" % (len(xi), b1, e1, S1, L1, rho[0], (rho * M.lam.numpy()).sum()))
        np.save("ustar.npy", np.array([b1] + list(u)))
        for (nx, nt) in ((560, 1400), (1400, 2000)):
            M2 = Model(nx=nx, nt=nt); u2 = np.interp(M2.xi.numpy(), xi, u)
            L2, u2 = lag_max(M2, b1, u2, 400); S, e, r = M2.SE(torch.tensor(u2))
            print("grid nx=%d at same beta: e=%.5f Sigma=%+.3e L=%.6f" % (nx, e.item(), S.item(), L2))
        # 2. Monte Carlo cross-check of the Laplace quadrature at the solution
        ElnW, se, rmc, pz = mc_check(M, u)
        _, rq, Eq = M.stats(torch.tensor(u))
        print("MC check: E ln W quad=%.5f MC=%.5f+-%.5f | max|rho_quad-rho_MC| (xi<6)=%.4f" % (Eq.item(), ElnW, se, np.abs(rq.numpy() - rmc)[xi < 6].max()))
        # 3. Hessian of L in u at the solution (shift direction is null)
        ut = torch.tensor(u)
        Hs = torch.autograd.functional.hessian(lambda v: M.SE(v)[0] - 3 * b1 * M.SE(v)[1], ut).numpy()
        ev = np.linalg.eigvalsh((Hs + Hs.T) / 2)
        print("Hessian eigenvalues of L in ln w: max=%.3e, 2nd max=%.3e, min=%.3e (max should be ~0: scale mode)" % (ev[-1], ev[-2], ev[0]))
        # rho-space criterion: lambda_max(R^-1/2 M R^-1/2) <= 3/2, M = d rho/d u (as a measure), R = rho(1-rho)lam
        J = torch.autograd.functional.jacobian(lambda v: M.SE(v)[2] * M.lam, ut).numpy()
        R = (rq.numpy() * (1 - rq.numpy()) * M.lam.numpy()); keep = xi < 9
        A = J[np.ix_(keep, keep)] / np.sqrt(np.outer(R[keep], R[keep]))
        print("lambda_max(R^-1/2 M R^-1/2)=%.4f (local max of Sigma-3 beta e in rho iff < 1.5)" % np.linalg.eigvalsh((A + A.T) / 2)[-1])
    if mode == "search":
        d = np.load("ustar.npy"); beta, ustar = d[0], d[1:]
        Ls = (M.SE(torch.tensor(ustar))[0] - 3 * beta * M.SE(torch.tensor(ustar))[1]).item()
        print("reference L*=%.6f at beta=%.5f" % (Ls, beta))
        rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
        knots = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 6, 9, 14.0])
        best = -1e9; nbetter = 0; ends = []
        for trial in range(int(sys.argv[3]) if len(sys.argv) > 3 else 40):
            kind = trial % 4
            if kind == 0: u0 = np.interp(xi, knots, rng.uniform(-30, 5, size=len(knots)))
            elif kind == 1: u0 = -rng.uniform(0.2, 6) * xi ** rng.choice([1, 2]) + rng.normal(0, 3, size=len(xi))
            elif kind == 2: c = rng.uniform(0.3, 3); u0 = np.where(xi < c, rng.uniform(0, 25), 0) - rng.uniform(0.3, 4) * xi
            else: u0 = np.where(np.floor(xi / rng.uniform(0.1, 1)) % 2 == 0, 0.0, -rng.uniform(2, 30)) - rng.uniform(0.3, 4) * xi
            L0 = (M.SE(torch.tensor(u0))[0] - 3 * beta * M.SE(torch.tensor(u0))[1]).item()
            L, u = lag_max(M, beta, u0, 1500); S, e, rho = M.SE(torch.tensor(u))
            dist = np.abs((u - u.max()) - (ustar - ustar.max()))[xi < 5].max()
            ends.append((L - Ls, e.item(), S.item(), dist)); best = max(best, L - Ls); nbetter += (L - Ls > 1e-7)
            print("trial %2d kind %d: L0-L*=%+.3e -> L-L*=%+.3e e=%.5f Sigma=%+.4f max|u-u*|(xi<5)=%.2e" % (trial, kind, L0 - Ls, L - Ls, e.item(), S.item(), dist), flush=True)
        print("SUMMARY: best L-L*=%+.3e ; trials ending above L*: %d ; trials ending at a different point (dist>1e-2): %d" % (best, nbetter, sum(1 for x in ends if x[3] > 1e-2)))
