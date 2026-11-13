"""V2 referee, front 3: (i) mass evaluation of arbitrary (non-stationary) profiles: is any (e, Sigma) above the
supporting line Sigma = 3 beta* (e - e*)?  is any profile with Sigma >= 0 cheaper than e*?
(ii) direct constrained minimisation of e subject to Sigma >= 0 (SLSQP, autograd) from perturbed starts."""
import sys, numpy as np, torch
from scipy.optimize import minimize
from var import Model
torch.set_default_dtype(torch.float64)
M = Model(); xi = M.xi.numpy(); d = np.load("ustar.npy"); beta, us = d[0], d[1:] - d[1:].max()
def ev(u):
    with torch.no_grad(): S, e, rho = M.SE(torch.tensor(u))
    return S.item(), e.item()
S0, e0 = ev(us); Ls = S0 - 3 * beta * e0
print("reference: beta*=%.5f e*=%.5f Sigma*=%.1e L*=%.6f" % (beta, e0, S0, Ls))
rng = np.random.default_rng(11); knots = np.array([0, 0.3, 0.6, 1, 1.5, 2, 2.5, 3, 4, 5, 7, 10, 14.0])
best = -1e9; emin = 1e9; n = int(sys.argv[1]) if len(sys.argv) > 1 else 3000; cnt = {}
for i in range(n):
    k = i % 6
    if k == 0: u = np.interp(xi, knots, rng.uniform(-30, 5, size=len(knots)))
    elif k == 1: u = us + rng.choice([0.01, 0.1, 0.5, 2, 8, 30]) * np.interp(xi, knots, rng.normal(size=len(knots)))
    elif k == 2: u = np.interp(xi * rng.uniform(0.5, 2.0), xi, us) * rng.uniform(0.5, 2.0)
    elif k == 3: u = us + np.where(np.floor(xi / rng.uniform(0.05, 1.5)) % 2 == 0, 0.0, rng.uniform(-30, 30))
    elif k == 4: c = rng.uniform(0.2, 4); u = us + np.where(xi < c, rng.uniform(-30, 30), 0.0)
    else: c = rng.uniform(0.2, 5); u = us + rng.uniform(-30, 30) * np.exp(-((xi - c) / rng.uniform(0.05, 1)) ** 2)
    S, e = ev(u); dL = S - 3 * beta * e - Ls
    if dL > best: best = dL; arg = (k, S, e)
    if S >= 0 and e < emin: emin = e; arg2 = (k, S, e)
    cnt[k] = max(cnt.get(k, -1e9), dL)
print("random profiles: %d ; max over all of L-L* = %+.3e (family %d, Sigma=%.4f, e=%.4f)" % (n, best, *arg))
print("  per-family max L-L*:", " ".join("%d:%+.2e" % (k, v) for k, v in sorted(cnt.items())))
print("  lowest e among profiles with Sigma>=0: %.5f (family %d, Sigma=%.4f)  vs e*=%.5f" % (emin, arg2[0], arg2[1], e0))
def fe(x):
    u = torch.tensor(x, requires_grad=True); e = M.SE(u)[1]; e.backward(); return e.item(), u.grad.numpy().copy()
def fS(x):
    u = torch.tensor(x, requires_grad=True); S = M.SE(u)[0]; S.backward(); return S.item(), u.grad.numpy().copy()
for j, amp in enumerate((0.5, 2.0, 5.0)):
    u0 = us + amp * np.interp(xi, knots, rng.normal(size=len(knots)))
    r = minimize(lambda x: fe(x)[0], u0, jac=lambda x: fe(x)[1], method="SLSQP",
                 constraints=[dict(type="ineq", fun=lambda x: fS(x)[0], jac=lambda x: fS(x)[1])], options=dict(maxiter=400, ftol=1e-12))
    S, e = ev(r.x); S_i, e_i = ev(u0)
    print("SLSQP start %d (amp %.1f: e=%.4f Sigma=%+.4f) -> e=%.5f Sigma=%+.2e  [%s, %d it]" % (j, amp, e_i, S_i, e, S, r.message[:30], r.nit), flush=True)
