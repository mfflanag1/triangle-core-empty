"""V16 independent check: numerical fixed point of X = min_k (xi_k - X'_k - X''_k),
Poisson intensity xi^2/2.  P(X > x) = exp(-G(x)), G(x) = E[(x+Y)_+^3]/6, Y = X'+X''.
Own code; no lane files read."""
import numpy as np

def G_from_py(x, py, y, K=0.0):
    """G_K(x) = sum_{y: x+y > K} py * ((x+y)^3 - K^3)/6 (K=0 gives G)."""
    # reverse cumulative sums of py*y^j
    r0 = np.concatenate([np.cumsum(py[::-1])[::-1], [0.0]])
    r1 = np.concatenate([np.cumsum((py*y)[::-1])[::-1], [0.0]])
    r2 = np.concatenate([np.cumsum((py*y*y)[::-1])[::-1], [0.0]])
    r3 = np.concatenate([np.cumsum((py*y**3)[::-1])[::-1], [0.0]])
    j = np.searchsorted(y, K - x, side='right')  # first index with y > K - x
    return (x**3*r0[j] + 3*x*x*r1[j] + 3*x*r2[j] + r3[j] - K**3*r0[j]) / 6.0

def solve(L=8.0, h=0.01, iters=400, damp=0.5, verbose=False):
    x = np.arange(-L, L + h/2, h)
    G = np.maximum(x + 0.5, 0)**3 / 6
    xm = (x[:-1] + x[1:]) / 2
    for it in range(iters):
        surv = np.exp(-G)
        pm = surv[:-1] - surv[1:]
        pm = np.maximum(pm, 0); pm /= pm.sum()
        py = np.convolve(pm, pm)
        y = 2*xm[0] + h*np.arange(len(py))
        Gn = G_from_py(x, py, y)
        err = np.max(np.abs(np.exp(-Gn) - np.exp(-G)))
        G = (1-damp)*G + damp*Gn
        if verbose and it % 20 == 0:
            print(it, err)
        if err < 1e-13:
            break
    surv = np.exp(-G)
    pm = surv[:-1] - surv[1:]; pm /= pm.sum()
    py = np.convolve(pm, pm)
    y = 2*xm[0] + h*np.arange(len(py))
    ps = np.convolve(py, pm)
    s = 3*xm[0] + h*np.arange(len(ps))
    return dict(x=x, G=G, xm=xm, pm=pm, y=y, py=py, s=s, ps=ps, err=err, it=it, h=h)

if __name__ == "__main__":
    for h in (0.02, 0.01, 0.005):
        d = solve(h=h, verbose=False)
        EX = (d['pm']*d['xm']).sum()
        sp = np.maximum(d['s'], 0)
        ES3 = (d['ps']*sp**3).sum(); ES4 = (d['ps']*sp**4).sum()
        # E[X G(X)] identity: E S_+^4 = 18 E[X G(X)]
        Gm = np.interp(d['xm'], d['x'], d['G'])
        EXG = (d['pm']*d['xm']*Gm).sum()
        print(f"h={h}: iters={d['it']} resid={d['err']:.2e}  EX={EX:.6f}  9EX={9*EX:.6f}  "
              f"ES+^3={ES3:.6f}  ES+^4/8={ES4/8:.6f}  18E[XG(X)]/8={18*EXG/8:.6f}  "
              f"P(X>0)={np.exp(-np.interp(0,d['x'],d['G'])):.5f} P(S>3)={(d['ps'][d['s']>3]).sum():.4e} "
              f"P(S>5)={(d['ps'][d['s']>5]).sum():.3e}")
