"""Large-q LP constant via near-tight truncation (xi <= K), pair-additive vs iid-team (Poissonised) model.
Restricted LP = upper bound on the LP deficit; K-convergence is checked by re-solving with smaller K.
Output: (sum xi x)/q = (3q-LP*)/q^(1/3). Floating point, evidence only."""
import numpy as np, sys, time, json
from scipy.sparse import csc_matrix
from scipy.optimize import linprog

def gen_pair(q, K, rng):
    s = K * q ** (-2 / 3)
    B = rng.random((q, q), dtype=np.float32)          # d(b,c)
    ta, tb, tc, xi = [], [], [], []
    for a in range(q):
        ra = rng.random(q, dtype=np.float32); rc = rng.random(q, dtype=np.float32)   # d(a,b), d(a,c)
        nb = np.nonzero(ra < s)[0]; nc = np.nonzero(rc < s)[0]
        if len(nb) == 0 or len(nc) == 0: continue
        tot = ra[nb][:, None] + rc[nc][None, :] + B[np.ix_(nb, nc)]
        ib, ic = np.nonzero(tot < s)
        ta.append(np.full(len(ib), a)); tb.append(nb[ib]); tc.append(nc[ic]); xi.append(tot[ib, ic])
    return np.concatenate(ta), np.concatenate(tb), np.concatenate(tc), np.concatenate(xi).astype(float) * q ** (2 / 3)

def gen_iid(q, K, rng):
    n = rng.binomial(q ** 3, K ** 3 / (6 * q ** 2)) if q < 2000 else rng.poisson(q * K ** 3 / 6)
    ta = rng.integers(0, q, n); tb = rng.integers(0, q, n); tc = rng.integers(0, q, n)
    xi = K * rng.random(n) ** (1 / 3)
    return ta, tb, tc, xi

METHOD = 'highs-ipm'
def solve(q, ta, tb, tc, xi):
    n = len(xi)
    rows = np.concatenate([ta, q + tb, 2 * q + tc]); cols = np.tile(np.arange(n), 3)
    A = csc_matrix((np.ones(3 * n), (rows, cols)), shape=(3 * q, n))
    r = linprog(xi, A_eq=A, b_eq=np.ones(3 * q), bounds=(0, None), method=METHOD)
    if r.status != 0: return None
    x = r.x; sup = x > 1e-7
    load_frac = np.zeros(q); np.add.at(load_frac, ta[sup & (x < 1 - 1e-7)], x[sup & (x < 1 - 1e-7)])
    return dict(val=float(r.fun / q), maxxi=float(xi[sup].max()), nsup=int(sup.sum()), frac=float((load_frac > 1e-7).mean()))

if __name__ == '__main__':
    model, q, K, nrep, seed = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    Ks = [float(k) for k in sys.argv[6].split(',')] if len(sys.argv) > 6 else []
    for i in range(nrep):
        rng = np.random.default_rng([seed, q, i]); t0 = time.time()
        ta, tb, tc, xi = (gen_pair if model == 'pair' else gen_iid)(q, K, rng)
        t1 = time.time(); out = dict(model=model, q=q, K=K, i=i, nteams=len(xi), gen_s=round(t1 - t0, 1))
        res = solve(q, ta, tb, tc, xi); out['res'] = res; out['lp_s'] = round(time.time() - t1, 1)
        for k in Ks:
            m = xi <= k; rk = solve(q, ta[m], tb[m], tc[m], xi[m]); out['val_K%g' % k] = None if rk is None else rk['val']
        print(json.dumps(out), flush=True)
