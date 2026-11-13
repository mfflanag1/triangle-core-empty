"""Monte Carlo volume of the core polytope by the radial formula Vol = vol(B_n) * E_theta[r(theta)^n] from the Chebyshev centre.
Rescaled units (q^{2/3} per coordinate). Heavy-tailed estimator: evidence only; validated against qhull at q = 3, 4."""
import sys, json, numpy as np
from scipy.optimize import linprog
from scipy.special import gammaln, logsumexp
import corepoly as cp

def logvol(q, idx, K=200000, chunk=20000, seed=0):
    cs = cp.core_system(q, idx)
    if cs is None: return None
    G, h, _, _ = cs; n = G.shape[1]; sc = q ** (2 / 3)
    nrm = np.linalg.norm(G, axis=1)
    cc = linprog(np.r_[np.zeros(n), -1.0], A_ub=np.c_[G, nrm], b_ub=h, bounds=[(None, None)] * n + [(0, None)], method='highs-ds')
    if cc.x[-1] <= 1e-12: return None
    sl = h - G @ cc.x[:n]
    rng = np.random.default_rng([seed, q, idx]); logs = []
    for _ in range(K // chunk):
        th = rng.normal(size=(chunk, n)); th /= np.linalg.norm(th, axis=1, keepdims=True)
        den = th @ G.T
        r = np.min(np.where(den > 1e-14, sl[None, :] / np.maximum(den, 1e-14), np.inf), axis=1)
        logs.append(n * np.log(r * sc))
    logs = np.concatenate(logs)
    lv = n / 2 * np.log(np.pi) - gammaln(n / 2 + 1) + logsumexp(logs) - np.log(len(logs))
    half = [n / 2 * np.log(np.pi) - gammaln(n / 2 + 1) + logsumexp(x) - np.log(len(x)) for x in np.array_split(logs, 2)]
    return float(lv), float(half[0] - half[1])

if __name__ == '__main__':
    q, N = int(sys.argv[1]), int(sys.argv[2]); K = int(sys.argv[3]) if len(sys.argv) > 3 else 200000
    out = []
    for i in range(N):
        r = logvol(q, i, K=K)
        if r: out.append((i,) + r)
    lv = np.array([o[1] for o in out]); sp = np.array([abs(o[2]) for o in out])
    print(json.dumps(dict(q=q, dim=2 * q - 2, cores=len(out), K=K, logvol_mean=round(float(lv.mean()), 3), se=round(float(lv.std() / np.sqrt(len(lv))), 3),
                          per_dim=round(float(lv.mean()) / (2 * q - 2), 3), split_half_abs_diff_mean=round(float(sp.mean()), 3))), flush=True)
    if q <= 4:   # validation against qhull volumes stored by corepoly (--enum 2)
        ref = {json.loads(l)['idx']: json.loads(l).get('logvol_scaled') for l in open(f'core_q{q:02d}.jsonl')}
        d = [o[1] - ref[o[0]] for o in out if ref.get(o[0]) is not None]
        print('validation vs qhull: mean diff', round(float(np.mean(d)), 4), 'max abs diff', round(float(np.max(np.abs(d))), 4), flush=True)
