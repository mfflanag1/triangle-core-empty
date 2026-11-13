"""Core polytope of a random triangle market when the core is nonempty (D = 0). Floating point; evidence only.
Gauge: relabel so M = identity, eliminate c_k = v_k - a_k - b_k, fix a_0, b_0 (kills the 2-dim lineality).
Free coordinates: a_1..a_{q-1}, b_1..b_{q-1} (dimension 2q-2). Constraints: a_i + b_j + c_k >= v(i,j,k) off M.
Outputs per instance: number of vertices, rescaled volume (units q^{-2/3} per coordinate), widths, tight counts."""
import os, sys, json, time, argparse
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.optimize import linprog
import scipy.sparse as sp
from scipy.spatial import HalfspaceIntersection, ConvexHull
from concurrent.futures import ProcessPoolExecutor

def draw(q, idx):
    rng = np.random.default_rng(np.random.SeedSequence([20261111, 3, 77, q, idx]))
    return rng.random((q, q)), rng.random((q, q)), rng.random((q, q))

def amat(q):
    idx = np.arange(q ** 3)
    a, b, c = np.unravel_index(idx, (q, q, q))
    cols = np.stack([a, q + b, 2 * q + c], 1).reshape(-1)
    rows = np.repeat(idx, 3)
    return sp.csr_matrix((np.ones(3 * q ** 3), (cols, rows)), shape=(3 * q, q ** 3))

def core_system(q, idx):
    """returns None if D>0, else (G, h, scale) with core = {y in R^{2q-2}: G y <= h} in gauge coordinates."""
    A, B, C = draw(q, idx)
    V = A[:, :, None] + B[:, None, :] + C[None, :, :]
    lp = linprog(-V.reshape(-1), A_eq=amat(q), b_eq=np.ones(3 * q), bounds=(0, None), method='highs-ds')
    x = lp.x
    if np.max(np.minimum(x, 1 - x)) > 1e-7:
        return None
    sel = np.nonzero(x > 0.5)[0]
    ta, tb, tc = np.unravel_index(sel, (q, q, q))
    # relabel so that M = identity: sector-2 agent tb[i] -> i, sector-3 agent tc[i] -> i (ta sorted = arange)
    o = np.argsort(ta); tb, tc = tb[o], tc[o]
    W = V[:, tb, :][:, :, tc]
    vd = W[np.arange(q), np.arange(q), np.arange(q)]
    p0 = -lp.eqlin.marginals
    a0, b0 = p0[0], p0[q + tb[0]]
    # variables y = (a_1..a_{q-1}, b_1..b_{q-1}); a_0, b_0 fixed. c_k = vd_k - a_k - b_k.
    n = 2 * q - 2
    I, J, K = np.meshgrid(np.arange(q), np.arange(q), np.arange(q), indexing='ij')
    I, J, K = I.ravel(), J.ravel(), K.ravel()
    off = ~((I == J) & (J == K))
    I, J, K = I[off], J[off], K[off]
    m = len(I)
    # constraint: a_I + b_J + vd_K - a_K - b_K >= W[I,J,K]  ->  -(a_I + b_J - a_K - b_K) <= vd_K - W
    Gf = np.zeros((m, 2 * q))
    r = np.arange(m)
    np.add.at(Gf, (r, I), -1.0); np.add.at(Gf, (r, q + J), -1.0)
    np.add.at(Gf, (r, K), 1.0); np.add.at(Gf, (r, q + K), 1.0)
    h = vd[K] - W[I, J, K]
    h = h - Gf[:, 0] * a0 - Gf[:, q] * b0
    G = np.delete(Gf, [0, q], axis=1)
    keep = np.abs(G).sum(1) > 0          # rows with no free variable are automatically satisfied (checked below)
    assert np.all(h[~keep] >= -1e-9)
    return G[keep], h[keep], (I[keep], J[keep], K[keep]), 3 * q - vd.sum()

def analyse(args):
    q, idx, do_enum = args
    t0 = time.time()
    cs = core_system(q, idx)
    if cs is None:
        return dict(q=q, idx=idx, core=False)
    G, h, (I, J, K), defic = cs
    n = G.shape[1]; sc = q ** (2 / 3)
    rec = dict(q=q, idx=idx, core=True, deficit=float(defic))
    # Chebyshev centre
    nrm = np.linalg.norm(G, axis=1)
    cc = linprog(np.r_[np.zeros(n), -1.0], A_ub=np.c_[G, nrm], b_ub=h, bounds=[(None, None)] * n + [(0, None)], method='highs-ds')
    rec['cheb_r_scaled'] = float(cc.x[-1] * sc)
    # widths of the free coordinates (a_i - a_0, b_i - b_0), rescaled
    wid = []; verts = set(); tightsets = []
    for j in range(n):
        e = np.zeros(n); e[j] = 1
        lo = linprog(e, A_ub=G, b_ub=h, bounds=(None, None), method='highs-ds')
        hi = linprog(-e, A_ub=G, b_ub=h, bounds=(None, None), method='highs-ds')
        wid.append((-hi.fun - lo.fun) * sc)
        for s in (lo, hi):
            tight = tuple(np.nonzero(h - G @ s.x < 1e-9)[0].tolist())
            verts.add(tight)
    rec['width_mean'] = float(np.mean(wid)); rec['width_geomean'] = float(np.exp(np.mean(np.log(np.maximum(wid, 1e-300)))))
    rec['lp_vertices_distinct'] = len(verts); rec['lp_solves'] = 2 * n
    rec['tight_union'] = len(set().union(*[set(t) for t in verts]))
    rec['tight_per_vertex'] = float(np.mean([len(t) for t in verts]))
    if do_enum and cc.x[-1] > 1e-11:
        try:
            hs = HalfspaceIntersection(np.c_[G, -h], cc.x[:n])
            P = hs.intersections
            P = P[np.all(np.isfinite(P), axis=1)]
            Pu = np.unique(np.round(P * 1e9).astype(np.int64), axis=0)
            rec['nvert'] = int(len(Pu))
            rec['logvol_scaled'] = float(np.log(ConvexHull(P * sc).volume)) if do_enum > 1 else None
        except Exception as ex:
            rec['enum_err'] = str(ex)[:80]
    rec['secs'] = time.time() - t0
    return rec

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--q', type=int); ap.add_argument('--n', type=int); ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--enum', type=int, default=0); ap.add_argument('--workers', type=int, default=2)
    a = ap.parse_args()
    fn = f'core_q{a.q:02d}.jsonl'
    with ProcessPoolExecutor(a.workers) as ex, open(fn, 'a') as f:
        for r in ex.map(analyse, [(a.q, i, a.enum) for i in range(a.start, a.start + a.n)]):
            f.write(json.dumps(r) + '\n'); f.flush()
