"""Core geometry when D = 0 (pair-additive model). Floating point; evidence only.
After relabelling so that M is the diagonal {(k,k,k)}, a core point is (p,p',p'') with p_k+p'_k+p''_k = v_kkk and
slack(x,y,z) = p_x - p_z + p'_y - p'_z + v_zzz - v_xyz >= 0 for every team. Sector shifts (c1,c2,-c1-c2) act trivially,
so the 'essential core' Q lives in the 2q-2 coordinates u = (p_1..p_{q-1}, p'_1..p'_{q-1}) with p_0 = p'_0 = 0; the
core of the brief (p >= 0) has dimension dim Q + 2.
Measured: eps* = max over Q of the least off-M slack (eps* > 0 <=> no off-M team tight at a relative-interior point
<=> dim Q = 2q-2); number of facets of Q (off-M teams with min slack 0 over Q); vertices of Q by cdd where feasible."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import sys, json, time, argparse
import numpy as np
from scipy.optimize import linprog
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling'))
import cg
TAG = 770033
VMAXQ = 6   # enumerate vertices up to this q

def run(args):
    q, idx = args
    t0 = time.time()
    A, B, C = cg.draw(q, idx, TAG)
    V = A[:, :, None] + B[:, None, :] + C[None, :, :]
    LPs, p, r, cols, x, it = cg.solve_lp(V, q)
    if np.max(np.minimum(x, 1 - x)) >= 1e-6:
        return dict(q=q, idx=idx, core=False)
    mx, my, mz = np.unravel_index(cols[x > 0.5], (q, q, q))
    o = np.argsort(mx); my, mz = my[o], mz[o]
    V = V[:, my, :][:, :, mz]                    # now M is the diagonal
    X, Y, Z = np.meshgrid(np.arange(q), np.arange(q), np.arange(q), indexing='ij')
    X, Y, Z = X.ravel(), Y.ravel(), Z.ravel()
    off = ~((X == Y) & (Y == Z))
    X, Y, Z = X[off], Y[off], Z[off]
    m = len(X); d = 2 * q - 2
    # slack = G u + h0, G row: +e_x - e_z (p block), +e_y - e_z (p' block); coordinate 0 of each block dropped
    G = np.zeros((m, 2 * q))
    np.add.at(G, (np.arange(m), X), 1); np.add.at(G, (np.arange(m), Z), -1)
    np.add.at(G, (np.arange(m), q + Y), 1); np.add.at(G, (np.arange(m), q + Z), -1)
    G = np.delete(G, [0, q], axis=1)
    dz = np.diagonal(np.diagonal(V, axis1=0, axis2=1), axis1=0, axis2=1) if False else V[np.arange(q), np.arange(q), np.arange(q)]
    h0 = dz[Z] - V[X, Y, Z]
    # eps*: max eps s.t. G u + h0 >= eps
    Aub = np.hstack([-G, np.ones((m, 1))]); c = np.zeros(d + 1); c[-1] = -1
    lp = linprog(c, A_ub=Aub, b_ub=h0, bounds=[(None, None)] * d + [(0, None)], method='highs-ds')
    assert lp.status == 0, lp.message
    eps = -lp.fun; uc = lp.x[:d]
    rec = dict(q=q, idx=idx, core=True, eps=float(eps), eps_scaled=float(eps * q ** (2 / 3)),
               tight_at_centre=int((G @ uc + h0 < 1e-9).sum()))
    # coordinate ranges
    lo = np.empty(d); hi = np.empty(d)
    for j in range(d):
        cj = np.zeros(d); cj[j] = 1
        a = linprog(cj, A_ub=-G, b_ub=h0, bounds=(None, None), method='highs-ds'); lo[j] = a.fun
        b = linprog(-cj, A_ub=-G, b_ub=h0, bounds=(None, None), method='highs-ds'); hi[j] = -b.fun
    rec['width_scaled'] = float((hi - lo).mean() * q ** (2 / 3))
    lb = np.where(G > 0, G * lo, G * hi).sum(1) + h0   # valid lower bound on each slack over Q
    cand = np.nonzero(lb <= 1e-9)[0]
    rec['ncand'] = int(len(cand))
    fac = []
    for f in cand:
        a = linprog(G[f], A_ub=-G, b_ub=h0, bounds=(None, None), method='highs-ds')
        if a.fun + h0[f] < 1e-8: fac.append(f)
    fac = np.array(fac, int)
    rec['nfacets'] = int(len(fac))
    # how many facet teams share exactly 2 / 1 / 0 pairs with M? (pairs of the team that are M pairs)
    same = (X[fac] == Y[fac]).astype(int) + (X[fac] == Z[fac]) + (Y[fac] == Z[fac])
    rec['facet_keep1'] = int((same == 1).sum()); rec['facet_keep0'] = int((same == 0).sum())
    rec['facet_def_scaled_mean'] = float((3 - V[X[fac], Y[fac], Z[fac]]).mean() * q ** (2 / 3))
    if q <= VMAXQ:
        import cdd
        Hm = np.hstack([h0[fac, None], G[fac]])
        try:
            mat = cdd.matrix_from_array(Hm, rep_type=cdd.RepType.INEQUALITY)
            gens = np.array(cdd.copy_generators(cdd.polyhedron_from_matrix(mat)).array)
        except AttributeError:   # older pycddlib API
            mat = cdd.Matrix(Hm.tolist(), number_type='float'); mat.rep_type = cdd.RepType.INEQUALITY
            gens = np.array(cdd.Polyhedron(mat).get_generators())
        rec['nvert'] = int((gens[:, 0] == 1).sum()); rec['nrays'] = int((gens[:, 0] == 0).sum())
    rec['secs'] = time.time() - t0
    return rec

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--n', type=int); ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--workers', type=int, default=1)
    a = ap.parse_args()
    os.makedirs(os.path.join(HERE, 'out'), exist_ok=True)
    for q in a.q:
        fn = os.path.join(HERE, 'out', f'core_q{q:03d}.jsonl')
        with ProcessPoolExecutor(a.workers) as ex, open(fn, 'a') as f:
            for rr in ex.map(run, [(q, i) for i in range(a.start, a.start + a.n)]):
                f.write(json.dumps(rr) + '\n'); f.flush()
        print('done', q, flush=True)
