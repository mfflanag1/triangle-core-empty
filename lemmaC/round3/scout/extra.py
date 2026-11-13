"""Supplementary pass on the same instances as scout.py (same seeds): tails of the used pair deficits, and the
Berge triangles in the fractional LP support. Also the iid-U[0,1] sister model by full LP. Floating point; evidence only."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import sys, json, time, argparse
import numpy as np
from scipy.optimize import linprog
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling'))
import cg, scout
TS = [1, 2, 3, 4, 5, 6]

def run(args):
    q, idx = args
    A, B, C = cg.draw(q, idx, scout.TAG)
    V = A[:, :, None] + B[:, None, :] + C[None, :, :]
    vf = V.reshape(-1)
    LPs, p, r, cols, x, it = cg.solve_lp(V, q)
    integral = bool(np.max(np.minimum(x, 1 - x)) < 1e-6)
    if integral: Msel = cols[x > 0.5]
    else: Msel, T = scout.integer_opt_sel(vf, q, LPs, r, 0.15)
    s = q ** (2 / 3)
    mx, my, mz = np.unravel_index(Msel, (q, q, q))
    dM = np.concatenate([1 - A[mx, my], 1 - B[mx, mz], 1 - C[my, mz]]) * s
    sup = cols[x > 1e-7]; w = x[x > 1e-7]
    sx, sy, sz = np.unravel_index(sup, (q, q, q))
    dL = np.concatenate([1 - A[sx, sy], 1 - B[sx, sz], 1 - C[sy, sz]]) * s; wL = np.concatenate([w, w, w])
    rec = dict(q=q, idx=idx, integral=integral,
               tailM=[int((dM > t).sum()) for t in TS], tailL=[float(wL[dL > t].sum()) for t in TS],
               sqM=[float((dM[dM > t] ** 2).sum()) for t in [0] + TS])
    # M pairs kept by the LP: load on each M pair
    Mset = set(zip(mx.tolist(), my.tolist(), mz.tolist()))
    xy = {(a, b) for a, b, c in Mset}; xz = {(a, c) for a, b, c in Mset}; yz = {(b, c) for a, b, c in Mset}
    keep = np.array([((a, b) in xy) + ((a, c) in xz) + ((b, c) in yz) for a, b, c in zip(sx.tolist(), sy.tolist(), sz.tolist())])
    rec['w_keep'] = [float(w[keep == k].sum()) for k in range(4)]   # LP weight on teams sharing k pairs with M (2 impossible)
    if not integral:
        fr = w < 1 - 1e-6
        F = list(zip(sx[fr].tolist(), sy[fr].tolist(), sz[fr].tolist()))
        byxy = {}; byyz = {}; byxz = {}
        for t in F:
            byxy.setdefault((t[0], t[1]), []).append(t); byyz.setdefault((t[1], t[2]), []).append(t); byxz.setdefault((t[0], t[2]), []).append(t)
        ntri = 0; ntri_strong = 0; centreM = 0; centreSup = 0
        Fs = set(F)
        for (a, b), L1 in byxy.items():
            for (b2, c), L2 in byyz.items():
                if b2 != b or (a, c) not in byxz: continue
                L3 = byxz[(a, c)]
                # triples (e1,e2,e3) of distinct teams, e1 has a,b ; e2 has b,c ; e3 has a,c
                for e1 in L1:
                    for e2 in L2:
                        for e3 in L3:
                            if len({e1, e2, e3}) < 3: continue
                            ntri += 1
                            if (a, b, c) not in (e1, e2, e3):
                                ntri_strong += 1
                                centreM += (a, b, c) in Mset
                                centreSup += (a, b, c) in Fs
        rec.update(ntri=ntri, ntri_strong=ntri_strong, centreM=centreM, centreSup=centreSup, nfrac=len(F))
    return rec

def run_u(args):
    q, idx = args
    rng = np.random.default_rng(np.random.SeedSequence([20261111, 770034, q, idx]))
    V = rng.random((q, q, q)); vf = V.reshape(-1)
    lp = linprog(-vf, A_eq=cg.amat(np.arange(q ** 3), q), b_eq=np.ones(3 * q), bounds=(0, None), method='highs-ds')
    return dict(model='U', q=q, idx=idx, LP=-lp.fun, T=None, D=None, integral=bool(np.max(np.minimum(lp.x, 1 - lp.x)) < 1e-6))

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--n', type=int); ap.add_argument('--u', action='store_true')
    a = ap.parse_args()
    for q in a.q:
        fn = os.path.join(HERE, 'out', (f'sister_U_q{q:03d}.jsonl' if a.u else f'extra_q{q:03d}.jsonl'))
        with ProcessPoolExecutor(1) as ex, open(fn, 'a') as f:
            for rr in ex.map(run_u if a.u else run, [(q, i) for i in range(a.n)]):
                f.write(json.dumps(rr) + '\n'); f.flush()
        print('done', q, flush=True)
