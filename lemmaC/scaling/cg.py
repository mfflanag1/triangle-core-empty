"""Gap D = LP* - T scaling on the random triangle market (uniform pair values). Floating point; evidence only.
LP* by column generation over teams: restricted LP (equality, highs-ds) + exact dense pricing over all q^3 teams
(r_f = v_f - p_x - p_y - p_z), stop when max r <= PTOL. The final dual is feasible for the full LP up to PTOL.
T by reduced-cost fixing (as in ../structure/mustruct.py): MILP on {r_f >= -delta}, double delta until
T_inc >= LP* - delta.
Seeds: tag=None -> gate-1 instance (validation); tag=k -> SeedSequence([20261111, 0, 0, q, idx, k])."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import json, sys, time, argparse
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
TAG = 7
PTOL = 1e-9

def draw(q, idx, tag):
    ent = [20261111, 0, 0, q, idx] + ([] if tag is None else [tag])
    rng = np.random.default_rng(np.random.SeedSequence(ent))
    return rng.random((q, q)), rng.random((q, q)), rng.random((q, q))

def teams_cols(idx, q):
    a, b, c = np.unravel_index(idx, (q, q, q))
    return np.stack([a, q + b, 2 * q + c], 1)

def amat(idx, q):
    cols = teams_cols(idx, q).reshape(-1)
    rows = np.repeat(np.arange(len(idx)), 3)
    return sp.csr_matrix((np.ones(3 * len(idx)), (cols, rows)), shape=(3 * q, len(idx)))

def solve_lp(V, q, k0=3, maxadd=None):
    n = 3 * q; vf = V.reshape(-1)
    maxadd = maxadd or 20 * q
    # initial columns: top-k0 third member for each pair in each of the three directions, plus the diagonal
    init = [np.ravel_multi_index((np.arange(q),) * 3, (q, q, q))]
    top = np.argsort(-V, axis=2)[:, :, :k0]                      # (x,y) -> z
    xs, ys = np.meshgrid(np.arange(q), np.arange(q), indexing='ij')
    init.append(np.ravel_multi_index((np.repeat(xs[..., None], k0, 2), np.repeat(ys[..., None], k0, 2), top), (q, q, q)).ravel())
    top = np.argsort(-V, axis=1)[:, :k0, :]                      # (x,z) -> y
    xs, zs = np.meshgrid(np.arange(q), np.arange(q), indexing='ij')
    init.append(np.ravel_multi_index((np.repeat(xs[:, None, :], k0, 1), top, np.repeat(zs[:, None, :], k0, 1)), (q, q, q)).ravel())
    top = np.argsort(-V, axis=0)[:k0]                            # (y,z) -> x
    ys, zs = np.meshgrid(np.arange(q), np.arange(q), indexing='ij')
    init.append(np.ravel_multi_index((top, np.repeat(ys[None], k0, 0), np.repeat(zs[None], k0, 0)), (q, q, q)).ravel())
    cols = np.unique(np.concatenate(init))
    it = 0
    while True:
        it += 1
        lp = linprog(-vf[cols], A_eq=amat(cols, q), b_eq=np.ones(n), bounds=(0, None), method='highs-ds')
        assert lp.status == 0, lp.message
        p = -lp.eqlin.marginals
        R = V - p[:q, None, None] - p[None, q:2 * q, None] - p[None, None, 2 * q:]
        r = R.reshape(-1)
        pos = np.nonzero(r > PTOL)[0]
        if len(pos) == 0:
            return -lp.fun, p, r, cols, lp.x, it
        if len(pos) > maxadd:
            pos = pos[np.argpartition(-r[pos], maxadd)[:maxadd]]
        cols = np.union1d(cols, pos)

def integer_opt(vf, q, LPs, r, delta0):
    n = 3 * q; delta = delta0; tries = []
    while True:
        cand = np.nonzero(r >= -delta)[0]
        mi = milp(-vf[cand], constraints=LinearConstraint(amat(cand, q), np.ones(n), np.ones(n)),
                  integrality=np.ones(len(cand)), bounds=Bounds(0, 1), options=dict(time_limit=900, mip_rel_gap=1e-12))
        tries.append((delta, len(cand), int(mi.status)))
        if mi.status == 0 and -mi.fun >= LPs - delta + 1e-9:
            return -mi.fun, tries
        if delta > 3: raise RuntimeError('delta blew up')
        delta *= 2

def run(args):
    q, idx, tag, full = args
    t0 = time.time()
    A, B, C = draw(q, idx, tag)
    V = A[:, :, None] + B[:, None, :] + C[None, :, :]
    vf = V.reshape(-1)
    LPs, p, r, cols, x, it = solve_lp(V, q)
    t1 = time.time()
    rec = dict(q=q, idx=idx, tag=tag, LP=LPs, ncols=int(len(cols)), cg_iters=it, lp_secs=t1 - t0)
    frac = np.max(np.minimum(x, 1 - x)) < 1e-6
    if frac:
        T = float(vf[cols[x > 0.5]].sum()); rec['tries'] = []
    else:
        T, rec['tries'] = integer_opt(vf, q, LPs, r, 0.15 if q <= 48 else 0.3)
    rec.update(T=T, gap=LPs - T, lp_integral=bool(frac), secs=time.time() - t0)
    if full:  # validation: full LP over all q^3 teams
        allc = np.arange(q ** 3)
        lpf = linprog(-vf, A_eq=amat(allc, q), b_eq=np.ones(3 * q), bounds=(0, None), method='highs-ds')
        rec['LP_full'] = -lpf.fun; rec['full_secs'] = time.time() - t0 - rec['secs']
    return rec

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--n', type=int); ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--workers', type=int, default=4); ap.add_argument('--gate1', action='store_true')
    ap.add_argument('--full', action='store_true'); ap.add_argument('--out', default=None)
    a = ap.parse_args()
    tag = None if a.gate1 else TAG
    for q in a.q:
        fn = a.out or os.path.join(HERE, 'out', f"{'val' if a.gate1 else 'scal'}_q{q:03d}.jsonl")
        with ProcessPoolExecutor(a.workers) as ex, open(fn, 'a') as f:
            for rr in ex.map(run, [(q, i, tag, a.full) for i in range(a.start, a.start + a.n)]):
                f.write(json.dumps(rr) + '\n'); f.flush()
                print(q, rr['idx'], round(rr['gap'], 4), round(rr['secs'], 1), flush=True)
