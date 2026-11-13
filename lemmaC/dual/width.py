"""Core width on nonempty-core U[0,1] triangle markets (gate-1 seeds).  The core is invariant under the
sector shifts p -> p + (aX 1_X + aY 1_Y + aZ 1_Z), aX+aY+aZ = 0 (while p >= 0), so width is measured
modulo shifts: for each agent i in sector S, the range over the core of p_i - mean_S(p).  Reports
max and median agent range, the dual depth delta* = max{t : p(e) >= v(e) + t for e not in M, p in
core-hyperplane}  (= -mu* of Lemma 2), and counts of non-M teams with slack <= eps at the analytic-ish
centre (midpoint of the depth LP).  Floating point; evidence only."""
import json, sys, os, argparse
import numpy as np
from scipy.optimize import linprog
from scipy import sparse as sp
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'gate1'))
from gate1 import instance, structure

def run(args):
    q, idx = args
    A, B, C = instance('triangle', 'unif', q, idx)
    v, Mm = structure(A, B, C); n = 3 * q; Mm = Mm.tocsr()
    lp = linprog(-v, A_eq=Mm.T.tocsr(), b_eq=np.ones(n), bounds=(0, None), method='highs-ds')
    x = lp.x
    if np.max(np.minimum(x, 1 - x)) > 1e-6:
        return dict(q=q, idx=idx, empty=True)
    inM = x > 0.5; T = float(v @ x)
    # core: p >= 0, p(e) >= v(e) all e, p(e) = v(e) on M (sum p = T then automatic)
    Aub = -Mm[~inM]; bub = -v[~inM]; Aeq = Mm[inM]; beq = v[inM]
    rng_ = []
    for s in range(3):
        for i in range(q):
            c = np.zeros(n); c[s * q:(s + 1) * q] = -1.0 / q; c[s * q + i] += 1.0
            lo = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=(0, None), method='highs')
            hi = linprog(-c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=(0, None), method='highs')
            rng_.append(-hi.fun - lo.fun)
    rng_ = np.array(rng_)
    # depth LP: max t s.t. p(e) - t >= v(e) off M, p(e) = v(e) on M, p >= 0
    Aub2 = sp.hstack([-Mm[~inM], np.ones(((~inM).sum(), 1))]).tocsr()
    Aeq2 = sp.hstack([Mm[inM], np.zeros((q, 1))]).tocsr()
    c2 = np.zeros(n + 1); c2[-1] = -1
    d = linprog(c2, A_ub=Aub2, b_ub=bub, A_eq=Aeq2, b_eq=beq, bounds=[(0, None)] * n + [(None, None)], method='highs')
    p = d.x[:n]; slack = (Mm[~inM] @ p) - v[~inM]
    eps = [0.01, 0.03, 0.1]
    return dict(q=q, idx=idx, empty=False, T=T, wmax=float(rng_.max()), wmed=float(np.median(rng_)),
                wmean=float(rng_.mean()), depth=float(-d.fun),
                nslack={str(e): int((slack <= e).sum()) for e in eps},
                gapM=float(3 - v[inM].min()), gapMmean=float(3 - v[inM].mean()))

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--n', type=int)
    a = ap.parse_args()
    for q in a.q:
        with ProcessPoolExecutor(4) as ex, open(os.path.join(HERE, 'out', f'width_q{q:02d}.jsonl'), 'w') as f:
            rs = list(ex.map(run, [(q, i) for i in range(a.n)]))
            for r in rs: f.write(json.dumps(r) + '\n')
        ne = [r for r in rs if not r['empty']]
        if ne:
            m = lambda k: np.mean([r[k] for r in ne]); ns = lambda e: np.mean([r['nslack'][e] for r in ne])
            print(q, 'nonempty %d/%d' % (len(ne), len(rs)), 'wmax %.4f wmed %.4f depth %.5f gapMmean %.4f nslack.01 %.1f .03 %.1f .1 %.1f'
                  % (m('wmax'), m('wmed'), m('depth'), m('gapMmean'), ns('0.01'), ns('0.03'), ns('0.1')), flush=True)
