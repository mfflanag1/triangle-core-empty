"""Which off-M teams carry the gap?  Floating point (HiGHS); evidence only.  Same seeds as clone.py.
D     = LP* - T                         (all teams)
D21   = LP(M + pair-keeping teams) - T  (off-teams sharing exactly one pair with M: the three-potential family of REPORT section 4)
D111  = LP(M + teams sharing no pair with M) - T
Also the steepest edge direction at M (cone LP, as in ../../structure/mustruct.py): k = M-teams touched, support, zmax/zmin ratio."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import sys, json, time, argparse
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
HERE = os.path.dirname(os.path.abspath(__file__)); TOL = 1e-7
from clone import amat

def run(q, idx):
    rng = np.random.default_rng(np.random.SeedSequence([20261111, 3, 1, q, idx]))
    A, B, C = rng.random((q, q)), rng.random((q, q)), rng.random((q, q))
    vf = (A[:, :, None] + B[None, :, :] + C[:, None, :]).reshape(-1); Ae = amat(q); n = 3 * q
    lp = linprog(-vf, A_eq=Ae, b_eq=np.ones(n), bounds=(0, None), method='highs-ds'); LPs = -lp.fun
    if np.max(np.minimum(lp.x, 1 - lp.x)) < 1e-6: x = np.round(lp.x); T = LPs
    else:
        mi = milp(-vf, constraints=LinearConstraint(Ae, np.ones(n), np.ones(n)), integrality=np.ones(q ** 3), bounds=Bounds(0, 1),
                  options=dict(mip_rel_gap=0.0)); x = np.round(mi.x); T = -mi.fun
    Mi = np.nonzero(x > 0.5)[0]; a, b, c = np.unravel_index(np.arange(q ** 3), (q, q, q))
    ta = np.empty(q, int); tb = np.empty(q, int); tc = np.empty(q, int)
    ma, mb, mc = np.unravel_index(Mi, (q, q, q)); ta[ma] = np.arange(q); tb[mb] = np.arange(q); tc[mc] = np.arange(q)
    i, j, k = ta[a], tb[b], tc[c]
    nd = (i != j).astype(int) + (j != k) + (i != k)          # 0: M-team, 2: pair-keeping, 3: no M-pair
    out = dict(q=q, idx=idx, D=LPs - T)
    for nm, mask in (('D21', nd <= 2), ('D111', nd != 2)):
        cols = np.nonzero(mask)[0]
        r = linprog(-vf[cols], A_eq=Ae[:, cols], b_eq=np.ones(n), bounds=(0, None), method='highs-ds'); out[nm] = -r.fun - T
    # steepest edge at M: max v.y - v_M.z, loads of off-teams at each agent = z_team, sum z = 1
    sel = x > 0.5; sgn = np.where(sel, -1.0, 1.0)
    Ac = sp.vstack([Ae @ sp.diags(sgn), sp.csr_matrix(sel.astype(float)[None, :])]).tocsr(); rhs = np.zeros(n + 1); rhs[-1] = 1
    cr = linprog(-np.where(sel, -vf, vf), A_eq=Ac, b_eq=rhs, bounds=(0, None), method='highs-ds')
    yz = cr.x; z = yz[sel]; y = yz[~sel]; nz = z[z > TOL]
    out.update(mu=-cr.fun, k=int(len(nz)), supp=int((y > TOL).sum()), zratio=float(nz.max() / nz.min()) if len(nz) else 0.0,
               keep_frac=float(yz[(nd == 2) & (yz > TOL)].sum() / max(y.sum(), 1e-300)))
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--n', type=int)
    a = ap.parse_args()
    for q in a.q:
        with open(os.path.join(HERE, f'fam_q{q:02d}.jsonl'), 'a') as f:
            for idx in range(a.n): f.write(json.dumps(run(q, idx)) + '\n')
        print('done', q, flush=True)
