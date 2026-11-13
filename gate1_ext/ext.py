"""Post-gate extension (2026-11-11), NOT part of the frozen gate 1.

Core-existence share only, triangle, U[0,1], larger q.  Same instances/seeds as gate1.py
(SeedSequence([20261111, 0, 0, q, idx])) and the same core test: dual-simplex LP vertex;
if integral -> core nonempty (a.s.-unique optimum); if fractional -> exact integer T by
scipy milp and core nonempty iff (LP - T)/T <= 1e-7.  No ratio/subsidy work.
Floating point (HiGHS): evidence, not certificates.
"""
import argparse, json, os, sys, time
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gate1'))
from gate1 import instance, structure, ok, INT_TOL, GAP_TOL


def solve(q, idx, milp_time=600.0, use_milp=True):
    t0 = time.time()
    v, M = structure(*instance('triangle', 'unif', q, idx))
    n = 3 * q
    lp = ok(linprog(-v, A_ub=M.T.tocsr(), b_ub=np.ones(n), bounds=(0, None), method='highs-ds'), 'primal LP')
    lam, V = lp.x, -lp.fun
    fmax = float(np.max(np.minimum(lam, 1 - lam)))
    rec = dict(q=q, idx=idx, V_lp=V, lp_integral=fmax < INT_TOL, lp_frac_max=fmax)
    if rec['lp_integral']:
        rec.update(T=float(v[lam > 0.5].sum()), rel_gap=0.0, core_nonempty=True, core_basis='integral_vertex')
    elif not use_milp:  # LP-only mode: a.s.-unique LP optimum is fractional => core empty (no gap computed)
        rec.update(T=None, core_nonempty=False, core_basis='fractional_vertex_only_lp_mode')
    else:
        mi = milp(-v, constraints=LinearConstraint(M.T.tocsr(), -np.inf, np.ones(n)),
                  integrality=np.ones(q ** 3), bounds=Bounds(0, 1),
                  options=dict(time_limit=milp_time, mip_rel_gap=1e-12))
        rec['milp_status'] = int(mi.status)
        if mi.status == 0:
            T = -mi.fun
            rec.update(T=T, rel_gap=(V - T) / T, core_nonempty=bool((V - T) / T <= GAP_TOL), core_basis='milp')
            if rec['core_nonempty']:
                rec['flag'] = 'fractional_vertex_but_zero_gap'
        else:  # a.s.-unique LP optimum is fractional => core empty
            rec.update(T=None, core_nonempty=False, core_basis='fractional_vertex_only')
    rec['secs'] = time.time() - t0
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--q', type=int, required=True)
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--lp-only', action='store_true', help='skip the MILP (no gap); core test = LP vertex integrality')
    a = ap.parse_args()
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out', f'triangle_unif_q{a.q:02d}.jsonl')
    done = set()
    if os.path.exists(path):
        with open(path) as f:
            done = {json.loads(l)['idx'] for l in f}
    todo = [i for i in range(a.n) if i not in done]
    t0 = time.time()
    with open(path, 'a') as f, ProcessPoolExecutor(a.workers) as ex:
        for rec in ex.map(solve, [a.q] * len(todo), todo, [600.0] * len(todo), [not a.lp_only] * len(todo), chunksize=1):
            f.write(json.dumps(rec) + '\n'); f.flush()
    print(f'{path}: {len(todo)} new in {time.time() - t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()
