"""Check of the 'global' route: compare sum of three bipartite optima S with LP* and T
(stored q=24..48 records; instances regenerated with gate1 seeds)."""
import json, sys, os, numpy as np
from scipy.optimize import linear_sum_assignment
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gate1'))
from gate1 import instance
def bip(M):
    r, c = linear_sum_assignment(M, maximize=True); return M[r, c].sum()
out = {}
for q in [24, 28, 32, 40, 48]:
    f = f'../gate1_ext/out/triangle_unif_q{q:02d}.jsonl'
    recs = [json.loads(l) for l in open(f)][:300]
    rows = []
    for r in recs:
        A, B, C = instance('triangle', 'unif', q, r['idx'])
        S = bip(A) + bip(B) + bip(C)
        rows.append((3*q - S, 3*q - r['V_lp'], (3*q - r['T']) if r.get('T') is not None else np.nan, r['V_lp'] - r['T'] if r.get('T') is not None else np.nan))
    a = np.array(rows)
    out[q] = dict(n=len(rows), deficit_bip_mean=float(a[:,0].mean()), deficit_LP_mean=float(a[:,1].mean()),
                  deficit_T_mean=float(np.nanmean(a[:,2])) if np.isfinite(a[:,2]).any() else None,
                  S_minus_LP_mean=float((a[:,1]-a[:,0]).mean()),
                  LP_minus_T_mean=float(np.nanmean(a[:,3])) if np.isfinite(a[:,3]).any() else None)
    print(q, out[q], flush=True)
json.dump(out, open('out/global_check.json', 'w'), indent=1)
