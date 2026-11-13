"""Ordinal ("lexicographic", infinitely-spread) triangle market.  The 3q^2 pair values are ranked by a
uniform random order and each dominates the sum of all lower ones.  Then the LP optimum is the lex-max
vector of pair weights, computed by sequential LPs: for pairs in decreasing rank, maximise the pair's
weight subject to unit loads and all earlier pair weights fixed.  Core nonempty iff every sequential
value is 0 or 1 (0/1 XY- and XZ-marginals force an integral lambda).  Stops at the first fractional
value.  Floating point (HiGHS), tolerance 1e-7; evidence only."""
import json, os, argparse
import numpy as np
from scipy.optimize import linprog
from scipy import sparse as sp
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))

def run(args):
    q, idx = args
    rng = np.random.default_rng(np.random.SeedSequence([20261111, 77, q, idx]))
    order = rng.permutation(3 * q * q)          # pair id: t*q*q + i*q + j, t in {XY, XZ, YZ}
    a, b, c = np.unravel_index(np.arange(q ** 3), (q, q, q))   # team (x=a, y=b, z=c)
    rows = np.repeat(np.arange(q ** 3), 3)
    loads = sp.csr_matrix((np.ones(3 * q ** 3), (rows, np.stack([a, q + b, 2 * q + c], 1).ravel())),
                          shape=(q ** 3, 3 * q)).T.tocsr()
    pair_of = [a * q + b, a * q + c, b * q + c]    # per type, pair index of each team
    def prow(pid):
        t, ij = divmod(pid, q * q)
        r = np.zeros(q ** 3); r[pair_of[t] == ij] = 1.0; return r
    fixed_rows, fixed_vals = [], []
    one = [np.zeros((q, q), bool) for _ in range(3)]  # pairs fixed at weight 1, per type
    for k, pid in enumerate(order):
        t, ij = divmod(pid, q * q); i, j = divmod(ij, q)
        if one[t][i].any() or one[t][:, j].any():     # an endpoint already has a full partner of this type
            val = 0.0
        else:
            r = prow(pid)
            Aeq = sp.vstack([loads] + ([sp.csr_matrix(np.array(fixed_rows))] if fixed_rows else [])).tocsr()
            beq = np.concatenate([np.ones(3 * q), np.array(fixed_vals)])
            res = linprog(-r, A_eq=Aeq, b_eq=beq, bounds=(0, None), method='highs-ds')
            assert res.success, res.message
            val = float(r @ res.x)
        if min(val, 1 - val) > 1e-7:
            return dict(q=q, idx=idx, empty=True, k=k, frac=k / (3 * q * q), val=val,
                        ones=int(sum(o.sum() for o in one)))
        val = round(val)
        if val == 1:
            one[t][i, j] = True
        if 1 or val == 1:
            fixed_rows.append(prow(pid)); fixed_vals.append(float(val))
        if sum(o.sum() for o in one) == 3 * q:
            return dict(q=q, idx=idx, empty=False, k=k, ones=3 * q)
    return dict(q=q, idx=idx, empty=False, k=len(order), ones=int(sum(o.sum() for o in one)))

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--n', type=int)
    a = ap.parse_args()
    for q in a.q:
        with ProcessPoolExecutor(4) as ex, open(os.path.join(HERE, 'out', f'lex_q{q:02d}.jsonl'), 'w') as f:
            rs = list(ex.map(run, [(q, i) for i in range(a.n)]))
            for r in rs: f.write(json.dumps(r) + '\n')
        print(q, 'P(empty)=%.3f' % np.mean([r['empty'] for r in rs]),
              'median k/3q^2 of first frac=%s' % (np.median([r['frac'] for r in rs if r['empty']]) if any(r['empty'] for r in rs) else None),
              'vals', sorted(set(round(r['val'], 3) for r in rs if r['empty'])), flush=True)
