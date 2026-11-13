"""Superadditivity probe. For a partition of M into blocks S_1..S_m (Lemma 1(i): M|S_j is the integer
optimum of G_{S_j}), gap(G) >= sum_j gap(G_{S_j}) with gap(G_S) = LP*(G_S) - v(M|S).
For random partitions of M into blocks of size s, record the share of blocks with positive sub-gap and
sum_j gap(G_Sj) / gap(G). Uses M from out/struct_q*.jsonl runs (recomputed). Floating point; evidence only."""
import json, sys, os, argparse
import numpy as np
from scipy.optimize import linprog
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'gate1')); sys.path.insert(0, HERE)
from gate1 import instance, structure, ok
from mustruct import integer_opt

def run(args):
    q, idx, sizes, reps = args
    A, B, C = instance('triangle', 'unif', q, idx)
    v, Mm = structure(A, B, C); n = 3 * q
    lp = ok(linprog(-v, A_eq=Mm.T.tocsr(), b_eq=np.ones(n), bounds=(0, None), method='highs-ds'), 'LP')
    lam, LPs = lp.x, -lp.fun; rlp = v - Mm @ (-lp.eqlin.marginals)
    if np.max(np.minimum(lam, 1 - lam)) < 1e-6: sel = lam > 0.5; T = LPs
    else: sel, T, _, _ = integer_opt(v, Mm, n, LPs, rlp)
    tm = np.array(np.unravel_index(np.nonzero(sel)[0], (q, q, q))).T
    rng = np.random.default_rng([11, q, idx]); out = dict(q=q, idx=idx, gap=LPs - T)
    for s in sizes:
        pos = 0; tot = 0; sums = []
        for _ in range(reps):
            perm = rng.permutation(q); S_sum = 0.0
            for j in range(q // s):
                S = perm[j * s:(j + 1) * s]; X, Y, Z = tm[S, 0], tm[S, 1], tm[S, 2]
                vs, Ms = structure(A[np.ix_(X, Y)], B[np.ix_(X, Z)], C[np.ix_(Y, Z)])
                r = ok(linprog(-vs, A_eq=Ms.T.tocsr(), b_eq=np.ones(3 * s), bounds=(0, None), method='highs-ds'), 'sub')
                g = -r.fun - (A[X, Y] + B[X, Z] + C[Y, Z]).sum()
                tot += 1; pos += g > 1e-9; S_sum += max(g, 0)
            sums.append(S_sum)
        out[f's{s}'] = dict(blocks=tot, pos=int(pos), sum_mean=float(np.mean(sums)), any_pos=float(np.mean(np.array(sums) > 1e-9)))
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--n', type=int)
    a = ap.parse_args()
    for q in a.q:
        sizes = sorted({s for s in (4, 6, 8, 12, q // 4, q // 2) if 4 <= s <= q // 2})
        with ProcessPoolExecutor(4) as ex, open(os.path.join(HERE, 'out', f'blocks_q{q:02d}.jsonl'), 'w') as f:
            for r in ex.map(run, [(q, i, sizes, 3) for i in range(a.n)]):
                f.write(json.dumps(r) + '\n'); f.flush()
        print('done', q, flush=True)
