"""Route-3 locality probe. Lemma 1(ii) (THEOREM_VANISHING) gives, for ANY partition of M into blocks S_j,
D(G) >= sum_j D(G_Sj).  Random partitions capture nothing (structure/blocks.py).  Here blocks are grown
greedily along (a) the support of the LP optimum lam* ('lam' mode, an oracle choice: best case for any
local rule) or (b) the near-tight hypergraph {f notin M : v_f - p_lp(f) > -q^(-2/3)} ('tight' mode).
Reports captured share sum_j D(G_Sj) / D(G).  Floating point; evidence only."""
import json, sys, os, argparse
import numpy as np
from scipy.optimize import linprog
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'gate1')); sys.path.insert(0, os.path.join(HERE, '..', 'structure'))
from gate1 import instance, structure, ok
from mustruct import integer_opt

def blocks_greedy(W, s, rng):
    q = W.shape[0]; un = set(range(q)); out = []
    while un:
        seed = rng.choice(sorted(un)); B = [seed]; un.remove(seed)
        while len(B) < s and un:
            cand = sorted(un); sc = W[np.ix_(B, cand)].sum(0) + 1e-9 * rng.random(len(cand))
            j = cand[int(np.argmax(sc))]; B.append(j); un.remove(j)
        out.append(np.array(B))
    return out

def subgap(A, B, C, tm, S):
    X, Y, Z = tm[S, 0], tm[S, 1], tm[S, 2]
    vs, Ms = structure(A[np.ix_(X, Y)], B[np.ix_(X, Z)], C[np.ix_(Y, Z)])
    r = ok(linprog(-vs, A_eq=Ms.T.tocsr(), b_eq=np.ones(3 * len(S)), bounds=(0, None), method='highs-ds'), 'sub')
    return max(-r.fun - (A[X, Y] + B[X, Z] + C[Y, Z]).sum(), 0.0)

def run(args):
    q, idx, sizes = args
    A, B, C = instance('triangle', 'unif', q, idx)
    v, Mm = structure(A, B, C); n = 3 * q
    lp = ok(linprog(-v, A_eq=Mm.T.tocsr(), b_eq=np.ones(n), bounds=(0, None), method='highs-ds'), 'LP')
    lam, LPs = lp.x, -lp.fun; rlp = v - Mm @ (-lp.eqlin.marginals)
    if np.max(np.minimum(lam, 1 - lam)) < 1e-6: return dict(q=q, idx=idx, gap=0.0)
    sel, T, _, _ = integer_opt(v, Mm, n, LPs, rlp)
    tm = np.array(np.unravel_index(np.nonzero(sel)[0], (q, q, q))).T
    team_of = np.zeros((3, q), int)
    for k, (x, y, z) in enumerate(tm): team_of[0, x] = k; team_of[1, y] = k; team_of[2, z] = k
    out = dict(q=q, idx=idx, gap=LPs - T)
    for mode in ('lam', 'tight'):
        W = np.zeros((q, q))
        if mode == 'lam': js = np.nonzero((lam > 1e-7) & ~sel)[0]; wts = lam[js]
        else: js = np.nonzero((rlp > -q ** (-2 / 3)) & ~sel)[0]; wts = np.ones(len(js))
        for j, w in zip(js, wts):
            x, y, z = np.unravel_index(j, (q, q, q)); ks = {team_of[0, x], team_of[1, y], team_of[2, z]}
            for a in ks:
                for b in ks:
                    if a != b: W[a, b] += w
        rng = np.random.default_rng([7, q, idx])
        for s in sizes:
            caps = []
            for _ in range(2):
                caps.append(sum(subgap(A, B, C, tm, S) for S in blocks_greedy(W, s, rng)))
            out[f'{mode}_s{s}'] = float(np.mean(caps)) / (LPs - T)
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--n', type=int)
    a = ap.parse_args()
    for q in a.q:
        sizes = sorted({s for s in (6, 12, q // 4, q // 2) if 4 <= s <= q // 2})
        with ProcessPoolExecutor(4) as ex, open(os.path.join(HERE, 'out', f'local_q{q:02d}.jsonl'), 'w') as f:
            for r in ex.map(run, [(q, i, sizes) for i in range(a.n)]):
                f.write(json.dumps(r) + '\n'); f.flush()
        print('done', q, flush=True)
