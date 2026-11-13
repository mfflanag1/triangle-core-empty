"""Theorem A numerics: for every pair (a,b) used by M (team (a,b,c)), the cheapest b-rotation of M - {(a,b,c)} + (a,b,c)
(any length; cost = sum of NEW pair deficits) via Dijkstra. By optimality of M it must be >= d(a,b). Floating point; evidence only."""
import os, sys
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.sparse.csgraph import dijkstra
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling'))
from cg import draw, solve_lp
from pd import int_teams, TAG
def run(q, idx):
    A, B, C = draw(q, idx, TAG); dA, dB, dC = 1 - A, 1 - B, 1 - C     # cg.py: V[a,b,c]=A[a,b]+B[a,c]+C[b,c]
    V = A[:, :, None] + B[:, None, :] + C[None, :, :]; vf = V.reshape(-1)
    LPs, p, r, cols, x, it = solve_lp(V, q)
    if np.max(np.minimum(x, 1 - x)) < 1e-6: Mi = cols[x > 0.5]
    else: Mi, _ = int_teams(vf, q, LPs, r)
    a, b, c = np.unravel_index(Mi, (q, q, q))
    W = dA[a[None, :], b[:, None]] + dC[b[:, None], c[None, :]]      # W[j,l]: b_j joins position (a_l,c_l)
    Rinf, R2, hops, slack = [], [], [], []
    for t in range(q):
        o = np.nonzero(np.arange(q) != t)[0]; m = len(o)
        G = np.full((m + 2, m + 2), np.inf)
        G[:m, :m] = W[np.ix_(o, o)]; np.fill_diagonal(G, np.inf)
        G[m, :m] = W[t, o]             # b = b_t joins team j
        G[:m, m + 1] = W[o, t]         # b_l joins (a_t, c_t)
        dist, pred = dijkstra(G, indices=m, return_predecessors=True)
        Rinf.append(dist[m + 1]); R2.append((W[t, o] + W[o, t]).min())
        h = 0; u = m + 1
        while pred[u] >= 0: u = pred[u]; h += 1
        hops.append(h); slack.append(dist[m + 1] - dA[a[t], b[t]])
    return max(Rinf), max(R2), float(np.mean(hops)), min(slack), float(dA[a, b].max())
if __name__ == '__main__':
    for q, n in ((16, 60), (24, 60), (32, 60), (48, 15)):
        R = np.array([run(q, i) for i in range(n)])
        print(f"q={q} n={n} max over M's ab-pairs of cheapest rotation: mean {R[:,0].mean():.3f} (x q^1/2: {R[:,0].mean()*q**.5:.2f}; x q^2/3: {R[:,0].mean()*q**(2/3):.2f}) | "
              f"2-swap only: {R[:,1].mean():.3f} (x q^1/4: {R[:,1].mean()*q**.25:.2f}) | mean moves {R[:,2].mean():.2f} | min slack (rotation - d_ab) {R[:,3].min():.4f} | max used d_ab x q^2/3 {R[:,4].mean()*q**(2/3):.2f}", flush=True)
