"""Direct test of the deterministic part of Theorem B (LP price chains). Floating point; evidence only.
For each support team f=(a,b,c) of the LP optimum and threshold s: BFS the sets B_i (sector b; A_i, C_i symmetric)
using ALL support teams, check the price inequality p_z > p_f-agent + (two deficits of f at that agent) - (depth+1)s
for every reached agent z, and whenever the three reached sets each exceed 2q/3 check d(f) < 3(maxdepth+1)s."""
import os, sys
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling'))
from cg import draw, solve_lp
TOL = 1e-7

def bfs(cost0, costpos, sup_sector, s, q):
    """cost0[z] = cost for z to join f's position; costpos[g, z] = cost for z to join the position of support team g;
    sup_sector[g] = the agent of this sector in support team g. Returns depth per agent (-1 if unreached)."""
    depth = np.full(q, -1); cur = np.nonzero(cost0 < s)[0]; depth[cur] = 0; i = 0
    while len(cur):
        i += 1
        gs = np.nonzero(np.isin(sup_sector, cur))[0]
        new = np.nonzero((costpos[gs] < s).any(0) & (depth < 0))[0]
        depth[new] = i; cur = new
    return depth

def run(q, idx, svals):
    A, B, C = draw(q, idx, 31); dA, dB, dC = 1 - A, 1 - B, 1 - C      # cg.py: V[a,b,c]=A[a,b]+B[a,c]+C[b,c]
    V = A[:, :, None] + B[:, None, :] + C[None, :, :]
    LPs, p, r, cols, x, it = solve_lp(V, q)
    pa, pb, pc = p[:q], p[q:2 * q], p[2 * q:]
    sup = cols[x > 1e-9]; sa, sb, sc = np.unravel_index(sup, (q, q, q)); dsup = 3 - V.reshape(-1)[sup]
    posA = dA[:, sb].T + dB[:, sc].T      # [g, a'] : a' joins position (b_g, c_g)
    posB = dA[sa, :] + dC[:, sc].T        # [g, b'] : b' joins position (a_g, c_g)
    posC = dB[sa, :] + dC[sb, :]          # [g, c'] : c' joins position (a_g, b_g)
    out = dict(price_checks=0, price_viol=0, close_cases=0, close_viol=0, worst=np.inf)
    for s in svals:
        for a, b, c, df in zip(sa, sb, sc, dsup):
            dep = [bfs(dA[:, b] + dB[:, c], posA, sa, s, q), bfs(dA[a, :] + dC[:, c], posB, sb, s, q), bfs(dB[a, :] + dC[b, :], posC, sc, s, q)]
            base = [pa[a] + dA[a, b] + dB[a, c], pb[b] + dA[a, b] + dC[b, c], pc[c] + dB[a, c] + dC[b, c]]
            for P, d_, b0 in zip((pa, pb, pc), dep, base):
                m = d_ >= 0
                slack = P[m] - b0 + (d_[m] + 1) * s
                out['price_checks'] += int(m.sum()); out['price_viol'] += int((slack < -TOL).sum())
                if m.any(): out['worst'] = min(out['worst'], float(slack.min()))
            if all((d_ >= 0).sum() > 2 * q / 3 for d_ in dep):
                out['close_cases'] += 1
                rmax = max(int(d_.max()) for d_ in dep)
                out['close_viol'] += int(df >= 3 * (rmax + 1) * s + TOL)
    return out

if __name__ == '__main__':
    for q, n in ((16, 30), (24, 30), (32, 20)):
        tot = dict(price_checks=0, price_viol=0, close_cases=0, close_viol=0); worst = np.inf
        for i in range(n):
            o = run(q, i, [0.5 * q ** -0.5, 1.0 * q ** -0.5, 2.0 * q ** -0.5, 3.0 * q ** -0.5])
            for k in tot: tot[k] += o[k]
            worst = min(worst, o['worst'])
        print(q, n, tot, 'worst price slack', worst, flush=True)
