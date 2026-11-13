"""Check of the LP price-chain step and max TEAM deficits (M and LP support). Floating point; evidence only. Same seeds as pd.py."""
import os, sys, json
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling'))
from cg import draw, solve_lp
from pd import int_teams, TAG
def run(q, idx):
    A, B, C = draw(q, idx, TAG); dA, dB, dC = 1 - A, 1 - B, 1 - C     # cg.py: V[a,b,c]=A[a,b]+B[a,c]+C[b,c]
    V = A[:, :, None] + B[:, None, :] + C[None, :, :]; vf = V.reshape(-1)
    LPs, p, r, cols, x, it = solve_lp(V, q)
    pa, pb, pc = p[:q], p[q:2 * q], p[2 * q:]
    if np.max(np.minimum(x, 1 - x)) < 1e-6: Mi = cols[x > 0.5]
    else: Mi, _ = int_teams(vf, q, LPs, r)
    sup = cols[x > 1e-9]
    sa, sb, sc = np.unravel_index(sup, (q, q, q))
    # initial step of the chain, sector b: p_b' - p_b >= d_ab + d_bc - d(a,b') - d(b',c) for every support team (a,b,c), every b'
    viol = 0.0; under = []
    for a, b, c in zip(sa, sb, sc):
        lhs = pb - pb[b] - (dA[a, b] + dC[b, c]) + dA[a, :] + dC[:, c]
        viol = min(viol, lhs.min())
    # how underpriced is each agent of a support team relative to the median price of its sector, in units of the team deficit
    dteam = 3 - vf
    spr = max(float(z.max() - np.quantile(z, 0.1)) for z in (pa, pb, pc)); rng = max(float(z.max() - z.min()) for z in (pa, pb, pc))
    return dict(q=q, idx=idx, maxteamM=float(dteam[Mi].max()), maxteamLP=float(dteam[sup].max()), chain_viol=float(viol), spr=spr, rng=rng,
                tight=float(np.abs(r[sup]).max()))
if __name__ == '__main__':
    for q, n in ((16, 100), (24, 100), (32, 100), (48, 20)):
        R = [run(q, i) for i in range(n)]
        g = lambda k: np.array([r_[k] for r_ in R])
        print(f"q={q} n={n} max team deficit M: mean {g('maxteamM').mean()*q**(2/3):.2f} max {g('maxteamM').max()*q**(2/3):.2f} (x q^2/3) | "
              f"LP support: mean {g('maxteamLP').mean()*q**(2/3):.2f} max {g('maxteamLP').max()*q**(2/3):.2f} | x q^1/2: M {g('maxteamM').mean()*q**.5:.2f} LP {g('maxteamLP').mean()*q**.5:.2f} | "
              f"price max-q10 {g('spr').mean()*q**(2/3):.2f}, max-min {g('rng').mean()*q**(2/3):.2f} (x q^2/3) | min chain-step slack {g('chain_viol').min():.2e} | max |reduced cost| on support {g('tight').max():.1e}", flush=True)
