"""Used pair deficits of the integer optimum M and of the LP optimum (pair loads). Floating point; evidence only.
Reuses solve_lp / amat from ../../scaling/cg.py. Own seeds: tag 31."""
import os, sys, json, time
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling'))
from cg import draw, solve_lp, amat
TAG = 31

def int_teams(vf, q, LPs, r):
    n = 3 * q; delta = 0.15
    while True:
        cand = np.nonzero(r >= -delta)[0]
        mi = milp(-vf[cand], constraints=LinearConstraint(amat(cand, q), np.ones(n), np.ones(n)),
                  integrality=np.ones(len(cand)), bounds=Bounds(0, 1), options=dict(time_limit=900, mip_rel_gap=1e-12))
        if mi.status == 0 and -mi.fun >= LPs - delta + 1e-9:
            return cand[mi.x > 0.5], -mi.fun
        delta *= 2

def run(args):
    q, idx = args
    t0 = time.time()
    A, B, C = draw(q, idx, TAG)
    dA, dB, dC = 1 - A, 1 - B, 1 - C
    V = A[:, :, None] + B[:, None, :] + C[None, :, :]
    vf = V.reshape(-1)
    LPs, p, r, cols, x, it = solve_lp(V, q)
    if np.max(np.minimum(x, 1 - x)) < 1e-6:
        Mi = cols[x > 0.5]; T = float(vf[Mi].sum())
    else:
        Mi, T = int_teams(vf, q, LPs, r)
    a, b, c = np.unravel_index(Mi, (q, q, q))
    dM = np.concatenate([dA[a, b], dB[a, c], dC[b, c]])  # cg.py convention: V[a,b,c]=A[a,b]+B[a,c]+C[b,c]
    # LP pair loads
    sup = x > 1e-9
    sa, sb, sc = np.unravel_index(cols[sup], (q, q, q)); w = x[sup]
    LA = np.zeros((q, q)); LB = np.zeros((q, q)); LC = np.zeros((q, q))
    np.add.at(LA, (sa, sb), w); np.add.at(LB, (sa, sc), w); np.add.at(LC, (sb, sc), w)
    load = np.concatenate([LA.ravel(), LB.ravel(), LC.ravel()])
    dall = np.concatenate([dA.ravel(), dB.ravel(), dC.ravel()])
    inM = np.zeros(3 * q * q)
    inM[np.concatenate([a * q + b, q * q + a * q + c, 2 * q * q + b * q + c])] = 1
    pos = load > 1e-7
    # 2-swap repair cost for every M pair of type A (pair ab): min over other teams j and both swap types of the
    # sum of the four NEW pair deficits (rigorous upper bound for the used deficit d_ab by optimality of M)
    R2 = []
    for t in range(q):
        at, bt, ct = a[t], b[t], c[t]
        o = np.arange(q) != t
        sb_ = dA[at, b[o]] + dC[b[o], ct] + dA[a[o], bt] + dC[bt, c[o]]      # swap the b's
        sa_ = dA[a[o], bt] + dB[a[o], ct] + dA[at, b[o]] + dB[at, c[o]]      # swap the a's
        R2.append(float(min(sb_.min(), sa_.min())))
    s23 = q ** (2 / 3)
    rec = dict(q=q, idx=idx, LP=LPs, T=T, D=LPs - T, defT=3 * q - T, defLP=3 * q - LPs,
               maxM=float(dM.max()), sum2M=float((dM ** 2).sum()),
               maxLP=float(dall[pos].max()), sum2LP=float((load * dall ** 2).sum()), sum22LP=float((load ** 2 * dall ** 2).sum()),
               gradD2=float((((load - inM) * dall) ** 2).sum()),
               nM=[int((dM * s23 > K).sum()) for K in (2, 3, 4, 6)],
               wLP=[float(load[dall * s23 > K].sum()) for K in (2, 3, 4, 6)],
               nsupLP=int(pos.sum()), maxR2=max(R2), maxdA=float(dA[a, b].max()),
               slackR2=float(min(np.array(R2) - dA[a, b])), secs=time.time() - t0)
    return rec

if __name__ == '__main__':
    q = int(sys.argv[1]); n = int(sys.argv[2]); workers = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    fn = os.path.join(HERE, f'pd_q{q:03d}.jsonl')
    with ProcessPoolExecutor(workers) as ex, open(fn, 'a') as f:
        for rr in ex.map(run, [(q, i) for i in range(n)]):
            f.write(json.dumps(rr) + '\n'); f.flush()
