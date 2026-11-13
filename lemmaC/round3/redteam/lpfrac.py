"""LP-only scaling: fraction of agents fractionally covered by the (basic) LP optimum, LP constant.
Reuses solve_lp/draw from ../../scaling/cg.py. Floating point; evidence only. Own seeds (tag=91)."""
import os, sys, json, time
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scaling'))
from cg import solve_lp, draw, teams_cols
def one(q, idx):
    t0 = time.time()
    A, B, C = draw(q, idx, 91)
    V = A[:, :, None] + B[:, None, :] + C[None, :, :]
    LPs, p, r, cols, x, it = solve_lp(V, q)
    fr = (x > 1e-7) & (x < 1 - 1e-7)
    tc = teams_cols(cols[fr], q)
    nfa = len(np.unique(tc))            # agents fractionally covered
    nint = int((x > 1 - 1e-7).sum())
    xs = x[fr]
    return dict(q=q, idx=idx, LPdef=3 * q - LPs, ffrac=nfa / (3 * q), nint=nint, nfrac_teams=int(fr.sum()),
                xmin=float(xs.min()) if len(xs) else 1.0, half=float(np.mean(np.abs(xs - 0.5) < 1e-6)) if len(xs) else 0.0,
                secs=time.time() - t0)
if __name__ == '__main__':
    out = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lpfrac.jsonl'), 'a')
    for q, n in [(16, 40), (24, 40), (32, 40), (48, 30), (64, 30), (96, 20), (128, 12), (160, 8), (200, 5)]:
        for i in range(n):
            r = one(q, i); out.write(json.dumps(r) + '\n'); out.flush()
        print('done', q, flush=True)
