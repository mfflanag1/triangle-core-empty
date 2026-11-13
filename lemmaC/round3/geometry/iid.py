"""Control: independent team values with the same marginal law (deficit of each team = sum of three fresh U[0,1], independently
per team), so no pair is shared by two teams.  This is the model the pair-additive market resembles locally when q is large
(two cheap teams through an agent share a pair with probability O(q^(-1/3))).  D = LP* - T by ../../scaling/cg.py.  Evidence only.
Seeds: SeedSequence([20261111, 3, 5, q, idx])."""
import os, sys, json, time
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling'))
from cg import integer_opt, solve_lp
def run(q, idx):
    t0 = time.time()
    rng = np.random.default_rng(np.random.SeedSequence([20261111, 3, 5, q, idx]))
    V = 3 - rng.random((3, q, q, q)).sum(0); vf = V.reshape(-1)
    LPs, p, r, cols, x, it = solve_lp(V, q)
    if np.max(np.minimum(x, 1 - x)) < 1e-6: T = LPs
    else: T, _ = integer_opt(vf, q, LPs, r, 0.15 if q <= 48 else 0.3)
    return dict(q=q, idx=idx, D=LPs - T, defT=3 * q - T, secs=time.time() - t0)
if __name__ == '__main__':
    q = int(sys.argv[1]); nn = int(sys.argv[2])
    with open(os.path.join(HERE, f'iid_q{q:02d}.jsonl'), 'a') as f:
        for idx in range(nn): f.write(json.dumps(run(q, idx)) + '\n'); f.flush()
    R = [json.loads(l) for l in open(os.path.join(HERE, f'iid_q{q:02d}.jsonl'))]; D = np.array([r['D'] for r in R])
    print(f"iid q={q} n={len(R)} P(D>0)={(D > 1e-7).mean():.3f} ED={D.mean():.4f} sd={D.std():.4f} E(3q-T)={np.mean([r['defT'] for r in R]):.3f} = {np.mean([r['defT'] for r in R]) / q ** (1 / 3):.3f} q^(1/3) secs/inst={np.mean([r['secs'] for r in R]):.2f}", flush=True)
