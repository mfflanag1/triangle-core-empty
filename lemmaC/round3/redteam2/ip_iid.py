"""Integer optimum in the iid-team (Poissonised) model on the K-truncated hypergraph (upper bound on true IP cost;
exact for the truncated instance when status says optimal). Evidence only."""
import numpy as np, sys, time, json
from scipy.sparse import csc_matrix
from scipy.optimize import milp, LinearConstraint, Bounds
from lp_big import gen_iid, solve
q, K, nrep, seed, tl = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
for i in range(nrep):
    rng = np.random.default_rng([seed, q, i]); ta, tb, tc, xi = gen_iid(q, K, rng); n = len(xi)
    lp = solve(q, ta, tb, tc, xi)
    A = csc_matrix((np.ones(3 * n), (np.concatenate([ta, q + tb, 2 * q + tc]), np.tile(np.arange(n), 3))), shape=(3 * q, n))
    t0 = time.time()
    r = milp(xi, constraints=LinearConstraint(A, 1, 1), integrality=np.ones(n), bounds=Bounds(0, 1), options=dict(time_limit=tl))
    print(json.dumps(dict(q=q, K=K, i=i, lp=None if lp is None else round(lp['val'], 4), ip=None if r.x is None else round(float(r.fun / q), 4),
          bound=None if getattr(r, 'mip_dual_bound', None) is None else round(float(r.mip_dual_bound / q), 4), status=int(r.status), s=round(time.time() - t0, 1))), flush=True)
