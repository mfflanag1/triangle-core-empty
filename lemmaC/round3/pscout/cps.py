import numpy as np, sys, time
from ortools.sat.python import cp_model
from sim import gen
SC = 10**6
def cpsat(q, K, seed, tlim, workers=1, hint=None):
    a, b, c, xi = gen(q, K, seed)
    n = len(xi)
    m = cp_model.CpModel()
    x = [m.NewBoolVar('') for _ in range(n)]
    for arr in (a, b, c):
        order = np.argsort(arr, kind='stable'); bounds = np.searchsorted(arr[order], np.arange(q + 1))
        for v in range(q):
            m.AddExactlyOne([x[j] for j in order[bounds[v]:bounds[v + 1]]])
    w = np.rint(xi * SC).astype(np.int64)
    m.Minimize(sum(int(w[j]) * x[j] for j in range(n)))
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = tlim; s.parameters.num_workers = workers
    t0 = time.time(); st = s.Solve(m)
    ok = st in (cp_model.OPTIMAL, cp_model.FEASIBLE)
    return dict(ub=(s.ObjectiveValue() / SC / q if ok else None), lb=s.BestObjectiveBound() / SC / q,
                closed=int(st == cp_model.OPTIMAL), tip=round(time.time() - t0, 1))
if __name__ == '__main__':
    q, seed, tlim, wk = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
    print(q, seed, wk, cpsat(q, 6.0, seed, tlim, wk), flush=True)
