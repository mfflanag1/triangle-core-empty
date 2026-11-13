"""Post-gate extension (2026-11-11): float search for a triangle market with rho > 1/3.
Small q, non-uniform laws (sparse 0/1, small integers, heavy tail, planted structure).
Reuses gate1.solve with gate1.instance monkeypatched.  Evidence only."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gate1'))
import gate1

LAWS = {
    'bern': lambda r, q: (r.random((q, q)) < 0.3).astype(float),
    'int3': lambda r, q: r.integers(0, 4, (q, q)).astype(float),
    'pareto': lambda r, q: r.pareto(1.0, (q, q)),
    'sq': lambda r, q: r.random((q, q)) ** 4,
}
best = (-1, None)
cnt = 0
for law, f in LAWS.items():
    for q in (2, 3, 4):
        for idx in range(150):
            r = np.random.default_rng([20261112, list(LAWS).index(law), q, idx])
            M = (f(r, q), f(r, q), f(r, q))
            gate1.instance = lambda *a, M=M: M
            try:
                rec = gate1.solve('triangle', 'unif', q, idx)
            except Exception as ex:
                continue
            if not rec.get('core_nonempty') or 'rho' not in rec:
                continue
            cnt += 1
            if rec['rho'] > best[0]:
                best = (rec['rho'], (law, q, idx))
print('nonempty instances examined', cnt, 'max rho', best)
