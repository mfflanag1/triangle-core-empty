"""Long population-dynamics runs with per-sweep energy, for error bars (batch means)."""
import numpy as np, sys, json
from lp_popdyn import update, energy
n, init, seed, N, sweeps, burn = (int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4]),
                                  int(sys.argv[5]), int(sys.argv[6]))
XI = 10.0; B = 5000
rng = np.random.default_rng(seed)
X = rng.normal(0.25, 0.86, N)
if init == 'linear':
    pop = np.repeat(-X[:, None], n, axis=1)
else:
    pert = np.sort(rng.normal(0, 1.0, (N, n)), axis=1)
    pop = -X[:, None] + pert - pert.mean(axis=1, keepdims=True)
rec = []
for sw in range(sweeps):
    for s in range(0, N, B):
        idx = rng.permutation(N)[:B] if False else rng.integers(0, N, B)
        pop[idx] = update(pop, B, rng, XI)
    if sw >= burn:
        a, nm, se = energy(pop, rng, M=200000)
        # fraction of fractional team marginals among teams with x*>0, and fractional root agents
        M2 = 100000
        aa, bb, cc = (pop[rng.integers(0, N, M2)] for _ in range(3))
        T = np.clip(-(aa + bb + cc), 0, None)
        # P(x* fractional | team in support), weighting by intensity: x* = k/n on xi in (T_{k+1}, T_k)
        Ts = np.sort(T, axis=1)[:, ::-1]  # decreasing; x* >= (k+1)/n iff xi < Ts[k]
        supp = Ts[:, 0]**3 / 6; full = Ts[:, -1]**3 / 6
        rec.append(dict(sweep=sw, a=a, norm=nm, kink=float(np.mean(pop[:, -1] - pop[:, 0])),
                        supp_mass=float(supp.mean()), full_mass=float(full.mean()), maxabs=float(np.abs(pop).max())))
        print(json.dumps(rec[-1]), flush=True)
json.dump(rec, open(f'out/precise_n{n}_{init}_s{seed}.json', 'w'))
np.save(f'out/precise_pop_n{n}_{init}_s{seed}.npy', pop)
