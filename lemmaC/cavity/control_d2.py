"""Control: the same n-fold relaxation population dynamics for ordinary bipartite matching
on the PWIT (team size 2, edge costs PPP of intensity 1).  There the LP is integral on every
bipartite graph and the integer limit is known: cost per vertex-pair pi^2/6 = 1.644934.
Recursion: G_k(x) = xi_k x + h_{j_k}(x);  energy a = int xi E[x*(xi)] dxi = E mean_c T_c^2/2,
T_c = -(s_a+s_b)_c, normalisation E mean_c T_c = 1."""
import numpy as np, sys, json
def update(pop, B, rng, XI):
    N, n = pop.shape
    K = rng.poisson(XI, size=B); Kmax = K.max()
    xi = XI * rng.random((B, Kmax))
    mask = np.arange(Kmax)[None, :] < K[:, None]
    j = rng.integers(0, N, (B, Kmax))
    cells = np.where(mask[:, :, None], xi[:, :, None] + pop[j], np.inf).reshape(B, Kmax * n)
    part = np.partition(cells, n - 1, axis=1)[:, :n]; part.sort(axis=1)
    return -part[:, ::-1]
def energy(pop, rng, M=200000):
    N = len(pop)
    T = np.clip(-(pop[rng.integers(0, N, M)] + pop[rng.integers(0, N, M)]), 0, None)
    return (T**2 / 2).mean(), T.mean()
n, init, N, sweeps, burn = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
rng = np.random.default_rng(7); XI = 14.0; B = 5000
X = rng.logistic(0, 1, N)
pop = np.repeat(-X[:, None], n, axis=1)
if init == 'kinked':
    pert = np.sort(rng.normal(0, 1.0, (N, n)), axis=1); pop += pert - pert.mean(axis=1, keepdims=True)
rec = []
for sw in range(sweeps):
    for s in range(0, N, B):
        pop[rng.integers(0, N, B)] = update(pop, B, rng, XI)
    if sw >= burn:
        a, m = energy(pop, rng)
        rec.append(dict(sweep=sw, a=a, norm=m, kink=float(np.mean(pop[:, -1] - pop[:, 0]))))
        print(json.dumps(rec[-1]), flush=True)
json.dump(rec, open(f'out/control_d2_n{n}_{init}.json', 'w'))
