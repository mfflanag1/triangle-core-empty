"""Population dynamics for the cavity recursion of the n-fold relaxation of the
local-limit perfect 3-matching problem (team weights x in {0, 1/n, ..., 1}, every agent
covered with total weight 1).  n = 1 is the integer problem; n -> infinity is the LP.

A cavity message h_j(x) = cost of agent j's subtree when the parent team takes weight x,
minus the cost at x = 0.  h_j is (discretely) convex with h_j(0) = 0, stored as its n
cell slopes s_j[0] <= ... <= s_j[n-1].  Integral RS fixed point: h_j(x) = -X_j x.

Recursion at agent i with children teams k (costs xi_k, PPP intensity xi^2/2, other agents
j_k, l_k):   G_k(x) = xi_k x + h_{j_k}(x) + h_{l_k}(x)   (cell slopes xi_k + s_j + s_l),
Phi = inf-convolution of the G_k (= the n smallest cell slopes over all k, sorted),
h_i(t) = Phi(1 - t) - Phi(1)  =>  s_i = - reversed(Phi slopes).

Edge marginal: team with cost xi and incoming messages h_a, h_b, h_c takes
x* = argmin_x xi x + h_a + h_b + h_c  = (# cells with xi + s_a + s_b + s_c < 0) / n.
Energy per team  a = int (xi^2/2) xi E[x*(xi)] dxi,  normalisation int (xi^2/2) E[x*] = 1.
"""
import numpy as np, sys, json

def update(pop, B, rng, XI):
    N, n = pop.shape
    lam = XI**3 / 6.0
    K = rng.poisson(lam, size=B)
    Kmax = K.max()
    # PPP on [0, XI] with density xi^2/2: xi = XI * U^(1/3)
    xi = XI * rng.random((B, Kmax)) ** (1 / 3)
    mask = np.arange(Kmax)[None, :] < K[:, None]
    j = rng.integers(0, N, (B, Kmax)); l = rng.integers(0, N, (B, Kmax))
    cells = xi[:, :, None] + pop[j] + pop[l]          # B x Kmax x n, sorted along last axis
    cells = np.where(mask[:, :, None], cells, np.inf).reshape(B, Kmax * n)
    part = np.partition(cells, n - 1, axis=1)[:, :n]
    part.sort(axis=1)
    return -part[:, ::-1]

def energy(pop, rng, M=400000, XI=None):
    """Exact integral over the team cost xi given the three incoming messages:
    per cell c with T_c = -(s_a+s_b+s_c)_c > 0, int_0^T xi^3/2 = T^4/8, int_0^T xi^2/2 = T^3/6."""
    N, n = pop.shape
    a, b, c = (pop[rng.integers(0, N, M)] for _ in range(3))
    T = np.clip(-(a + b + c), 0, None)
    e = (T**4 / 8).mean(axis=1); m = (T**3 / 6).mean(axis=1)
    return e.mean(), m.mean(), e.std() / np.sqrt(M)

def run(n, N=20000, sweeps=60, init='kinked', seed=0, XI=9.0, B=4000, log=None):
    rng = np.random.default_rng(seed)
    X = rng.normal(0.25, 0.86, N)
    if init == 'linear':
        pop = np.repeat(-X[:, None], n, axis=1)
    else:   # random convex perturbation of the linear fixed point ('tiny': scale 1e-4)
        pert = np.sort(rng.normal(0, 1e-4 if init == 'tiny' else 1.0, (N, n)), axis=1)
        pop = -X[:, None] + pert - pert.mean(axis=1, keepdims=True)
    hist = []
    for sw in range(sweeps):
        for s in range(0, N, B):
            idx = rng.integers(0, N, B)
            pop[idx] = update(pop, B, rng, XI)
        kink = np.mean(pop[:, -1] - pop[:, 0])
        frac_kinked = np.mean(pop[:, -1] - pop[:, 0] > 1e-9)
        if sw % 5 == 4 or sw == sweeps - 1 or init == 'tiny':
            a, nm, se = energy(pop, rng, XI=XI + 3)
            hist.append(dict(sweep=sw + 1, a=a, norm=nm, se=se, kink=kink, frac_kinked=frac_kinked,
                             mean_slope=float(pop.mean())))
            if log: print(n, init, json.dumps(hist[-1]), file=log, flush=True)
    return pop, hist

if __name__ == '__main__':
    n = int(sys.argv[1]); init = sys.argv[2]; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    N = int(sys.argv[4]) if len(sys.argv) > 4 else 20000
    sweeps = int(sys.argv[5]) if len(sys.argv) > 5 else 60
    pop, hist = run(n, N=N, sweeps=sweeps, init=init, seed=seed, log=sys.stdout)
    np.save(f'out/pop_n{n}_{init}_s{seed}.npy', pop)
