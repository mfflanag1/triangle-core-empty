"""Energy estimators from a saved population of cavity messages (n cell slopes each).
a_direct = E mean_c T_c^4/8 (sensitive to population drift: d a/d shift = 9);
a_stat   = 3 E[-mean(s_fresh)] + E mean_c T_c^4/12, with s_fresh = one recursion step applied to
           the same population (Bethe form; stationary under a uniform shift of the messages,
           so the leading finite-population drift cancels).  Also: norm, fractional diagnostics."""
import numpy as np, sys, glob
from lp_popdyn import update
def est(pop, seed=0, M=400000, F=60000, XI=10.0):
    rng = np.random.default_rng(seed); N, n = pop.shape
    T = np.clip(-(pop[rng.integers(0, N, M)] + pop[rng.integers(0, N, M)] + pop[rng.integers(0, N, M)]), 0, None)
    m4 = (T**4).mean(axis=1); m3 = (T**3 / 6).mean(axis=1)
    fresh = np.vstack([update(pop, 5000, rng, XI) for _ in range(F // 5000)])
    node = -fresh.mean(axis=1)
    a_dir = m4.mean() / 8; a_stat = 3 * node.mean() + m4.mean() / 12
    se_stat = np.sqrt(9 * node.var() / len(node) + m4.var() / 144 / M)
    # expected number of teams at an agent with 0 < x* < 1 (fractional) and with x* = 1
    Ts = np.sort(T, axis=1)
    frac_teams = (Ts[:, -1]**3 - Ts[:, 0]**3).mean() / 6
    kinked_agents = np.mean(fresh[:, -1] - fresh[:, 0] > 1e-9)
    return dict(n=n, N=N, a_direct=a_dir, a_stat=a_stat, se_stat=se_stat, norm=m3.mean(),
                frac_teams_per_agent=frac_teams, kink=float(np.mean(pop[:, -1] - pop[:, 0])))
if __name__ == '__main__':
    for f in sys.argv[1:]:
        r = est(np.load(f)); print(f, ' '.join(f'{k}={v:.4f}' if isinstance(v, float) else f'{k}={v}' for k, v in r.items()), flush=True)
