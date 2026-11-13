"""V3 referee test of Theorem B's deterministic steps on LP optima, q=8..16. Floating point, evidence only.
Follows the proof literally: B_0=N({(a,c)}), S_i = max matching among (a,c)-positions of support teams with b in B_i,
check |S_i|>=|B_i|, B_{i+1}=N(S_i), price inequality at each depth; Close for any support team in A_i x B_j x C_k.
Also repeats with a second optimal dual (random objective on the optimal dual face)."""
import os, sys
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, coo_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
TOL = 1e-7

def build(q):
    a, b, c = np.unravel_index(np.arange(q ** 3), (q, q, q))
    rows = np.concatenate([a, q + b, 2 * q + c]); cols = np.tile(np.arange(q ** 3), 3)
    return a, b, c, csr_matrix((np.ones(3 * q ** 3), (rows, cols)), shape=(3 * q, q ** 3))

def chain(q, s, pos0_cost, D1, D2, sup, k1, k2, kz, P, base, x, stats):
    """sector Z growing. pos0_cost[z]: cost of z joining f's position. D1[u,z], D2[w,z]: deficits of z to the two
    other sectors. sup: support teams array (n,3) columns; k1,k2 = columns of the position, kz = column of sector Z.
    Returns list of sets by depth."""
    sets = []
    cur = np.nonzero(pos0_cost < s)[0]
    depth = 0
    while len(cur) and depth <= 6:
        sl = P[cur] - base + (depth + 1) * s
        stats['price_checks'] += len(cur); stats['price_viol'] += int((sl <= -TOL).sum())
        stats['worst'] = min(stats['worst'], float(sl.min()))
        sets.append(set(cur.tolist()))
        g = sup[np.isin(sup[:, kz], cur)]
        # mass carried
        # bipartite matching among positions (g[:,k1], g[:,k2])
        G = coo_matrix((np.ones(len(g)), (g[:, k1], g[:, k2])), shape=(q, q)).tocsr()
        G.data[:] = 1
        mt = maximum_bipartite_matching(G, perm_type='column')   # mt[row] = col or -1
        rowsm = np.nonzero(mt >= 0)[0]
        stats['growth_checks'] += 1
        if len(rowsm) < len(cur): stats['growth_viol'] += 1
        # N(S_i)
        costs = D1[rowsm, :] + D2[mt[rowsm], :]                  # [pos, z]
        cur = np.nonzero((costs < s).any(0))[0]
        depth += 1
    return sets

def run(q, seed, svals, stats, alt_dual=False):
    rng = np.random.default_rng(seed)
    dAB, dBC, dAC = rng.random((q, q)), rng.random((q, q)), rng.random((q, q))
    a, b, c, Aeq = build(q)
    d = dAB[a, b] + dBC[b, c] + dAC[a, c]
    res = linprog(d, A_eq=Aeq, b_eq=np.ones(3 * q), bounds=(0, None), method='highs')
    x = res.x; u = res.eqlin.marginals.copy()
    if alt_dual:
        # another optimal dual: max w.u  s.t. A^T u <= d, sum u = LP value
        w = rng.random(3 * q)
        w = w - 0.5
        from scipy.sparse import vstack
        Aub = vstack([Aeq.T.tocsr(), csr_matrix(-np.ones((1, 3 * q)))]).tocsr()
        r2 = linprog(-w, A_ub=Aub, b_ub=np.concatenate([d, [-(res.fun - 1e-10)]]), bounds=list(zip(u - 1, u + 1)), method='highs')
        stats['alt_status'] = stats.get('alt_status', 0) + int(r2.status != 0)
        if r2.status == 0:
            stats['alt_dual_dist'] = max(stats.get('alt_dual_dist', 0), float(np.abs(r2.x - u).max())); u = r2.x
    red = d - Aeq.T @ u
    assert red.min() > -1e-6, red.min()
    supi = np.nonzero(x > 1e-9)[0]
    assert np.abs(red[supi]).max() < 1e-6
    p = 1 - u; pa, pb, pc = p[:q], p[q:2 * q], p[2 * q:]
    sup = np.stack([a[supi], b[supi], c[supi]], 1)
    stats['frac'] += int(np.max(np.minimum(x, 1 - x)) > 1e-6); stats['inst'] += 1
    supset = set(map(tuple, sup.tolist()))
    for s in svals:
        for (fa, fb, fc) in sup:
            df = dAB[fa, fb] + dBC[fb, fc] + dAC[fa, fc]
            # sector 2 (b): position (a,c); z=b' cost dAB[a_,z]+dBC[z,c_]
            SB = chain(q, s, dAB[fa, :] + dBC[:, fc], dAB, dBC.T, sup, 0, 2, 1, pb, pb[fb] + dAB[fa, fb] + dBC[fb, fc], x, stats)
            # sector 1 (a): position (b,c); z=a' cost dAB[z,b_]+dAC[z,c_]
            SA = chain(q, s, dAB[:, fb] + dAC[:, fc], dAB.T, dAC.T, sup, 1, 2, 0, pa, pa[fa] + dAB[fa, fb] + dAC[fa, fc], x, stats)
            # sector 3 (c): position (a,b); z=c' cost dAC[a_,z]+dBC[b_,z]
            SC = chain(q, s, dAC[fa, :] + dBC[fb, :], dAC, dBC, sup, 0, 1, 2, pc, pc[fc] + dAC[fa, fc] + dBC[fb, fc], x, stats)
            # Close: any support team t in A_i x B_j x C_k  =>  d(f) < (i+j+k+3)s
            for (ta, tb, tc) in sup:
                ia = next((i for i, S in enumerate(SA) if ta in S), None)
                ib = next((i for i, S in enumerate(SB) if tb in S), None)
                ic = next((i for i, S in enumerate(SC) if tc in S), None)
                if ia is None or ib is None or ic is None: continue
                stats['close_checks'] += 1
                marg = (ia + ib + ic + 3) * s - df
                stats['close_worst'] = min(stats['close_worst'], marg)
                if marg <= -TOL: stats['close_viol'] += 1
            # mass claim: if all three final sets >= 0.9q there must be a support team inside the product
            for r in range(min(len(SA), len(SB), len(SC))):
                if min(len(SA[r]), len(SB[r]), len(SC[r])) >= 0.9 * q:
                    stats['mass_checks'] += 1
                    if not any(t[0] in SA[r] and t[1] in SB[r] and t[2] in SC[r] for t in supset): stats['mass_viol'] += 1

if __name__ == '__main__':
    for alt in (True,):
        for q, n in ((8, 40), (10, 30), (12, 20), (14, 12), (16, 8)):
            stats = dict(inst=0, frac=0, price_checks=0, price_viol=0, worst=np.inf, growth_checks=0, growth_viol=0,
                         close_checks=0, close_viol=0, close_worst=np.inf, mass_checks=0, mass_viol=0)
            for i in range(n):
                run(q, 77000 + 100 * q + i, [0.5 * q ** -0.5, 1.0 * q ** -0.5, 1.5 * q ** -0.5, 2.5 * q ** -0.5], stats, alt_dual=alt)
            print('alt_dual' if alt else 'highs_dual', q, n, stats, flush=True)
