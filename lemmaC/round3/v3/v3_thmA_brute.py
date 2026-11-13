"""V3 referee test of Theorem A's deterministic steps by brute force (q=3,4,5). Floating point, evidence only.
Checks for every pair e=(a,b):
 (1) M uses e  iff  d(e) < W1 - W0   (W0: best assignment using e with d(e) counted 0; W1: best avoiding e)
 (2) W0 = min_c [d(a,c)+d(b,c)+cost(M_c)]
 (3) every rotation of every M_c is a perfect assignment avoiding e with cost <= d(a,c)+d(b,c)+cost(M_c)+N
 (4) rho_e <= N for every rotation of M_{c*}
"""
import itertools, sys
import numpy as np
TOL = 1e-9

def all_assign(q):
    perms = list(itertools.permutations(range(q)))
    return perms

def run(q, seed):
    rng = np.random.default_rng(seed)
    dAB, dBC, dAC = rng.random((q, q)), rng.random((q, q)), rng.random((q, q))
    perms = all_assign(q)
    P = np.array(perms)                                  # (q!, q)
    ar = np.arange(q)
    # cost[s,t] : team i = (i, P[s][i], P[t][i])
    cAB = dAB[ar, P].sum(1)                              # (q!,)
    cAC = dAC[ar, P].sum(1)
    cBC = dBC[P[:, None, :], P[None, :, :]].sum(2)       # (q!, q!)
    cost = cAB[:, None] + cAC[None, :] + cBC
    s0, t0 = np.unravel_index(cost.argmin(), cost.shape)
    Tcost = cost[s0, t0]
    stats = dict(pairs=0, thr_viol=0, w0_viol=0, rot=0, rot_viol=0, rho_viol=0, invalid=0, min_gap=np.inf)
    for a in range(q):
        for b in range(q):
            uses = (P[:, a] == b)                        # sigma uses (a,b)
            W0 = cost[uses, :].min() - dAB[a, b]
            W1 = cost[~uses, :].min()
            rho = W1 - W0
            Muses = (P[s0, a] == b)
            stats['pairs'] += 1
            if Muses != (dAB[a, b] < rho): stats['thr_viol'] += 1
            # decoupled instances
            best = np.inf; cstar = None; Mc = {}
            A_ = [i for i in range(q) if i != a]; B_ = [i for i in range(q) if i != b]
            for c in range(q):
                C_ = [i for i in range(q) if i != c]
                bc = np.inf; arg = None
                for sg in itertools.permutations(B_):
                    for tg in itertools.permutations(C_):
                        cc = sum(dAB[A_[i], sg[i]] + dBC[sg[i], tg[i]] + dAC[A_[i], tg[i]] for i in range(q - 1))
                        if cc < bc: bc = cc; arg = (sg, tg)
                Mc[c] = (bc, [(A_[i], arg[0][i], arg[1][i]) for i in range(q - 1)])
                val = dAC[a, c] + dBC[b, c] + bc
                if val < best: best = val; cstar = c
            if abs(best - W0) > TOL: stats['w0_viol'] += 1
            for c in range(q):
                bc, teams = Mc[c]
                base = dAC[a, c] + dBC[b, c] + bc
                for k in range(1, q):
                    for js in itertools.permutations(range(q - 1), k):
                        new = list(teams)
                        N = 0.0
                        # b into team j1
                        aj, bj, cj = teams[js[0]]
                        new[js[0]] = (aj, b, cj); N += dAB[aj, b] + dBC[b, cj]
                        for i in range(k - 1):
                            src = teams[js[i]][1]; al, bl, cl = teams[js[i + 1]]
                            new[js[i + 1]] = (al, src, cl); N += dAB[al, src] + dBC[src, cl]
                        last = teams[js[-1]][1]
                        new.append((a, last, c)); N += dAB[a, last] + dBC[last, c]
                        # validity
                        ok = (sorted(t[0] for t in new) == list(range(q)) and sorted(t[1] for t in new) == list(range(q))
                              and sorted(t[2] for t in new) == list(range(q)) and (a, b) not in [(t[0], t[1]) for t in new])
                        if not ok: stats['invalid'] += 1
                        cst = sum(dAB[x, y] + dBC[y, z] + dAC[x, z] for x, y, z in new)
                        stats['rot'] += 1
                        if cst > base + N + TOL: stats['rot_viol'] += 1
                        if W1 > cst + TOL: stats['rot_viol'] += 1      # W1 must be <= any avoiding assignment
                        if c == cstar:
                            stats['min_gap'] = min(stats['min_gap'], N - rho)
                            if rho > N + TOL: stats['rho_viol'] += 1
    return stats

if __name__ == '__main__':
    for q, n in ((3, 200), (4, 60), (5, 6)):
        tot = None
        for i in range(n):
            st = run(q, 1000 * q + i)
            if tot is None: tot = st
            else:
                for k in st: tot[k] = min(tot[k], st[k]) if k == 'min_gap' else tot[k] + st[k]
        print(q, n, tot, flush=True)
