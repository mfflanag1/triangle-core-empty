"""Round-3 scout: per-instance measurements on the pair-additive triangle market. Floating point; evidence only.
Reuses ../../scaling/cg.py (column-generation LP, reduced-cost-fixing MILP). Seeds: cg.draw(q, idx, TAG) with TAG below.
Team (x,y,z) is worth A[x,y] + B[x,z] + C[y,z] (cg.py's labelling; same model as the brief up to renaming sectors)."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import sys, json, time, argparse
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scaling'))
import cg
TAG = 770031
XTOL = 1e-7
KMAX = 15          # longest odd cycle searched
BUDGET = 3_000_000 # DFS node budget per (instance, kind)

def integer_opt_sel(vf, q, LPs, r, delta0):
    n = 3 * q; delta = delta0
    while True:
        cand = np.nonzero(r >= -delta)[0]
        mi = milp(-vf[cand], constraints=LinearConstraint(cg.amat(cand, q), np.ones(n), np.ones(n)),
                  integrality=np.ones(len(cand)), bounds=Bounds(0, 1), options=dict(time_limit=900, mip_rel_gap=1e-12))
        if mi.status == 0 and -mi.fun >= LPs - delta + 1e-9:
            return cand[mi.x > 0.5], -mi.fun
        if delta > 3: raise RuntimeError('delta blew up')
        delta *= 2

def shortest_odd_cycles(teams, n):
    """teams: list of 3-tuples of agent ids (0..n-1). Returns (weak, strong): length (number of teams) of the shortest
    odd Berge cycle v1,e1,...,vk,ek (distinct agents, distinct teams, v_i,v_{i+1} in e_i), and of the shortest odd
    cycle in which no team of the cycle contains a third cycle agent (the cycles that witness non-balancedness).
    0 = none of length <= KMAX; -1 = budget exhausted."""
    inc = [[] for _ in range(n)]
    for t, tr in enumerate(teams):
        for a in tr: inc[a].append(t)
    INF = 10 ** 9
    dist = np.full((n, n), INF, dtype=np.int64)   # Berge distance between agents
    for s in range(n):
        if not inc[s]: continue
        d = dist[s]; d[s] = 0; fr = [s]
        while fr:
            nf = []
            for u in fr:
                for t in inc[u]:
                    for w in teams[t]:
                        if d[w] == INF: d[w] = d[u] + 1; nf.append(w)
            fr = nf
    out = []
    for strong in (False, True):
        res = 0
        for k in range(3, KMAX + 1, 2):
            budget = [BUDGET]; found = [False]
            vis = np.zeros(n, bool); used = np.zeros(len(teams), bool); forb = np.zeros(n, np.int64)
            def dfs(s, v, j):
                # at agent v having used j teams; need k teams in total and return to s
                if budget[0] <= 0 or found[0]: return
                budget[0] -= 1
                ds = dist[s]
                for t in inc[v]:
                    if used[t]: continue
                    tr = teams[t]
                    if j == k - 1:
                        if s in tr:
                            if strong:
                                u = [w for w in tr if w != v and w != s][0]
                                if vis[u]: continue
                            found[0] = True; return
                        continue
                    for w in tr:
                        if w == v or w < s or vis[w] or ds[w] > k - j - 1: continue
                        if strong:
                            if forb[w]: continue
                            u = [z for z in tr if z != v and z != w][0]
                            if vis[u]: continue
                            forb[u] += 1
                        used[t] = True; vis[w] = True
                        dfs(s, w, j + 1)
                        used[t] = False; vis[w] = False
                        if strong: forb[u] -= 1
                        if found[0]: return
            for s in range(n):
                if not inc[s]: continue
                vis[s] = True
                dfs(s, s, 0)
                vis[s] = False
                if found[0] or budget[0] <= 0: break
            if found[0]: res = k; break
            if budget[0] <= 0: res = -1; break
        out.append(res)
    return out

def run(args):
    q, idx = args
    t0 = time.time()
    A, B, C = cg.draw(q, idx, TAG)
    V = A[:, :, None] + B[:, None, :] + C[None, :, :]
    vf = V.reshape(-1)
    LPs, p, r, cols, x, it = cg.solve_lp(V, q)
    integral = bool(np.max(np.minimum(x, 1 - x)) < 1e-6)
    if integral:
        Msel = cols[x > 0.5]; T = float(vf[Msel].sum())
    else:
        Msel, T = integer_opt_sel(vf, q, LPs, r, 0.15 if q <= 48 else 0.3)
    s = q ** (2 / 3)
    dA, dB, dC = 1 - A, 1 - B, 1 - C
    mx, my, mz = np.unravel_index(Msel, (q, q, q))
    dM = np.concatenate([dA[mx, my], dB[mx, mz], dC[my, mz]])
    sup = cols[x > XTOL]; w = x[x > XTOL]
    sx, sy, sz = np.unravel_index(sup, (q, q, q))
    LA = np.zeros((q, q)); LB = np.zeros((q, q)); LC = np.zeros((q, q))
    np.add.at(LA, (sx, sy), w); np.add.at(LB, (sx, sz), w); np.add.at(LC, (sy, sz), w)
    MA = np.zeros((q, q)); MB = np.zeros((q, q)); MC = np.zeros((q, q))
    MA[mx, my] = 1; MB[mx, mz] = 1; MC[my, mz] = 1
    load = np.concatenate([LA.ravel(), LB.ravel(), LC.ravel()])
    ind = np.concatenate([MA.ravel(), MB.ravel(), MC.ravel()])
    d = np.concatenate([dA.ravel(), dB.ravel(), dC.ravel()])
    pos = load > XTOL
    tdM = 3 - vf[Msel]; tdL = 3 - vf[sup]
    rec = dict(q=q, idx=idx, LP=LPs, T=T, D=LPs - T, integral=integral,
               defT=3 * q - T, defLP=3 * q - LPs,
               maxdM=float(dM.max() * s), maxdL=float(d[pos].max() * s),
               sumsqM=float((dM ** 2).sum()), sumsqL1=float((load * d ** 2).sum()), sumsqL2=float((load ** 2 * d ** 2).sum()),
               gradD2=float(((load - ind) ** 2 * d ** 2).sum()),
               maxteamM=float(tdM.max() * s), maxteamL=float(tdL.max() * s),
               npairsL=int(pos.sum()), nsup=int(len(sup)))
    if not integral:
        fr = (w < 1 - 1e-6)
        teams = [(int(a), int(q + b), int(2 * q + c)) for a, b, c in zip(sx[fr], sy[fr], sz[fr])]
        ag = set(a for t in teams for a in t)
        rec['nfrac_teams'] = len(teams); rec['nfrac_agents'] = len(ag)
        rec['minw'] = float(w[fr].min())
        rec['half'] = bool(np.all(np.abs(w[fr] - 0.5) < 1e-6))
        tc = time.time()
        rec['odd_weak'], rec['odd_strong'] = shortest_odd_cycles(teams, 3 * q)
        rec['cyc_secs'] = time.time() - tc
    rec['secs'] = time.time() - t0
    return rec

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--n', type=int); ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--workers', type=int, default=2)
    a = ap.parse_args()
    for q in a.q:
        fn = os.path.join(HERE, 'out', f'main_q{q:03d}.jsonl')
        os.makedirs(os.path.dirname(fn), exist_ok=True)
        with ProcessPoolExecutor(a.workers) as ex, open(fn, 'a') as f:
            for rr in ex.map(run, [(q, i) for i in range(a.start, a.start + a.n)]):
                f.write(json.dumps(rr) + '\n'); f.flush()
        print('done', q, flush=True)
