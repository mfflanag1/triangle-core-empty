"""Gate 1 of the typical-case protection lane (frozen in ../PLAN.md, 2026-11-09).

Random pair-additive three-sector markets, q agents per sector.
graph 'triangle': matrices A (1x2), B (1x3), C (2x3) iid from the value law.
graph 'path':     C == 0 (control; P1 tree theorem says core nonempty, rho >= 1/3, s = 0).

Per instance (floating point, scipy HiGHS; evidence, not certificates):
  1. primal LP  max sum_e lam_e v(e), each player in total weight <= 1, lam >= 0,
     solved by dual simplex so the returned point is a vertex; integrality recorded.
     Exact integer T by scipy milp: always for q <= 6, and for q > 6 whenever the LP
     vertex is fractional (to get the gap).  Core nonempty  <=>  T == LP value.
  2. if the core is nonempty: H_i = max p_i over the core (3q LPs),
     rho = max t s.t. p in core, p_i >= t H_i for H_i > 0,
     s   = min y(N) - T s.t. y(e) >= v(e), y_i >= H_i / 3,
     spread at the rho-optimal point: (a) at the vertex HiGHS returns for the rho LP,
     (b) minimum over the rho-optimal face of the largest within-sector range.
Core: C = {p >= 0 : p(N) = T, p(e) >= v(e) for every full team e} (P1 Sec. 2, P7 Sec. 2).

Seeds: numpy SeedSequence([BASE_SEED, graph_code, dist_code, q, idx]).
"""
import argparse, json, os, sys, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from concurrent.futures import ProcessPoolExecutor

BASE_SEED = 20261111
GRAPH = {'triangle': 0, 'path': 1}
DIST = {'unif': 0, 'exp': 1}
INT_TOL = 1e-6      # LP vertex coordinate counts as integral
GAP_TOL = 1e-7      # relative: LP value - T below this counts as zero gap
H_TOL = 1e-9        # H_i > H_TOL counts as H_i > 0
RHO_SLACK = 1e-9    # slack on t when fixing the rho-optimal face


def instance(graph, dist, q, idx):
    rng = np.random.default_rng(np.random.SeedSequence([BASE_SEED, GRAPH[graph], DIST[dist], q, idx]))
    draw = (lambda: rng.random((q, q))) if dist == 'unif' else (lambda: rng.exponential(1.0, (q, q)))
    A, B, C = draw(), draw(), draw()
    if graph == 'path':
        C = np.zeros((q, q))
    return A, B, C


def structure(A, B, C):
    q = A.shape[0]
    v = (A[:, :, None] + B[:, None, :] + C[None, :, :]).reshape(-1)  # team (a,b,c) -> a*q*q+b*q+c
    a, b, c = np.unravel_index(np.arange(q ** 3), (q, q, q))
    rows = np.repeat(np.arange(q ** 3), 3)
    cols = np.stack([a, q + b, 2 * q + c], axis=1).reshape(-1)
    M = sp.csr_matrix((np.ones(3 * q ** 3), (rows, cols)), shape=(q ** 3, 3 * q))  # team x player
    return v, M


def ok(r, what):
    if not r.success:
        raise RuntimeError(f'{what}: {r.message}')
    return r


def spread(p, q):
    return float(max(np.ptp(p[r * q:(r + 1) * q]) for r in range(3)))


def solve(graph, dist, q, idx, milp_always_upto=6, milp_time=300.0):
    t0 = time.time()
    A, B, C = instance(graph, dist, q, idx)
    v, M = structure(A, B, C)
    n = 3 * q
    rec = dict(graph=graph, dist=dist, q=q, idx=idx)

    # 1. primal LP (dual simplex -> vertex)
    lp = ok(linprog(-v, A_ub=M.T.tocsr(), b_ub=np.ones(n), bounds=(0, None), method='highs-ds'), 'primal LP')
    lam = lp.x
    V = -lp.fun
    integral = bool(np.max(np.minimum(lam, 1 - lam)) < INT_TOL)
    rec.update(V_lp=V, lp_integral=integral, lp_frac_max=float(np.max(np.minimum(lam, 1 - lam))))
    T_round = float(v[lam > 0.5].sum()) if integral else None

    T_milp = None
    if q <= milp_always_upto or not integral:
        mi = milp(-v, constraints=LinearConstraint(M.T.tocsr(), -np.inf, np.ones(n)),
                  integrality=np.ones(q ** 3), bounds=Bounds(0, 1),
                  options=dict(time_limit=milp_time, mip_rel_gap=1e-12))
        rec['milp_status'] = int(mi.status)
        if mi.status == 0:
            T_milp = -mi.fun
        elif mi.x is not None:
            rec['milp_incumbent'] = float(-mi.fun)
    rec['T_milp'] = T_milp

    T = T_milp if T_milp is not None else T_round
    rec['T'] = T
    if T is None:  # fractional vertex and MILP unsolved: a.s.-unique optimum => empty core
        rec['core_nonempty'] = False
        rec['core_basis'] = 'fractional_vertex_only'
        rec['secs'] = time.time() - t0
        return rec
    gap = (V - T) / T
    rec['rel_gap'] = gap
    nonempty = bool(gap <= GAP_TOL)
    rec['core_nonempty'] = nonempty
    rec['core_basis'] = 'milp' if T_milp is not None else 'integral_vertex'
    if T_milp is not None and integral:
        rec['milp_vs_vertex_diff'] = abs(T_milp - T_round)
    if not integral and nonempty:
        rec['flag'] = 'fractional_vertex_but_zero_gap'
    if not nonempty:
        rec['secs'] = time.time() - t0
        return rec

    # 2. core quantities
    negM = -M
    Aeq = np.ones((1, n))
    feas = linprog(np.ones(n), A_ub=negM, b_ub=-v, bounds=(0, None), method='highs')
    rec['dual_min'] = float(ok(feas, 'core min').fun)  # LP duality: equals V
    H = np.empty(n)
    for i in range(n):
        c = np.zeros(n); c[i] = -1
        H[i] = -ok(linprog(c, A_ub=negM, b_ub=-v, A_eq=Aeq, b_eq=[T], bounds=(0, None), method='highs'), f'H_{i}').fun
    rec['H'] = H.tolist()
    pos = H > H_TOL

    # rho LP: variables (p, t)
    Hm = sp.csr_matrix((H[pos], (np.arange(pos.sum()), np.arange(n)[pos])), shape=(pos.sum(), n))
    Ip = sp.csr_matrix((-np.ones(pos.sum()), (np.arange(pos.sum()), np.arange(n)[pos])), shape=(pos.sum(), n))
    A_ub = sp.vstack([sp.hstack([negM, sp.csr_matrix((q ** 3, 1))]),
                      sp.hstack([Ip, sp.csr_matrix(H[pos][:, None])])]).tocsr()
    b_ub = np.concatenate([-v, np.zeros(pos.sum())])
    c = np.zeros(n + 1); c[-1] = -1
    r = ok(linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=np.concatenate([np.ones(n), [0]])[None, :], b_eq=[T],
                   bounds=(0, None), method='highs'), 'rho')
    rho = float(r.x[-1])
    p_rho = r.x[:n]
    rec['rho'] = rho
    rec['spread_vertex'] = spread(p_rho, q)
    rec['p_rho_vertex'] = p_rho.tolist()

    # min over rho-optimal face of the largest within-sector range: vars (p, u0..2, l0..2, z)
    m = n + 7
    rows = [sp.hstack([negM, sp.csr_matrix((q ** 3, 7))])]
    rhs = [-v]
    floor = np.zeros((pos.sum(), m)); floor[np.arange(pos.sum()), np.arange(n)[pos]] = -1
    rows.append(sp.csr_matrix(floor)); rhs.append(-(rho - RHO_SLACK) * H[pos])
    E = np.zeros((2 * n + 3, m)); k = 0
    for s_ in range(3):
        for j in range(q):
            i = s_ * q + j
            E[k, i] = 1; E[k, n + s_] = -1; k += 1        # p_i <= u_s
            E[k, i] = -1; E[k, n + 3 + s_] = 1; k += 1    # l_s <= p_i
        E[k, n + s_] = 1; E[k, n + 3 + s_] = -1; E[k, n + 6] = -1; k += 1  # u_s - l_s <= z
    rows.append(sp.csr_matrix(E)); rhs.append(np.zeros(2 * n + 3))
    c = np.zeros(m); c[-1] = 1
    eq = np.zeros((1, m)); eq[0, :n] = 1
    r2 = ok(linprog(c, A_ub=sp.vstack(rows).tocsr(), b_ub=np.concatenate(rhs), A_eq=eq, b_eq=[T],
                    bounds=[(0, None)] * n + [(None, None)] * 7, method='highs'), 'min spread')
    rec['spread_min'] = float(r2.fun)

    # P7 subsidy
    rs = ok(linprog(np.ones(n), A_ub=negM, b_ub=-v, bounds=[(h / 3, None) for h in H], method='highs'), 's')
    rec['s'] = float(rs.fun - T)
    rec['T_over_sumH'] = float(T / H.sum())
    rec['secs'] = time.time() - t0
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--graph', required=True, choices=list(GRAPH))
    ap.add_argument('--dist', required=True, choices=list(DIST))
    ap.add_argument('--q', type=int, required=True)
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--out', default=os.path.join(os.path.dirname(__file__), 'out'))
    a = ap.parse_args()
    path = os.path.join(a.out, f'{a.graph}_{a.dist}_q{a.q:02d}.jsonl')
    done = set()
    if os.path.exists(path):
        with open(path) as f:
            done = {json.loads(l)['idx'] for l in f}
    todo = [i for i in range(a.start, a.start + a.n) if i not in done]
    t0 = time.time()
    with open(path, 'a') as f, ProcessPoolExecutor(a.workers) as ex:
        for rec in ex.map(solve, [a.graph] * len(todo), [a.dist] * len(todo), [a.q] * len(todo), todo, chunksize=1):
            f.write(json.dumps(rec) + '\n'); f.flush()
    print(f'{path}: {len(todo)} new instances in {time.time() - t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()
