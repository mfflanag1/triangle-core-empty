"""Structure of the stability-margin maximiser and of the LP optimum relative to the integer optimum M.
Random triangle market, gate1 seeds. Floating point (scipy HiGHS); evidence only.

mu* is computed by ONE LP (tangent cone of the perfect axial polytope P at M, normalised by d):
  max  v.y - v_M.z   s.t. for every agent i: sum_{f ni i, f notin M} y_f = z_{M(i)},  sum_e z_e = 1,  y,z >= 0.
Dual: mu* = min mu s.t. exists p (free) with p(f) >= v_f (f notin M), p(e) <= v_e + mu (e in M).
(Equivalent to Lemma 2: w - M for w in P spans the tangent cone; ratio maximised on an extreme ray = edge at M.)
Integer optimum M: reduced-cost fixing with the full-LP dual p_lp (v.x = LP* + sum_f r_f x_f for perfect x,
r_f = v_f - p_lp(f) <= 0), so every optimal team has r_f >= T_inc - LP*; MILP on the restricted set,
enlarged until T_inc >= LP* - delta."""
import json, sys, os, time, argparse
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'gate1'))
from gate1 import instance, structure, ok
TOL = 1e-7

def integer_opt(v, Mm, n, LPs, rlp):
    delta = 0.15
    while True:
        cand = np.nonzero(rlp >= -delta)[0]
        A = Mm[cand].T.tocsr()
        mi = milp(-v[cand], constraints=LinearConstraint(A, np.ones(n), np.ones(n)), integrality=np.ones(len(cand)),
                  bounds=Bounds(0, 1), options=dict(time_limit=600, mip_rel_gap=1e-12))
        if mi.status == 0 and -mi.fun >= LPs - delta + 1e-9:
            sel = np.zeros(len(v), bool); sel[cand[mi.x > 0.5]] = True
            return sel, -mi.fun, len(cand), delta
        delta *= 2

def uf_components(teams, q):
    par = list(range(3 * q))
    def f(a):
        while par[a] != a:
            par[a] = par[par[a]]; a = par[a]
        return a
    for (x, y, z) in teams:
        for b in (q + y, 2 * q + z):
            ra, rb = f(x), f(b)
            if ra != rb: par[ra] = rb
    return [f(x) for (x, y, z) in teams]

def describe(dirv, sel, tm_of_agent, q, v):
    """dirv = w - M (full team vector). Returns structural stats incl. components."""
    sup = np.nonzero(np.abs(dirv) > TOL)[0]
    teams = [tuple(int(t) for t in np.unravel_index(j, (q, q, q))) for j in sup]
    roots = uf_components(teams, q)
    comps = {}
    for j, t, r in zip(sup, teams, roots):
        comps.setdefault(r, []).append((j, t))
    out = []
    for r, lst in comps.items():
        js = np.array([j for j, _ in lst])
        offs = [(j, t) for j, t in lst if not sel[j]]
        Ms = [j for j, t in lst if sel[j]]
        gain = float(v[js] @ dirv[js]); d = float(-dirv[Ms].sum()) if Ms else 0.0
        # auxiliary graph on M-teams: each off team links the M-teams of its agents
        nkeep = 0; edges = set(); newpairs = set()
        for j, (x, y, z) in offs:
            ms = (tm_of_agent[x], tm_of_agent[q + y], tm_of_agent[2 * q + z])
            if len(set(ms)) == 2: nkeep += 1
            u = sorted(set(ms))
            for a in range(len(u)):
                for b in range(a + 1, len(u)): edges.add((u[a], u[b]))
            for pr, (a, b) in (('xy', (x, y)), ('xz', (x, z)), ('yz', (y, z))):
                ia, ib = {'xy': (x, q + y), 'xz': (x, 2 * q + z), 'yz': (q + y, 2 * q + z)}[pr]
                if tm_of_agent[ia] != tm_of_agent[ib]: newpairs.add((pr, a, b))
        nodes = sorted({tm_of_agent[i] for _, (x, y, z) in offs for i in (x, q + y, 2 * q + z)})
        # bipartiteness of the 2-section
        adj = {u: [] for u in nodes}
        for a, b in edges: adj[a].append(b); adj[b].append(a)
        col = {}; bip = True
        for s0 in nodes:
            if s0 in col: continue
            col[s0] = 0; st = [s0]
            while st:
                u = st.pop()
                for w_ in adj[u]:
                    if w_ not in col: col[w_] = 1 - col[u]; st.append(w_)
                    elif col[w_] == col[u]: bip = False
        out.append(dict(kM=len(Ms), nodesM=len(nodes), noff=len(offs), nkeep=nkeep, nedges_aux=len(edges),
                        newpairs=len(newpairs), gain=gain, d=d, aux_bipartite=bip))
    return out

def denoms(x):
    fr = x[(x > TOL) & (x < 1 - TOL)]
    ds = set()
    for w in fr:
        for k in range(2, 25):
            if abs(w * k - round(w * k)) < 1e-6: ds.add(k); break
        else: ds.add(-1)
    return sorted(ds)

def run(args):
    q, idx = args
    t0 = time.time()
    A, B, C = instance('triangle', 'unif', q, idx)
    v, Mm = structure(A, B, C); n = 3 * q; Ae = Mm.T.tocsr()
    lp = ok(linprog(-v, A_eq=Ae, b_eq=np.ones(n), bounds=(0, None), method='highs-ds'), 'LP')
    lam, LPs = lp.x, -lp.fun
    p_lp = -lp.eqlin.marginals
    rlp = v - Mm @ p_lp
    integral = np.max(np.minimum(lam, 1 - lam)) < 1e-6
    if integral:
        sel = lam > 0.5; T = LPs; ncand = None
    else:
        sel, T, ncand, _ = integer_opt(v, Mm, n, LPs, rlp)
    Mi = np.nonzero(sel)[0]
    tm = np.array(np.unravel_index(Mi, (q, q, q))).T
    tm_of_agent = np.empty(3 * q, int)
    for k, (x, y, z) in enumerate(tm):
        tm_of_agent[x] = k; tm_of_agent[q + y] = k; tm_of_agent[2 * q + z] = k
    # cone LP for mu*
    sgn = np.where(sel, -1.0, 1.0)
    Acone = sp.vstack([sp.diags(sgn) @ Mm, ]).T.tocsr()
    Acone = sp.vstack([Acone, sp.csr_matrix(sel.astype(float)[None, :])]).tocsr()
    b = np.zeros(n + 1); b[-1] = 1
    cobj = np.where(sel, -v, v)
    cr = ok(linprog(-cobj, A_eq=Acone, b_eq=b, bounds=(0, None), method='highs-ds'), 'cone')
    mu = -cr.fun
    lamd = -cr.eqlin.marginals
    p = lamd[:n]; mu_d = lamd[n]
    yz = cr.x
    zmax = yz[sel].max()
    w = sel.astype(float) + np.where(sel, -yz, yz) / zmax
    dir_w = w - sel
    rec = dict(q=q, idx=idx, LP=LPs, T=T, gap=LPs - T, mu=mu, mu_dual=float(mu_d), lp_integral=bool(integral), ncand=ncand,
               w_ratio_check=float(v @ dir_w / (-(dir_w[sel]).sum())), w_denoms=denoms(w), lam_denoms=denoms(lam))
    rec['eff'] = describe(dir_w, sel, tm_of_agent, q, v)
    if not integral:
        rec['lpcomp'] = describe(lam - sel, sel, tm_of_agent, q, v)
    # reduced costs w.r.t. cone dual p (p(f) >= v_f off M; p(e) <= v_e + mu on M)
    rc = v - Mm @ p
    off = ~sel
    rec['dual_feas_violation'] = float(max(rc[off].max(), (-(rc[sel] + mu)).max()))
    supw_off = (np.abs(dir_w) > TOL) & off
    rec['rc_supp_off_max_abs'] = float(np.abs(rc[supw_off]).max())
    rec['rc_M'] = sorted(float(x) for x in rc[sel])  # = v_e - p(e) >= -mu
    # near-tight counts per agent: off teams with rc > -eps (cone dual) and rlp > -eps (full-LP dual)
    eps_list = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
    nt = {}
    for nm, r in (('cone', rc), ('lp', rlp)):
        cnts = {}
        for e in eps_list:
            m = off & (r > -e)
            per_agent = np.asarray(Mm[m].sum(axis=0)).ravel()
            cnts[str(e)] = [float(per_agent.mean()), float(per_agent.min()), float(per_agent.max()), int(m.sum())]
        nt[nm] = cnts
    # scaled eps = c q^{-2/3}
    for c in (0.25, 0.5, 1.0):
        e = c * q ** (-2 / 3)
        m = off & (rlp > -e)
        pa = np.asarray(Mm[m].sum(axis=0)).ravel()
        nt['lp'][f'c{c}'] = [float(pa.mean()), float(pa.min()), float(pa.max()), int(m.sum())]
    rec['neartight'] = nt
    # M-team deficits under the full-LP dual (sum = T - LP*)
    rec['rlp_M'] = sorted(float(x) for x in rlp[sel])
    # pair values of M: deficits 1 - pair value
    pv = np.concatenate([A[tm[:, 0], tm[:, 1]], B[tm[:, 0], tm[:, 2]], C[tm[:, 1], tm[:, 2]]])
    rec['M_pair_deficit_mean'] = float((1 - pv).mean())
    rec['secs'] = time.time() - t0
    return rec

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--q', type=int, nargs='+'); ap.add_argument('--n', type=int)
    ap.add_argument('--start', type=int, default=0); ap.add_argument('--workers', type=int, default=4)
    a = ap.parse_args()
    for q in a.q:
        with ProcessPoolExecutor(a.workers) as ex, open(os.path.join(HERE, 'out', f'struct_q{q:02d}.jsonl'), 'a') as f:
            for r in ex.map(run, [(q, i) for i in range(a.start, a.start + a.n)]):
                f.write(json.dumps(r) + '\n'); f.flush()
        print('done', q, flush=True)
