"""V16 independent check: recursive tree process on truncated Poisson hypertrees.
Forest of nroots trees, depth d; agent at level l keeps child teams of cost <= Ks[l]
and gets a tail message W with P(W > x) = exp(-G_K(x)) (the claim of Prop. A4, tested here
by checking that the root message keeps law mu).  Boundary messages iid mu.
Checks: (i) root message law; (ii) M defined by xi < sum of three messages agrees with
'argmin at a', 'argmin at b', 'argmin at c'; (iii) each interior agent in exactly one team of
M or matched in its tail; (iv) Palm law of outward messages; (v) mass transport root <-> depth-1
for the tree WITH messages; (vi) cost of the root's team."""
import sys, numpy as np
from v16_fixed_point import solve, G_from_py

def make_sampler(x, Gk, rng):
    i0 = np.argmax(Gk > 1e-12)
    gx = Gk[max(i0-1, 0):]; xx = x[max(i0-1, 0):]
    assert np.all(np.diff(gx) > 0)
    def samp(n):
        return np.interp(rng.exponential(size=n), gx, xx)
    return samp

def run(nroots, depth, Ks, d, rng, x0=0.3):
    x, G, y, py = d['x'], d['G'], d['y'], d['py']
    samp_mu = make_sampler(x, G, rng)
    sampW = {K: make_sampler(x, G_from_py(x, py, y, K=K), rng) for K in set(Ks)}
    lev = []
    n = nroots
    for l in range(depth):
        K = Ks[l]
        deg = rng.poisson(K**3/6, size=n)
        nt = int(deg.sum())
        par = np.repeat(np.arange(n), deg)
        xi = K * rng.random(nt)**(1/3)
        lev.append(dict(n=n, deg=deg, par=par, xi=xi, W=sampW[K](n)))
        n = 2*nt
    Xin = samp_mu(n)
    for l in reversed(range(depth)):
        L = lev[l]
        L['c1'], L['c2'] = Xin[0::2], Xin[1::2]
        L['Zin'] = L['xi'] - L['c1'] - L['c2']
        X = L['W'].copy()
        np.minimum.at(X, L['par'], L['Zin'])
        L['Xin'] = X
        Xin = X
    # outward pass
    Zpar = np.full(nroots, np.inf)
    out = dict(mism=0, nteams_checked=0, bad_agents=0, agents_checked=0, tail_matched=0)
    for l in range(depth):
        L = lev[l]; n = L['n']; par = L['par']; Zin = L['Zin']; nt = len(Zin)
        W = L['W']
        m1 = np.minimum(W, Zpar); m2 = np.maximum(W, Zpar)
        barg = np.where(W <= Zpar, -1, -2)
        order = np.lexsort((Zin, par))
        ps = par[order]
        first = np.ones(nt, bool); first[1:] = ps[1:] != ps[:-1]
        c1 = np.full(n, np.inf); c2 = np.full(n, np.inf); a1 = np.full(n, -3)
        fi = np.nonzero(first)[0]
        c1[ps[fi]] = Zin[order[fi]]; a1[ps[fi]] = order[fi]
        sec = fi + 1
        ok = (sec < nt); sec = sec[ok]
        ok2 = ps[sec] == ps[sec-1]
        sec = sec[ok2]
        c2[ps[sec]] = Zin[order[sec]]
        F = np.minimum(m1, c1)
        second = np.where(m1 < c1, np.minimum(m2, c1), np.minimum(m1, c2))
        arg = np.where(m1 < c1, barg, a1)
        Xout = np.where(arg[par] == np.arange(nt), second[par], F[par])
        inM = L['xi'] < Xout + L['c1'] + L['c2']
        L.update(F=F, arg=arg, Xout=Xout, inM=inM, Zpar=Zpar)
        # consistency at the parent agent
        out['mism'] += int(np.sum(inM != (arg[par] == np.arange(nt))))
        # agent-level count using the SUM definition of M
        nM = np.bincount(par[inM], minlength=n).astype(int)
        if l > 0:
            P = lev[l-1]
            pin = np.repeat(P['inM'], 2)       # parent team in M (sum definition, evaluated above)
            nM = nM + pin
            # consistency at the two children of each level-(l-1) team
            out['mism'] += int(np.sum(pin != (arg == -2)))
            out['nteams_checked'] += len(P['inM'])
        out['bad_agents'] += int(np.sum(nM + (arg == -1) != 1))
        out['agents_checked'] += n
        out['tail_matched'] += int(np.sum(arg == -1))
        L['nM'] = nM
        # children parent-candidates
        Zp = np.empty(2*nt)
        Zp[0::2] = L['xi'] - Xout - L['c2']
        Zp[1::2] = L['xi'] - Xout - L['c1']
        Zpar = Zp
    R = lev[0]
    res = dict(out)
    res['Xroot'] = R['F']
    K0 = Ks[0]
    res['nM_root'] = R['nM']
    res['cost_root'] = np.bincount(R['par'], weights=R['xi']*R['inM'], minlength=nroots)
    # Palm law of outward message at root: sum_t 1{Xout>x0}, joint with children
    res['palm0'] = np.bincount(R['par'], weights=(R['Xout'] > x0), minlength=nroots)
    res['palm0j'] = np.bincount(R['par'], weights=(R['Xout'] > x0)&(R['c1'] > x0)&(R['c2'] > x0), minlength=nroots)
    if depth >= 2:
        L1 = lev[1]
        root_of_b = np.repeat(R['par'], 2)             # root of each depth-1 agent
        root_of_t1 = root_of_b[L1['par']]
        res['palm1'] = np.bincount(root_of_t1, weights=(L1['Xout'] > x0), minlength=nroots)
        # mass transport tests (need Ks[0]==Ks[1]): g(v) = cost of v's M-team in ball; phi = 1{F_v > x0}
        costb = np.bincount(L1['par'], weights=L1['xi']*L1['inM'], minlength=L1['n']) \
                + np.repeat(R['xi']*R['inM'], 2)
        res['mt_cost_lhs'] = np.bincount(root_of_b, weights=costb, minlength=nroots)
        res['mt_cost_rhs'] = res['cost_root'] * 2 * R['deg']
        res['mt_F_lhs'] = np.bincount(root_of_b, weights=(L1['F'] > x0), minlength=nroots)
        res['mt_F_rhs'] = (R['F'] > x0) * 2 * R['deg']
        res['mt_n_lhs'] = np.bincount(root_of_b, weights=L1['nM'], minlength=nroots)
        res['mt_n_rhs'] = R['nM'] * 2 * R['deg']
        # non-symmetric transport: f(o,b) = xi_e * 1{F_o > x0} * 1{b's M-team is in b's ball and is not e}
        other = (L1['nM'] == 1) & ~np.repeat(R['inM'], 2)
        xie = np.repeat(R['xi'], 2)
        res['mt_a_lhs'] = (R['F'] > x0) * np.bincount(root_of_b, weights=xie*other, minlength=nroots)
        rootother = (R['nM'] == 1)
        # reversed roles: xi_e * 1{F_b > x0} * 1{o's M-team in ball and not e}
        res['mt_a_rhs'] = np.bincount(root_of_b, weights=xie*(L1['F'] > x0)*(~np.repeat(R['inM'], 2)), minlength=nroots) * rootother
    return res

def ms(a):
    a = np.asarray(a, float); return a.mean(), a.std(ddof=1)/np.sqrt(len(a))

if __name__ == "__main__":
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 20261112)
    d = solve(h=0.005)
    x, G = d['x'], d['G']
    EX = (d['pm']*d['xm']).sum(); sp = np.maximum(d['s'], 0)
    x0 = 0.3
    pX = float(np.exp(-np.interp(x0, x, G)))
    def pred(K):
        m = np.minimum(sp, K)
        return (d['ps']*m**3).sum()/6, (d['ps']*m**4).sum()/8
    configs = [
        ("A: depth 2, K=(5,3), 40000 roots (cost at the root)", 40000, 2, [5.0, 3.0]),
        ("B: depth 3, K=3 uniform, 20000 roots (mass transport)", 20000, 3, [3.0, 3.0, 3.0]),
        ("C: depth 6, K=(3,3,2.5,2.5,2.5,2.5), 400 roots", 400, 6, [3.0, 3.0, 2.5, 2.5, 2.5, 2.5]),
        ("D: depth 8, K=2.2 uniform, 300 roots", 300, 8, [2.2]*8),
    ]
    for name, nr, dep, Ks in configs:
        r = run(nr, dep, Ks, d, rng, x0=x0)
        print("==", name)
        print(f"  agents checked {r['agents_checked']}, teams checked at children {r['nteams_checked']}, "
              f"argmin/sum-definition mismatches {r['mism']}, agents not in exactly one (M-team or tail) {r['bad_agents']}, "
              f"tail-matched share {r['tail_matched']/r['agents_checked']:.4f}")
        m, s = ms(r['Xroot']); print(f"  root message mean {m:.4f} +- {s:.4f}   (EX = {EX:.4f});  "
              f"P(X_root > {x0}) {np.mean(r['Xroot'] > x0):.4f} (mu: {pX:.4f}); "
              f"P(X_root>-0.5) {np.mean(r['Xroot'] > -0.5):.4f} (mu: {np.exp(-np.interp(-0.5, x, G)):.4f}); "
              f"P(X_root>1) {np.mean(r['Xroot'] > 1):.4f} (mu: {np.exp(-np.interp(1.0, x, G)):.4f})")
        p3, p4 = pred(Ks[0])
        m, s = ms(r['nM_root']); print(f"  P(root matched within K0) {m:.4f} +- {s:.4f}  (pred E min(S+,K)^3/6 = {p3:.4f})")
        m, s = ms(r['cost_root']); print(f"  cost of root's team within K0 {m:.4f} +- {s:.4f}  (pred E min(S+,K)^4/8 = {p4:.4f}; full 2.2495)")
        m, s = ms(r['palm0']); print(f"  Palm outward at root: {m:.4f} +- {s:.4f}  (pred K^3/6 * P(X>x0) = {Ks[0]**3/6*pX:.4f})")
        m, s = ms(r['palm0j']); print(f"  Palm joint three messages > x0: {m:.4f} +- {s:.4f}  (pred {Ks[0]**3/6*pX**3:.4f})")
        if dep >= 2:
            m, s = ms(r['palm1']); print(f"  Palm outward at depth 1: {m:.4f} +- {s:.4f}  (pred {2*Ks[0]**3/6*Ks[1]**3/6*pX:.4f})")
            if Ks[0] == Ks[1]:
                for k in ("cost", "F", "n", "a"):
                    dlt = r[f'mt_{k}_lhs'] - r[f'mt_{k}_rhs']
                    ml, sl = ms(r[f'mt_{k}_lhs']); mr, sr = ms(r[f'mt_{k}_rhs']); md, sd = ms(dlt)
                    print(f"  MT[{k}]: lhs {ml:.4f} rhs {mr:.4f}  diff {md:.4f} +- {sd:.4f}")
