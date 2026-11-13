"""V16: batched reruns of v16_tree_sim.run with more roots; corrected two-step Palm prediction."""
import sys, numpy as np
from v16_fixed_point import solve, G_from_py
from v16_tree_sim import run, ms

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
d = solve(h=0.005)
x, G = d['x'], d['G']
EX = (d['pm']*d['xm']).sum(); sp = np.maximum(d['s'], 0)
x0 = 0.3
pX = float(np.exp(-np.interp(x0, x, G)))
configs = [
    ("A: depth 2, K=(6,3)", 20000, 10, 2, [6.0, 3.0]),
    ("B: depth 3, K=3 uniform", 20000, 10, 3, [3.0]*3),
    ("C: depth 6, K=(3,3,2.5,2.5,2.5,2.5)", 400, 12, 6, [3.0, 3.0, 2.5, 2.5, 2.5, 2.5]),
    ("D: depth 8, K=2.2 uniform", 400, 12, 8, [2.2]*8),
]
keys = ['Xroot', 'nM_root', 'cost_root', 'palm0', 'palm0j', 'palm1'] + [f'mt_{k}_{s}' for k in ("cost", "F", "n", "a") for s in ("lhs", "rhs")]
for name, nr, nb, dep, Ks in configs:
    acc = {k: [] for k in keys}; tot = dict(mism=0, bad_agents=0, agents_checked=0, nteams_checked=0)
    for b in range(nb):
        r = run(nr, dep, Ks, d, rng, x0=x0)
        for k in keys: acc[k].append(r[k])
        for k in tot: tot[k] += r[k]
    r = {k: np.concatenate(v) for k, v in acc.items()}
    K0, K1 = Ks[0], Ks[1]
    m = np.minimum(sp, K0); p3 = (d['ps']*m**3).sum()/6; p4 = (d['ps']*m**4).sum()/8
    GK0 = float(np.interp(x0, x, G_from_py(x, d['py'], d['y'], K=K0)))
    palm1_pred = 2*(K1**3/6)*pX*(K0**3/6 - (-np.log(pX) - GK0))
    print("==", name, f"roots {nr*nb}")
    print(f"  agents {tot['agents_checked']}, teams checked from children {tot['nteams_checked']}, mismatches {tot['mism']}, bad agents {tot['bad_agents']}")
    mu_, s = ms(r['Xroot']); print(f"  root message mean {mu_:.4f} +- {s:.4f} (EX {EX:.4f})")
    for t in (-0.5, 0.0, 0.3, 1.0, 1.5):
        emp = np.mean(r['Xroot'] > t); th = np.exp(-np.interp(t, x, G)); se = np.sqrt(th*(1-th)/len(r['Xroot']))
        print(f"    P(X_root > {t}) {emp:.4f}  mu {th:.4f}  z = {(emp-th)/se:+.2f}")
    for lab, key, pr in (("P(root matched within K0)", 'nM_root', p3), ("cost of root team within K0", 'cost_root', p4),
                         ("Palm outward at root", 'palm0', K0**3/6*pX), ("Palm joint three > x0", 'palm0j', K0**3/6*pX**3),
                         ("two-step outward at depth 1", 'palm1', palm1_pred)):
        mu_, s = ms(r[key]); print(f"  {lab}: {mu_:.4f} +- {s:.4f}  pred {pr:.4f}  z = {(mu_-pr)/s:+.2f}")
    if K0 == K1:
        for k in ("cost", "F", "n", "a"):
            md, sd = ms(r[f'mt_{k}_lhs'] - r[f'mt_{k}_rhs']); ml = r[f'mt_{k}_lhs'].mean()
            print(f"  MT[{k}]: lhs {ml:.4f}  lhs-rhs {md:+.4f} +- {sd:.4f}  z = {md/sd:+.2f}")
    print(f"  full-cost prediction E S+^4/8 = {(d['ps']*sp**4).sum()/8:.4f}, 9EX = {9*EX:.4f}", flush=True)
