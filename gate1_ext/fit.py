"""Post-gate extension (2026-11-11): share table with Clopper-Pearson 95% intervals over the full
curve (gate-1 rows from ../gate1/out, extension rows from out/), and binomial-likelihood fits of
simple shapes to the nonempty-core share as a function of q.  Evidence only."""
import json, os, glob
import numpy as np
from scipy.stats import beta
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
rows = []
for src, pat in (('gate1', os.path.join(HERE, '..', 'gate1', 'out', 'triangle_unif_q*.jsonl')),
                 ('ext', os.path.join(HERE, 'out', 'triangle_unif_q*.jsonl'))):
    for f in sorted(glob.glob(pat)):
        R = [json.loads(l) for l in open(f)]
        q = R[0]['q']; n = len(R); k = sum(r['core_nonempty'] for r in R)
        lo = beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
        hi = beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
        gaps = [r['rel_gap'] for r in R if not r['core_nonempty'] and r.get('rel_gap') is not None]
        unsolved = sum(r.get('core_basis') == 'fractional_vertex_only' for r in R)
        flags = sum('flag' in r for r in R)
        lab = src + ('-lp' if any(r.get('core_basis') == 'fractional_vertex_only_lp_mode' for r in R) else '')  # LP-only: no gap
        rows.append(dict(src=lab, q=q, n=n, k=k, share=k / n, lo=lo, hi=hi,
                         gap_med_empty=float(np.median(gaps)) if gaps else None,
                         unsolved_milp=unsolved, zero_gap_fractional=flags))
rows.sort(key=lambda r: r['q'])
q = np.array([r['q'] for r in rows], float); n = np.array([r['n'] for r in rows]); k = np.array([r['k'] for r in rows])

def nll(pfun, th):
    p = np.clip(pfun(q, th), 1e-12, 1 - 1e-12)
    return -np.sum(k * np.log(p) + (n - k) * np.log(1 - p))

# all shapes anchored to share(q=2)=1 being allowed; shift q -> q-2 so share(2)=1 for exp shapes
sig = lambda x: 1 / (1 + np.exp(-x))  # floor f in (0,1); reported params are raw (f = sig(th[2]))
models = {
    'logistic in q: 1/(1+exp(b(q-q0)))': (lambda q, th: 1 / (1 + np.exp(th[0] * (q - th[1]))), [0.2, 12]),
    'exp(-c q^a) (unshifted)': (lambda q, th: np.exp(-abs(th[0]) * q ** abs(th[1])), [0.01, 1.5]),
    'exp(-c (q-2))': (lambda q, th: np.exp(-th[0] * (q - 2)), [0.05]),
    'exp(-c (q-2)^a)': (lambda q, th: np.exp(-abs(th[0]) * (q - 2) ** abs(th[1])), [0.01, 1.5]),
    'logistic in log q: 1/(1+(q/q0)^b)': (lambda q, th: 1 / (1 + (q / abs(th[1])) ** abs(th[0])), [3, 12]),
    'floor + stretched exp: f+(1-f)exp(-c (q-2)^a)': (lambda q, th: sig(th[2]) + (1 - sig(th[2])) * np.exp(-abs(th[0]) * (q - 2) ** abs(th[1])), [0.01, 1.5, -3.0]),
    'floor + logistic: f+(1-f)/(1+exp(b(q-q0)))': (lambda q, th: sig(th[2]) + (1 - sig(th[2])) / (1 + np.exp(th[0] * (q - th[1]))), [0.3, 14, -3.0]),
}
sat = -np.sum([ki * np.log(ki / ni) + (ni - ki) * np.log(1 - ki / ni) for ki, ni in zip(k, n) if 0 < ki < ni])
out = dict(rows=rows, fits={})
for name, (f, th0) in models.items():
    best = min((minimize(lambda th: nll(f, th), np.array(th0) * s, method='Nelder-Mead',
                         options=dict(maxiter=20000, xatol=1e-10, fatol=1e-10)) for s in (0.5, 1, 2, -1, 3)),
               key=lambda r: r.fun)
    dev = 2 * (best.fun - sat)
    pred = f(q, best.x)
    out['fits'][name] = dict(params=best.x.tolist(), floor=float(sig(best.x[2])) if 'floor' in name else None, deviance=dev, dof=len(q) - len(th0), aic=2 * best.fun + 2 * len(th0),
                             pred=dict(zip(q.astype(int).tolist(), np.round(pred, 4).tolist())),
                             extrap={Q: float(f(np.array([Q], float), best.x)[0]) for Q in (50, 64, 100)})
json.dump(out, open(os.path.join(HERE, 'out', 'ext_summary.json'), 'w'), indent=1)
print(f"{'q':>3} {'src':>8} {'n':>5} {'ne':>5} {'share':>6}  CP95            gap_med(empty) unsolved zero-gap-frac")
for r in rows:
    print(f"{r['q']:>3} {r['src']:>8} {r['n']:>5} {r['k']:>5} {r['share']:.3f}  [{r['lo']:.3f},{r['hi']:.3f}]  "
          f"{r['gap_med_empty'] if r['gap_med_empty'] is None else round(r['gap_med_empty'], 5)!s:>12} {r['unsolved_milp']:>4} {r['zero_gap_fractional']:>4}")
for name, d in sorted(out['fits'].items(), key=lambda kv: kv[1]['aic']):
    print(f"\n{name}: params {np.round(d['params'], 4).tolist()} floor {d['floor']}  deviance {d['deviance']:.1f} on {d['dof']} dof  AIC {d['aic']:.1f}")
    print('  pred', d['pred']); print('  extrap', {k_: round(v, 4) for k_, v in d['extrap'].items()})
