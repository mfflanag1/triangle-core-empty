"""Summarise gate 1 outputs (out/*.jsonl) into out/summary.json and a printed table."""
import glob, json, os, sys
import numpy as np
from scipy.stats import beta

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
TOL = 1e-7


def cp(k, n, a=0.05):
    lo = 0.0 if k == 0 else beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - a / 2, k + 1, n - k)
    return float(lo), float(hi)


def stats(x):
    x = np.asarray(x, float)
    if len(x) == 0:
        return None
    return dict(mean=float(x.mean()), median=float(np.median(x)), min=float(x.min()), max=float(x.max()))


def summarise(recs):
    n = len(recs)
    ne = [r for r in recs if r['core_nonempty']]
    em = [r for r in recs if not r['core_nonempty']]
    k = len(ne)
    d = dict(n=n, nonempty=k, share=k / n, ci95=cp(k, n),
             lp_integral=sum(r['lp_integral'] for r in recs),
             core_basis={b: sum(r['core_basis'] == b for r in recs) for b in {r['core_basis'] for r in recs}},
             flags=[(r['idx'], r['flag']) for r in recs if 'flag' in r],
             milp_unsolved=sum(r.get('milp_status') not in (None, 0) for r in recs),
             rel_gap_all=stats([r['rel_gap'] for r in recs if r.get('rel_gap') is not None]),
             rel_gap_empty=stats([r['rel_gap'] for r in em if r.get('rel_gap') is not None]),
             max_milp_vs_vertex=max([r.get('milp_vs_vertex_diff', 0.0) for r in recs]),
             max_dual_min_err=max([abs(r['dual_min'] - r['V_lp']) / r['V_lp'] for r in ne], default=None),
             secs_total=sum(r['secs'] for r in recs))
    if k:
        rho = np.array([r['rho'] for r in ne]); s = np.array([r['s'] for r in ne]); T = np.array([r['T'] for r in ne])
        d.update(rho=stats(rho), rho_ge_third=int(np.sum(rho >= 1 / 3 - TOL)), rho_lt_third=int(np.sum(rho < 1 / 3 - TOL)),
                 s=stats(s), s_over_T=stats(s / T), s_zero=int(np.sum(s <= TOL * T)),
                 link_violations=int(np.sum((s <= TOL * T) != (rho >= 1 / 3 - TOL))),
                 rho_above_T_over_sumH=int(np.sum(rho > np.array([r['T_over_sumH'] for r in ne]) + TOL)),
                 T_over_sumH=stats([r['T_over_sumH'] for r in ne]),
                 spread_min=stats([r['spread_min'] for r in ne]),
                 spread_vertex=stats([r['spread_vertex'] for r in ne]))
    return d


def main():
    groups = {}
    for f in sorted(glob.glob(os.path.join(OUT, '*.jsonl'))):
        with open(f) as fh:
            recs = [json.loads(l) for l in fh]
        r0 = recs[0]
        groups[f"{r0['graph']}_{r0['dist']}_q{r0['q']:02d}"] = summarise(recs)
    with open(os.path.join(OUT, 'summary.json'), 'w') as fh:
        json.dump(groups, fh, indent=1)
    hdr = f"{'group':22s} {'n':>5s} {'ne':>5s} {'share':>6s} {'ci95':>15s} {'gap_empty_med':>13s} {'rho_mean':>9s} {'rho_min':>8s} {'rho<1/3':>7s} {'s=0':>5s} {'s/T_max':>8s} {'sprMin_med':>10s} {'sprVtx_med':>10s} {'T/sumH':>7s} {'flags':>5s}"
    print(hdr)
    for g, d in groups.items():
        f = lambda x, key='mean', w=9, p=4: (f"{x[key]:{w}.{p}f}" if x else ' ' * (w - 1) + '-')
        print(f"{g:22s} {d['n']:5d} {d['nonempty']:5d} {d['share']:6.3f} [{d['ci95'][0]:.3f},{d['ci95'][1]:.3f}] "
              f"{f(d['rel_gap_empty'], 'median', 13, 5)} {f(d.get('rho'))} {f(d.get('rho'), 'min', 8)} "
              f"{d.get('rho_lt_third', 0):7d} {d.get('s_zero', 0):5d} {f(d.get('s_over_T'), 'max', 8)} "
              f"{f(d.get('spread_min'), 'median', 10)} {f(d.get('spread_vertex'), 'median', 10)} {f(d.get('T_over_sumH'), 'mean', 7, 3)} {len(d['flags']):5d}")


if __name__ == '__main__':
    main()
