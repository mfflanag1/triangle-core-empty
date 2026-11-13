"""Tables for STRUCTURE.md from out/struct_q*.jsonl. Floating point; evidence only."""
import json, glob, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
R = {}
for f in sorted(glob.glob(os.path.join(HERE, 'out', 'struct_q*.jsonl'))):
    rs = [json.loads(l) for l in open(f)]
    R[rs[0]['q']] = rs
Q = sorted(R)
def fit(xs, ys):
    b, a = np.polyfit(np.log(xs), np.log(ys), 1); return b
def boot_fit(qs, arrs, stat, B=400, seed=0):
    rng = np.random.default_rng(seed); bs = []
    for _ in range(B):
        ys = [stat(a[rng.integers(0, len(a), len(a))]) for a in arrs]
        if min(ys) > 0: bs.append(fit(qs, ys))
    return np.percentile(bs, [5, 95])
print('## T1 scalars')
print('| q | n | P(mu*>0) | mean mu* | sd mu* | mean/sd | mean(LP*-T) | mean(LP*-T given >0) | mean(3q-T) | sd(3q-T) | mean(3q-LP*) | M pair deficit x q^(2/3) | secs |')
rows = {}
for q in Q:
    rs = R[q]; mu = np.array([r['mu'] for r in rs]); gap = np.array([r['gap'] for r in rs]); T = np.array([r['T'] for r in rs]); L = np.array([r['LP'] for r in rs])
    pd = np.array([r['M_pair_deficit_mean'] for r in rs])
    rows[q] = dict(mu=mu, gap=gap, def_T=3*q-T, def_L=3*q-L)
    print(f"| {q} | {len(rs)} | {np.mean(mu>1e-9):.3f} | {mu.mean():.5f} | {mu.std(ddof=1):.5f} | {mu.mean()/mu.std(ddof=1):.2f} | {gap.mean():.4f} | {gap[gap>1e-7*3*q].mean():.4f} | {(3*q-T).mean():.3f} | {(3*q-T).std(ddof=1):.3f} | {(3*q-L).mean():.3f} | {pd.mean()*q**(2/3):.3f} | {np.mean([r['secs'] for r in rs]):.1f} |")
print('\n## T1b fits (log-log slope vs q; 90% bootstrap band); q >= 16')
qq = [q for q in Q if q >= 16]
for nm, st in [('mean mu*', lambda a: a.mean()), ('sd mu*', lambda a: a.std(ddof=1)), ('mean/sd mu*', lambda a: a.mean() / a.std(ddof=1)), ('median mu*', np.median)]:
    ys = [st(rows[q]['mu']) for q in qq]
    print(f'- {nm}: values {[round(y,5) for y in ys]}, slope {fit(qq, ys):.2f}, band {boot_fit(qq, [rows[q]["mu"] for q in qq], st).round(2)}')
for key, nm in [('gap', 'mean LP*-T'), ('def_T', 'mean 3q-T'), ('def_L', 'mean 3q-LP*')]:
    ys = [rows[q][key].mean() for q in qq]
    print(f'- {nm}: slope {fit(qq, ys):.3f}, band {boot_fit(qq, [rows[q][key] for q in qq], np.mean).round(3)}; ratio to q^(1/3): {[round(y / q ** (1/3), 3) for y, q in zip(ys, qq)]}')
ys = [rows[q]['def_T'].std(ddof=1) for q in qq]; print(f'- sd 3q-T: {[round(y,3) for y in ys]} slope {fit(qq, ys):.2f}')
ys = [rows[q]['gap'].std(ddof=1) for q in qq]; print(f'- sd LP*-T: {[round(y,4) for y in ys]} slope {fit(qq, ys):.2f}')
for q in Q:
    m = rows[q]['mu']; print(f'  q={q} mu* quantiles 5/10/25/50/75/90%:', np.percentile(m, [5, 10, 25, 50, 75, 90]).round(4))

print('\n## T2 mu*-maximiser w (adjacent vertex of P at M), all instances (empty and nonempty)')
print('| q | ncomp=1 share | kM mean | kM/q mean | kM median | noff/kM | keep-M-pair share of off teams | newpairs/kM | aux 2-section non-bipartite | w half-integral | w denominators seen (top) |')
from collections import Counter
for q in Q:
    rs = R[q]
    nc = np.array([len(r['eff']) for r in rs]); E = [r['eff'][0] for r in rs if len(r['eff']) == 1]
    k = np.array([e['kM'] for e in E]); no = np.array([e['noff'] for e in E]); nk = np.array([e['nkeep'] for e in E]); npairs = np.array([e['newpairs'] for e in E])
    hi = np.mean([r['w_denoms'] == [2] for r in rs])
    den = Counter(d for r in rs for d in r['w_denoms']).most_common(6)
    print(f"| {q} | {np.mean(nc==1):.2f} | {k.mean():.2f} | {np.mean(k/q):.2f} | {np.median(k):.0f} | {np.mean(no/k):.2f} | {nk.sum()/no.sum():.2f} | {np.mean(npairs/k):.2f} | {np.mean([not e['aux_bipartite'] for e in E]):.2f} | {hi:.2f} | {den} |")
print('\n(empty-core subset of T2)')
for q in Q:
    rs = [r for r in R[q] if r['mu'] > 1e-9]
    if not rs: continue
    E = [r['eff'][0] for r in rs]; k = np.array([e['kM'] for e in E])
    print(f'  q={q} n_empty={len(rs)} kM/q mean {np.mean(k/q):.2f} median kM {np.median(k):.0f} d(w) mean {np.mean([e["d"] for e in E]):.2f} min kM {k.min()} max {k.max()} half-integral {np.mean([r["w_denoms"]==[2] for r in rs]):.2f}')

print('\n## T3 LP-optimal vertex lam* - M, empty-core instances')
print('| q | n_empty | lam* half-integral | ncomp mean | ncomp dist | largest comp kM/q | 2nd comp kM (mean) | largest comp share of gain | any comp with gain<=0 | best comp ratio gain/d >0 share |')
for q in Q:
    rs = [r for r in R[q] if 'lpcomp' in r and r['gap'] > 1e-7 * 3 * q]
    if not rs: continue
    nc = np.array([len(r['lpcomp']) for r in rs])
    big = []; sec = []; share = []; nonpos = 0; bestpos = 0
    for r in rs:
        cs = sorted(r['lpcomp'], key=lambda c: -c['kM'])
        big.append(cs[0]['kM'] / q); sec.append(cs[1]['kM'] if len(cs) > 1 else 0)
        g = np.array([c['gain'] for c in cs]); share.append(g.max() / g.sum())
        nonpos += (g <= 1e-9).any(); bestpos += (g.max() > 1e-9)
    print(f"| {q} | {len(rs)} | {np.mean([r['lam_denoms']==[2] for r in rs]):.2f} | {nc.mean():.2f} | {dict(Counter(nc))} | {np.mean(big):.2f} | {np.mean(sec):.2f} | {np.mean(share):.3f} | {nonpos} | {bestpos/len(rs):.2f} |")
print('\nlam* denominators (empty instances):')
for q in Q:
    rs = [r for r in R[q] if r['gap'] > 1e-7 * 3 * q]
    print(f'  q={q}', Counter(d for r in rs for d in r['lam_denoms']).most_common(8), ' lcm>2 share', np.mean([r['lam_denoms'] != [2] for r in rs]).round(2))

print('\n## T4 near-tight off-M teams per agent (full-LP dual prices p_lp; reduced cost r_f = v_f - p_lp(f) > -eps)')
eps = ['0.005', '0.01', '0.02', '0.05', '0.1', '0.2']
print('| q | ' + ' | '.join(f'eps={e}' for e in eps) + ' | eps=0.5q^-2/3 | eps=q^-2/3 | min per agent (eps=q^-2/3) |')
NT = {}
for q in Q:
    rs = R[q]
    vals = [np.mean([r['neartight']['lp'][e][0] for r in rs]) for e in eps]
    sc = [np.mean([r['neartight']['lp'][c][0] for r in rs]) for c in ('c0.5', 'c1.0')]
    mn = np.mean([r['neartight']['lp']['c1.0'][1] for r in rs])
    NT[q] = vals + sc
    print(f'| {q} | ' + ' | '.join(f'{x:.2f}' for x in vals) + f' | {sc[0]:.2f} | {sc[1]:.2f} | {mn:.2f} |')
print('slopes in q (q>=16):', {e: round(fit(qq, [NT[q][i] for q in qq]), 2) for i, e in enumerate(eps + ['0.5q^-2/3', 'q^-2/3'])})
print('slope in eps at q=40 (or max q):', round(fit([float(e) for e in eps], NT[Q[-1]][:6]), 2))
print('\n(same with the cone dual p of the mu* LP)')
for q in Q:
    rs = R[q]
    print(f'  q={q}', [round(np.mean([r['neartight']['cone'][e][0] for r in rs]), 2) for e in eps])

print('\n## T5 M-team reduced costs')
print('| q | empty: #M-teams with r_lp<0 (mean) / q | empty: share of deficit on top 25% M teams | cone dual: #M-teams at r=-mu* (= supp z) mean / q |')
for q in Q:
    rs = [r for r in R[q] if r['gap'] > 1e-7 * 3 * q]
    if not rs: continue
    neg = np.mean([np.sum(np.array(r['rlp_M']) < -1e-9) / q for r in rs])
    sh = np.mean([(-np.array(r['rlp_M'])[: max(1, q // 4)]).sum() / (-np.array(r['rlp_M'])).sum() for r in rs])
    tight = np.mean([np.sum(np.abs(np.array(r['rc_M']) + r['mu']) < 1e-9) / q for r in rs])
    print(f'| {q} | {neg:.2f} | {sh:.2f} | {tight:.2f} |')
print('\nchecks: max dual feas violation', max(r['dual_feas_violation'] for q in Q for r in R[q]),
      ' max |mu - mu_dual|', max(abs(r['mu'] - r['mu_dual']) for q in Q for r in R[q]),
      ' max |mu - ratio(w)|', max(abs(r['mu'] - r['w_ratio_check']) for q in Q for r in R[q]),
      ' max |rc on supp(w)\\M|', max(r['rc_supp_off_max_abs'] for q in Q for r in R[q]))

print('\n## T6 gap as its own concentrating scalar; maximiser size')
print('| q | n | mean gap | sd gap | mean/sd gap | P(gap>0) | mean mu* x q | maximiser d(w)/q (empty) | maximiser gain v.(w-M) (empty) |')
for q in Q:
    rs = R[q]; g = np.array([r['gap'] for r in rs]); mu = np.array([r['mu'] for r in rs])
    E = [r['eff'][0] for r in rs if r['mu'] > 1e-9]
    print(f"| {q} | {len(rs)} | {g.mean():.4f} | {g.std(ddof=1):.4f} | {g.mean()/g.std(ddof=1):.2f} | {np.mean(g>1e-7*3*q):.3f} | {mu.mean()*q:.3f} | {np.mean([e['d']/q for e in E]):.3f} | {np.mean([e['gain'] for e in E]):.4f} |")
ys = [R[q] for q in qq]
print('- mean/sd gap slope (q>=16):', round(fit(qq, [np.mean([r['gap'] for r in R[q]]) / np.std([r['gap'] for r in R[q]], ddof=1) for q in qq]), 2),
      'band', boot_fit(qq, [np.array([r['gap'] for r in R[q]]) for q in qq], lambda a: a.mean() / a.std(ddof=1)).round(2))

for nm in ('blocks', 'influence'):
    fs = sorted(glob.glob(os.path.join(HERE, 'out', f'{nm}_q*.jsonl')))
    if not fs: continue
    print(f'\n## {nm}')
    for f in fs:
        rs = [json.loads(l) for l in open(f)]
        if not rs: continue
        q = rs[0]['q']
        if nm == 'blocks':
            ks = [k for k in rs[0] if k.startswith('s')]
            emp = [r for r in rs if r['gap'] > 1e-7 * 3 * q]
            print(f'  q={q} n={len(rs)} (empty {len(emp)}):', {k: dict(share_pos_blocks=round(sum(r[k]['pos'] for r in rs) / sum(r[k]['blocks'] for r in rs), 4),
                  P_any_pos_given_empty=round(np.mean([r[k]['any_pos'] for r in emp]), 3) if emp else None,
                  sum_over_gap=round(np.mean([r[k]['sum_mean'] / r['gap'] for r in emp]), 3) if emp else None) for k in ks})
        else:
            g = np.array([r['gap'] for r in rs]); mu = np.array([r['mu'] for r in rs])
            print(f'  q={q} n={len(rs)}: mean ES_gap {np.mean([r["ES_gap"] for r in rs]):.4g} (sample Var gap {g.var(ddof=1):.4g}); mean ES_mu {np.mean([r["ES_mu"] for r in rs]):.4g} (sample Var mu {mu.var(ddof=1):.4g});'
                  f' S1 size {np.mean([r["S1"]["size"] for r in rs]):.0f}, S1 E dgap^2 {np.mean([r["S1"]["msq_gap"] for r in rs]):.3g}, S2 E dgap^2 {np.mean([r["S2"]["msq_gap"] for r in rs]):.3g},'
                  f' S1 E dmu^2 {np.mean([r["S1"]["msq_mu"] for r in rs]):.3g}, S2 E dmu^2 {np.mean([r["S2"]["msq_mu"] for r in rs]):.3g}, S2 share nonzero dgap {np.mean([r["S2"]["frac_nonzero_gap"] for r in rs]):.3f}')

fs = sorted(glob.glob(os.path.join(HERE, 'out', 'influence2_q*.jsonl')))
if fs:
    print('\n## influence2 (Efron-Stein estimates for T, LP*, gap, mu*; Var columns from the 200-instance struct runs)')
    print('| q | n | ES_T | Var T | ES_LP | Var LP* | ES_gap | Var gap | (mean gap)^2 | ES_mu | Var mu* | (mean mu*)^2 |')
    for f in fs:
        rs = [json.loads(l) for l in open(f)]
        if not rs: continue
        q = rs[0]['q']; S = R[q]
        T = np.array([r['T'] for r in S]); L = np.array([r['LP'] for r in S]); g = np.array([r['gap'] for r in S]); mu = np.array([r['mu'] for r in S])
        print(f"| {q} | {len(rs)} | {np.mean([r['ES_T'] for r in rs]):.3g} | {T.var(ddof=1):.3g} | {np.mean([r['ES_LP'] for r in rs]):.3g} | {L.var(ddof=1):.3g} | {np.mean([r['ES_gap'] for r in rs]):.3g} | {g.var(ddof=1):.3g} | {g.mean()**2:.3g} | {np.mean([r['ES_mu'] for r in rs]):.3g} | {mu.var(ddof=1):.3g} | {mu.mean()**2:.3g} |")
