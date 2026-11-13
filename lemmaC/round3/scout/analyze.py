import json, glob, numpy as np, sys
def load(pat):
    out = {}
    for fn in sorted(glob.glob(pat)):
        R = [json.loads(l) for l in open(fn)]
        if R: out[R[0]['q']] = R
    return out
f = lambda x, k=3: f'{x:.{k}f}'
M = load('out/main_q*.jsonl')
print('## 1. Used deficits (pair deficit x q^(2/3); S2 = sum of squared pair deficits)')
print('| q | n | max d, M: mean / med / 95% / max | max d, LP: mean / med / 95% / max | max team def M (mean) | S2(M) | S2(M) q^(1/3) | sum load d^2 (LP) | sum load^2 d^2 (LP) | sum (load-1_M)^2 d^2 | 3q-T |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for q, R in M.items():
    g = lambda k: np.array([r[k] for r in R])
    a, b = g('maxdM'), g('maxdL')
    print(f"| {q} | {len(R)} | {f(a.mean(),2)} / {f(np.median(a),2)} / {f(np.quantile(a,.95),2)} / {f(a.max(),2)} | {f(b.mean(),2)} / {f(np.median(b),2)} / {f(np.quantile(b,.95),2)} / {f(b.max(),2)} | {f(g('maxteamM').mean(),2)} | {f(g('sumsqM').mean())} | {f(g('sumsqM').mean()*q**(1/3),2)} | {f(g('sumsqL1').mean())} | {f(g('sumsqL2').mean())} | {f(g('gradD2').mean())} | {f(g('defT').mean(),2)} |")
print()
print('Variance check (Poincare, constant 4): bound on sd = 2 sqrt(E S2).')
print('| q | sd T observed | 2 sqrt(E S2(M)) | 2 sqrt(E(3q-T)) (fact 1) | sd D observed | 2 sqrt(E sum (load-1_M)^2 d^2) | sqrt(8 E(3q-T)) (fact 2) | E D | Chebyshev bound on P(D=0) with grad D | with fact 2 |')
print('|---|---|---|---|---|---|---|---|---|---|')
for q, R in M.items():
    g = lambda k: np.array([r[k] for r in R])
    ED = g('D').mean()
    print(f"| {q} | {f(g('T').std(ddof=1))} | {f(2*np.sqrt(g('sumsqM').mean()),2)} | {f(2*np.sqrt(g('defT').mean()),2)} | {f(g('D').std(ddof=1))} | {f(2*np.sqrt(g('gradD2').mean()),2)} | {f(np.sqrt(8*g('defT').mean()),2)} | {f(ED)} | {f(4*g('gradD2').mean()/ED**2,1)} | {f(8*g('defT').mean()/ED**2,0)} |")
print()
print('## 3. Odd cycles in the fractional part of the LP support (instances with D > 0)')
print('| q | n (D>0) | weak: share 3 / 5 / 7 / >=9 | weak mean | strong: share 3 / 5 / 7 / >=9 | strong mean | frac teams (mean) | frac agents / 3q | all weights 1/2 | min weight (median) |')
print('|---|---|---|---|---|---|---|---|---|---|')
for q, R in M.items():
    F = [r for r in R if not r['integral']]
    if not F: continue
    row = [str(q), str(len(F))]
    for k in ('odd_weak', 'odd_strong'):
        a = np.array([r[k] for r in F]); bad = (a <= 0).sum()
        sh = [np.mean(a == 3), np.mean(a == 5), np.mean(a == 7), np.mean(a >= 9)]
        row.append(' / '.join(f(s, 2) for s in sh) + (f' (none<=15 or budget: {bad})' if bad else ''))
        row.append(f(a[a > 0].mean(), 2))
    g = lambda k: np.array([r[k] for r in F])
    row += [f(g('nfrac_teams').mean(), 1), f((g('nfrac_agents') / (3 * q)).mean(), 2), f(g('half').mean(), 2), f(np.median(g('minw')), 4)]
    print('| ' + ' | '.join(row) + ' |')
print()
print('## 5. Upper tail of Y = (3q - T)/q^(1/3)')
print('| q | n | mean | sd | sd of T | skew | 50% | 90% | 95% | 99% | max | min | E D | sd D | P(D>0) | (3q-LP*)/q^(1/3) |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for q, R in M.items():
    y = np.array([r['defT'] for r in R]) / q ** (1 / 3); D = np.array([r['D'] for r in R])
    z = (y - y.mean()) / y.std()
    qs = np.quantile(y, [.5, .9, .95, .99])
    print(f"| {q} | {len(R)} | {f(y.mean())} | {f(y.std(ddof=1))} | {f(y.std(ddof=1)*q**(1/3))} | {f((z**3).mean(),2)} | " + ' | '.join(f(v, 2) for v in qs) + f" | {f(y.max(),2)} | {f(y.min(),2)} | {f(D.mean())} | {f(D.std(ddof=1))} | {f(np.mean([not r['integral'] for r in R]))} | {f(np.mean([r['defLP'] for r in R])/q**(1/3))} |")
print()
print('## 2. Sister model: P(LP optimum integral)')
S = {}
for fn in sorted(glob.glob('out/sister_*.jsonl')):
    R = [json.loads(l) for l in open(fn)]
    R = list({r['idx']: r for r in R[::-1]}.values())   # dedupe by idx (keep first)
    if R: S[(R[0]['model'], R[0]['q'])] = R
qs = sorted({k[1] for k in S}); ms = [m for m in ('P', 'IH', 'U') if any(k[0] == m for k in S)]
print('| q | ' + ' | '.join(str(q) for q in qs) + ' |'); print('|---|' + '---|' * len(qs))
for m in ms:
    print(f'| {m}: P(integral) | ' + ' | '.join((f(np.mean([r['integral'] for r in S[(m, q)]]), 3) if (m, q) in S else '') for q in qs) + ' |')
    print(f'| {m}: n | ' + ' | '.join((str(len(S[(m, q)])) if (m, q) in S else '') for q in qs) + ' |')
    row = []
    for q in qs:
        R = [r for r in S.get((m, q), []) if r['D'] is not None and (r['integral'] or r['T'] is not None)]
        R = [r for r in R if r['idx'] >= 400] if (m, q) in S and any(r['D'] is None for r in S[(m, q)]) else R
        row.append(f(np.mean([r['D'] for r in R]), 4) + (f' (n={len(R)})' if len(R) != len(S.get((m, q), [])) else '') if R else '')
    print(f'| {m}: E D | ' + ' | '.join(row) + ' |')
print()
print('## 4. Core geometry when D = 0 (essential core Q, dimension 2q-2; brief\'s core has dimension 2q)')
C = load('out/core_q*.jsonl')
print('| q | n | core nonempty | off-M teams tight at centre (max) | eps* q^(2/3): mean / med | eps*>0 share | facets: mean / med / max | facets sharing 1 pair with M / 0 pairs (mean) | facets/(2q-2) | vertices: med / mean / max | mean width x q^(2/3) |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for q, R in C.items():
    K = [r for r in R if r['core']]
    if not K: continue
    g = lambda k: np.array([r[k] for r in K])
    nv = np.array([r['nvert'] for r in K if 'nvert' in r])
    print(f"| {q} | {len(R)} | {len(K)} | {g('tight_at_centre').max()} | {f(g('eps_scaled').mean())} / {f(np.median(g('eps_scaled')))} | {f(np.mean(g('eps')>1e-9),2)} | {f(g('nfacets').mean(),1)} / {np.median(g('nfacets')):.0f} / {g('nfacets').max()} | {f(g('facet_keep1').mean(),1)} / {f(g('facet_keep0').mean(),1)} | {f(g('nfacets').mean()/(2*q-2),2)} | " + (f"{np.median(nv):.0f} / {nv.mean():.0f} / {nv.max()}" if len(nv) else 'n/a') + f" | {f(g('width_scaled').mean(),2)} |")

print()
print('## 1b. Tails of used pair deficits (same instances): share of the 3q used pairs with d q^(2/3) > t')
E = load('out/extra_q*.jsonl')
print('| q | n | M: t=1 | 2 | 3 | 4 | 5 | 6 | LP load: t=1 | 2 | 3 | 4 | 5 | 6 | share of S2(M) from d q^(2/3) > 2 | > 3 | > 4 |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for q, R in E.items():
    tm = np.array([r['tailM'] for r in R]).mean(0) / (3 * q); tl = np.array([r['tailL'] for r in R]).mean(0) / (3 * q)
    sq = np.array([r['sqM'] for r in R]).mean(0)
    print(f'| {q} | {len(R)} | ' + ' | '.join(f'{v:.4f}' for v in tm) + ' | ' + ' | '.join(f'{v:.4f}' for v in tl) + ' | ' + ' | '.join(f(sq[k] / sq[0], 3) for k in (2, 3, 4)) + ' |')
print()
print('Extreme-value collapse test from the main runs: -ln P(max_M <= t)/(3q) (q-independent if the used deficits behave like 3q iid draws from one rescaled law)')
print('| q | t=3 | 3.5 | 4 | 4.5 | 5 |'); print('|---|---|---|---|---|---|')
for q, R in M.items():
    a = np.array([r['maxdM'] for r in R]); row = []
    for t in (3, 3.5, 4, 4.5, 5):
        pr = np.mean(a <= t); row.append(f'{-np.log(pr) / (3 * q):.4f}' if 0 < pr else 'inf')
    print(f'| {q} | ' + ' | '.join(row) + ' |')
print()
print('## 3b. LP weight by pairs shared with M, and Berge triangles in the fractional support (D > 0 instances)')
print('| q | n | LP weight/q on teams sharing 0 pairs with M | 1 pair | 3 pairs (M teams) | n (D>0) | triangles per instance | strong triangles | share of strong triangles whose centre (a,b,c) is an M team | centre in LP support |')
print('|---|---|---|---|---|---|---|---|---|---|')
for q, R in E.items():
    wk = np.array([r['w_keep'] for r in R if not r['integral']]).mean(0) / q
    F = [r for r in R if not r['integral']]
    g = lambda k: np.array([r[k] for r in F]).sum()
    print(f"| {q} | {len(R)} | {f(wk[0])} | {f(wk[1])} | {f(wk[3])} | {len(F)} | {f(g('ntri')/len(F),2)} | {f(g('ntri_strong')/len(F),2)} | {f(g('centreM')/max(g('ntri_strong'),1),3)} | {f(g('centreSup')/max(g('ntri_strong'),1),3)} |")
