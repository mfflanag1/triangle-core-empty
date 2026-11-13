import json, glob, numpy as np
from collections import defaultdict
def load(pat):
    out = []
    for f in glob.glob(pat):
        for l in open(f):
            try: out.append(json.loads(l))
            except Exception: pass
    return out
lpo = load('lp.jsonl'); ipr = load('ip*.jsonl')
def ms(v):
    v = np.asarray(v, float); return v.mean(), (v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else float('nan'))
# LP table: union over seeds (dedupe by q,K,seed)
lp = {}
for r in lpo + ipr:
    if r.get('feas'): lp[(r['q'], r['K'], r['seed'])] = r
infeas = sum(1 for r in lpo + ipr if not r.get('feas'))
print('infeasible instances:', infeas)
byq = defaultdict(list)
for (q, K, s), r in lp.items(): byq[(K, q)].append(r)
print('LP  K   q     n   mean    se     sd   P(int)  max-used  secs')
LP = {}
for (K, q) in sorted(byq):
    rs = byq[(K, q)]; m, se = ms([r['lp'] for r in rs]); LP[(K, q)] = (m, se, len(rs))
    print(f"   {K:.0f} {q:5d} {len(rs):5d} {m:.4f} {se:.4f} {np.std([r['lp'] for r in rs]):.3f} {np.mean([r['lpint'] for r in rs]):.4f}  {max(r['lpmax'] for r in rs):.2f}  {np.mean([r['tlp'] for r in rs]):.2f}")
ib = defaultdict(list)
for r in ipr:
    if r.get('feas'): ib[(r['K'], r['q'])].append(r)
print('IP  K   q    n closed  lp(sub) ub      lb      gap_ub  se     gap_lb  P(D=0) IPcv_ub IPcv_lb max-used secs')
G = {}
for (K, q) in sorted(ib):
    rs = [r for r in ib[(K, q)] if r.get('ub') is not None]; nn = len(ib[(K, q)]) - len(rs)
    gu = [r['ub'] - r['lp'] for r in rs]; gl = [r['lb'] - r['lp'] for r in rs]
    mgu, seu = ms(gu); mgl, _ = ms(gl); mu, seU = ms([r['ub'] for r in rs]); ml, _ = ms([r['lb'] for r in rs])
    G[(K, q)] = (mgu, seu, mgl, mu, seU)
    print(f"   {K:.0f} {q:4d} {len(rs):4d} {np.mean([r['closed'] for r in rs]):.2f}   {np.mean([r['lp'] for r in rs]):.4f} {mu:.4f} {ml:.4f} {mgu:.4f} {seu:.4f} {mgl:.4f} {np.mean([g < 1e-7 for g in gu]):.3f}  {LP[(K,q)][0]+mgu:.4f} {LP[(K,q)][0]+mgl:.4f} {max(r['ipmax'] for r in rs):.2f} {np.mean([r['tip'] for r in rs]):.1f} nosol={nn} se_ub={seU:.4f}")
# coupled K=6 vs 8
for q in sorted({q for (K, q) in byq if K == 8}):
    d = [lp[(q, 8.0, s)]['lp'] - lp[(q, 6.0, s)]['lp'] for (qq, K, s) in lp if qq == q and K == 8 and (q, 6.0, s) in lp]
    if d: print(f'LP(K=8)-LP(K=6) coupled q={q}: n={len(d)} mean={np.mean(d):.2e} min={np.min(d):.2e} frac nonzero={np.mean(np.abs(d) > 1e-9):.3f}')
i6 = {(r['q'], r['seed']): r for r in ipr if r.get('feas') and r['K'] == 6 and r.get('closed')}
for q in sorted({r['q'] for r in ipr if r['K'] == 8}):
    d = [r['ub'] - i6[(q, r['seed'])]['ub'] for r in ipr if r['K'] == 8 and r['q'] == q and r.get('closed') and (q, r['seed']) in i6]
    if d: print(f'IP(K=8)-IP(K=6) coupled q={q}: n={len(d)} mean={np.mean(d):.2e} min={np.min(d):.2e} frac nonzero={np.mean(np.abs(np.array(d)) > 1e-7):.3f}')
def fit(qs, y, se, p, label):
    qs = np.array(qs, float); y = np.array(y); se = np.array(se)
    X = np.c_[np.ones(len(qs)), qs ** (-p)]; W = 1 / se ** 2
    C = np.linalg.inv(X.T @ (W[:, None] * X)); beta = C @ X.T @ (W * y)
    chi = float(((y - X @ beta) ** 2 * W).sum())
    print(f'  {label}: a + b q^-{p:.3f}: a={beta[0]:.4f} +- {np.sqrt(C[0,0]):.4f}, b={beta[1]:.3f}, chi2={chi:.1f}/{len(qs)-2}')
for qmin in (16, 32, 64):
    qs = [q for (K, q) in sorted(LP) if K == 6 and q >= qmin]
    if len(qs) >= 3:
        print(f'LP fits (K=6, q>={qmin}):')
        for p in (1/3, 2/3, 1.0): fit(qs, [LP[(6.0, q)][0] for q in qs], [LP[(6.0, q)][1] for q in qs], p, 'LP')
qs = [q for (K, q) in sorted(G) if K == 6 and all(r.get('closed') for r in ib[(K, q)])]
if len(qs) >= 3:
    print('gap (IP-LP) fits, closed sizes', qs)
    for p in (1/3, 2/3, 1.0): fit(qs, [G[(6.0, q)][0] for q in qs], [G[(6.0, q)][1] for q in qs], p, 'gap')
    print('IP (= LP mean over all seeds + mean gap) fits')
    for p in (1/3, 2/3, 1.0): fit(qs, [LP[(6.0, q)][0] + G[(6.0, q)][0] for q in qs], [np.hypot(LP[(6.0, q)][1], G[(6.0, q)][1]) for q in qs], p, 'IP')
