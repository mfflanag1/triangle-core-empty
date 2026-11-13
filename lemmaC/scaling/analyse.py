"""Tables and fits for SCALING.md. q<=48 from ../structure/out/struct_q*.jsonl (gate-1 seeds, full LP);
q>=56 from out/scal_q*.jsonl (tag-7 seeds, column generation). Floating point."""
import json, glob, os, re
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
data = {}
for fn in glob.glob(os.path.join(HERE, '..', 'structure', 'out', 'struct_q*.jsonl')):
    q = int(re.findall(r'q(\d+)', fn)[-1])
    if q >= 16: data[q] = [json.loads(l) for l in open(fn)]
for fn in glob.glob(os.path.join(HERE, 'out', 'scal_q*.jsonl')):
    q = int(re.findall(r'q(\d+)', fn)[-1]); data[q] = [json.loads(l) for l in open(fn)]
qs = sorted(q for q in data if len(data[q]) >= 3)
arr = {q: np.array([[r['gap'], 3 * q - r['T'], 3 * q - r['LP']] for r in data[q]]) for q in qs}
print('| q | n | E D | sd D | se E D | P(D>0) | E D / q^(1/3) | E D / q^(1/4) | E D / q^(1/6) | E(3q-T) | /q^(1/3) | E(3q-LP*) | ratio (3q-LP*)/(3q-T) |')
print('|' + '---|' * 13)
for q in qs:
    D, a, b = arr[q].T; n = len(D)
    print(f'| {q} | {n} | {D.mean():.4f} | {D.std(ddof=1):.4f} | {D.std(ddof=1)/np.sqrt(n):.4f} | {(D>1e-7).mean():.3f} | '
          f'{D.mean()/q**(1/3):.4f} | {D.mean()/q**.25:.4f} | {D.mean()/q**(1/6):.4f} | {a.mean():.3f} | {a.mean()/q**(1/3):.3f} | '
          f'{b.mean():.3f} | {b.mean()/a.mean():.4f} |')
rng = np.random.default_rng(1)
def fit(sub, col=0, B=4000, f=lambda x: x):
    lq = np.log(sub)
    def slope(ms): return np.polyfit(lq, np.log(ms), 1)[0]
    est = slope([f(arr[q][:, col]).mean() for q in sub])
    bs = [slope([f(arr[q][rng.integers(0, len(arr[q]), len(arr[q])), col]).mean() for q in sub]) for _ in range(B)]
    return est, np.percentile(bs, [5, 95]), np.percentile(bs, [2.5, 97.5])
print()
for lo in (16, 24, 32, 40):
    sub = [q for q in qs if q >= lo]
    if len(sub) < 3: continue
    e, ci90, ci95 = fit(sub)
    print(f'E D ~ q^b over q in [{lo},{sub[-1]}]: b = {e:.3f}, 90% [{ci90[0]:.3f}, {ci90[1]:.3f}], 95% [{ci95[0]:.3f}, {ci95[1]:.3f}]')
print('local slopes of E D:', ', '.join(f'{a}->{b}: {np.log(arr[b][:,0].mean()/arr[a][:,0].mean())/np.log(b/a):.2f}' for a, b in zip(qs, qs[1:])))
for col, nm in ((1, '3q-T'), (2, '3q-LP*')):
    sub = [q for q in qs if q >= 16]; e, ci90, _ = fit(sub, col)
    print(f'{nm} ~ q^b over [16,{sub[-1]}]: b = {e:.3f}, 90% [{ci90[0]:.3f}, {ci90[1]:.3f}]')
# ratio trend: 1 - ratio = E D / E(3q-T)
print('1 - ratio = E D/E(3q-T):', ', '.join(f'{q}: {arr[q][:,0].mean()/arr[q][:,1].mean():.4f}' for q in qs))
# saturation model: E D = c q^(1/3) - d (offset) fit vs pure power, weighted least squares
x = np.array(qs, float); y = np.array([arr[q][:, 0].mean() for q in qs]); se = np.array([arr[q][:, 0].std(ddof=1) / np.sqrt(len(arr[q])) for q in qs])
for nm, X in (('c q^(1/3) + d', np.c_[x ** (1/3), np.ones_like(x)]), ('c log q + d', np.c_[np.log(x), np.ones_like(x)]),
              ('c q^(1/4) + d', np.c_[x ** .25, np.ones_like(x)])):
    W = 1 / se; coef, res, *_ = np.linalg.lstsq(X * W[:, None], y * W, rcond=None)
    chi2 = float(((X @ coef - y) / se) ** 2 @ np.ones_like(y))
    print(f'WLS {nm}: coef {coef.round(4)}, chi2 {chi2:.1f} on {len(y)-2} dof')
e = np.polyfit(np.log(x), np.log(y), 1, w=y / se)
print(f'WLS pure power on all q>=16: b = {e[0]:.3f}, chi2', round(float((((np.exp(np.polyval(e, np.log(x))) - y) / se) ** 2).sum()), 1))
