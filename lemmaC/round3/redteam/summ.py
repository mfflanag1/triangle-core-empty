import json, numpy as np, collections
d = collections.defaultdict(list)
for l in open('lpfrac.jsonl'):
    r = json.loads(l); d[r['q']].append(r)
print('q  n  LPdef/q^(1/3) (se)  ffrac (se)  P(frac)  ffrac|frac  half-share  xmin_med  secs')
for q in sorted(d):
    R = d[q]; a = np.array([r['LPdef'] for r in R]) / q ** (1 / 3); f = np.array([r['ffrac'] for r in R])
    hs = np.mean([r['half'] for r in R if r['ffrac'] > 0]) if (f > 0).any() else float('nan')
    print(q, len(R), round(a.mean(), 4), round(a.std(ddof=1) / len(R) ** .5, 4), round(f.mean(), 3), round(f.std(ddof=1) / len(R) ** .5, 3),
          round((f > 0).mean(), 2), round(f[f > 0].mean(), 3) if (f > 0).any() else None, round(hs, 2),
          round(float(np.median([r['xmin'] for r in R])), 4), round(np.mean([r['secs'] for r in R]), 1))
