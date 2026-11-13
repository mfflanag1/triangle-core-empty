import json, glob, numpy as np
M = {}
for fn in sorted(glob.glob('out/main_q*.jsonl')):
    R = [json.loads(l) for l in open(fn)]
    if R: M[R[0]['q']] = R
qs = np.array(sorted(M)); lq = np.log(qs)
for k in ('maxdM', 'maxdL'):
    m = np.array([np.mean([r[k] for r in M[q]]) for q in qs]); se = np.array([np.std([r[k] for r in M[q]]) / np.sqrt(len(M[q])) for q in qs])
    b, a = np.polyfit(lq, m, 1, w=1 / se); pb, pa = np.polyfit(lq, np.log(m), 1, w=m / se)
    b3, a3 = np.polyfit(np.log(3 * qs) ** (2 / 3), m, 1, w=1 / se)
    print(k, 'means', m.round(3), 'se', se.round(3), f'| log fit {a:.2f} + {b:.2f} ln q, resid {np.abs(m - a - b * lq).max():.3f} | power fit exponent {pb:.3f}, resid {np.abs(m - np.exp(pa) * qs ** pb).max():.3f} | (ln 3q)^(2/3) fit {a3:.2f} + {b3:.2f} x, resid {np.abs(m - a3 - b3 * np.log(3 * qs) ** (2 / 3)).max():.3f}')
print('q', qs)
sd = np.array([np.std([r['T'] for r in M[q]], ddof=1) for q in qs]); n = np.array([len(M[q]) for q in qs])
print('sd T', sd.round(3), 'x q^(1/6)', (sd * qs ** (1 / 6)).round(3), 'se of sd', (sd / np.sqrt(2 * n)).round(3), 'loglog slope', np.polyfit(lq, np.log(sd), 1, w=np.sqrt(2 * n))[0].round(3))
S2 = np.array([np.mean([r['sumsqM'] for r in M[q]]) for q in qs]); print('E S2(M) q^(1/3)', (S2 * qs ** (1 / 3)).round(3), 'loglog slope', np.polyfit(lq, np.log(S2), 1)[0].round(3))
g = np.array([np.mean([r['gradD2'] for r in M[q]]) for q in qs]); print('E gradD2', g.round(3))
L2 = np.array([np.mean([r['sumsqL2'] for r in M[q]]) for q in qs]); print('E sum load^2 d^2 x q^(1/3)', (L2 * qs ** (1 / 3)).round(3))
ED = np.array([np.mean([r['D'] for r in M[q]]) for q in qs]); print('E D', ED.round(4), 'se', np.array([np.std([r['D'] for r in M[q]]) / np.sqrt(len(M[q])) for q in qs]).round(4), 'ED/q^(1/3)', (ED / qs ** (1 / 3)).round(4))
b, a = np.polyfit(qs ** (1 / 3), ED, 1); print(f'E D fit {b:.3f} q^(1/3) + {a:.3f}')
print('corr(D, defT)', [np.corrcoef([r['D'] for r in M[q]], [r['defT'] for r in M[q]])[0, 1].round(2) for q in qs])
S = {}
for fn in sorted(glob.glob('out/sister_*.jsonl')):
    R = [json.loads(l) for l in open(fn)]; R = list({r['idx']: r for r in R[::-1]}.values())
    S[(R[0]['model'], R[0]['q'])] = np.mean([r['integral'] for r in R])
for m in ('P', 'IH', 'U'):
    q_ = np.array(sorted(q for (mm, q) in S if mm == m)); p = np.array([S[(m, q)] for q in q_])
    print(m, 'q', q_, '-ln P(int)', (-np.log(p)).round(2), 'slope per unit q over q>=10:', np.polyfit(q_[q_ >= 10], -np.log(p[q_ >= 10]), 1)[0].round(3))
