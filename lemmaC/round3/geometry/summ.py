import json, glob, numpy as np
print('== families ==')
for fn in sorted(glob.glob('fam_q*.jsonl')):
    R = [json.loads(l) for l in open(fn)]; q = R[0]['q']
    D = np.array([r['D'] for r in R]); D21 = np.array([r['D21'] for r in R]); D111 = np.array([r['D111'] for r in R]); pos = D > 1e-7
    f = lambda key: np.median([r[key] for r in R if r['D'] > 1e-7] or [0])
    print(f"q={q:2d} n={len(R):3d} P(D>0)={pos.mean():.3f} P(D21>0)={(D21 > 1e-7).mean():.3f} P(D111>0)={(D111 > 1e-7).mean():.3f} "
          f"ED={D.mean():.4f} ED21={D21.mean():.4f} ED111={D111.mean():.4f} | steepest ray: k/q={f('k') / q:.2f} supp/k={np.median([r['supp'] / r['k'] for r in R if r['D'] > 1e-7] or [0]):.2f} "
          f"zmax/zmin={f('zratio'):.1f} pair-keeping weight share={np.mean([r['keep_frac'] for r in R if r['D'] > 1e-7] or [0]):.2f}")
print('== clone (hub) and half-integral ==')
for fn in sorted(glob.glob('clone_q*.jsonl')):
    R = [json.loads(l) for l in open(fn)]; q = R[0]['q']
    D = np.array([r['D'] for r in R]); gh = np.array([r['ghub'] for r in R]); Dh = np.array([r.get('Dhalf', 0) for r in R]); pos = D > 1e-7
    nh = np.array([r['nhub'] for r in R]); ks = [k for r in R for k in r['khub']]
    print(f"q={q:2d} n={len(R):3d} P(D>0)={pos.mean():.3f} P(hub>0)={(gh > 1e-7).mean():.3f} P(half>0)={(Dh > 1e-7).mean():.3f} "
          f"P(hub>0|D>0)={(gh[pos] > 1e-7).mean() if pos.any() else 0:.3f} P(half>0|D>0)={(Dh[pos] > 1e-7).mean() if pos.any() else 0:.3f} "
          f"ED={D.mean():.4f} E[ghub/2]={gh.mean() / 2:.4f} E[Dhalf]={Dh.mean():.4f} improving hubs per instance (given D>0)={nh[pos].mean() if pos.any() else 0:.2f} of {q} "
          f"median M-teams changed by improving hub structure={np.median(ks) if ks else 0}")
print('== kinds of pair-keeping teams added to the pair-disjoint teams: (P(gap>0), E gap) ==')
for fn in sorted(glob.glob('kinds_q*.json')):
    R = json.load(open(fn))
    print(f"q={R[0]['q']:2d} n={len(R)} " + ' '.join(f"{k}: {np.mean([r[k] > 1e-7 for r in R]):.3f}, {np.mean([r[k] for r in R]):.4f} |" for k in ('D', 'AB', 'AB+BC', 'AB+AC', 'BC+AC')))
print('== thinning: keep each pair-keeping team with probability phi ==')
for fn in sorted(glob.glob('thin_q*.jsonl')):
    R = [json.loads(l) for l in open(fn)]
    print(f"q={R[0]['q']:2d} n={len(R)} " + ' '.join(f"{k}: P={np.mean([r[k] > 1e-7 for r in R]):.3f} E={np.mean([r[k] for r in R]):.4f} |" for k in ('D', 'phi0.5', 'phi0.25', 'phi0.125')))
print('== column generation (big.py): D, D111 (pair-disjoint only), D21 (pair-keeping only) ==')
for fn in sorted(glob.glob('big_q*.jsonl')):
    R = [json.loads(l) for l in open(fn)]
    D = np.array([r['D'] for r in R]); a = np.array([r['D111'] for r in R]); b = np.array([r['D21'] for r in R])
    print(f"q={R[0]['q']:2d} n={len(R)} P(D>0)={(D > 1e-7).mean():.3f} ED={D.mean():.4f} | D111: n>0={(a > 1e-7).sum()} E={a.mean():.5f} values={np.round(a[a > 1e-7], 4).tolist()} | "
          f"D21: n>0={(b > 1e-7).sum()} E={b.mean():.5f} | E[D111]/E[D]={a.mean() / D.mean():.4f}")
print('== control: independent team values, same marginal law (iid.py) ==')
for fn in sorted(glob.glob('iid_q*.jsonl')):
    R = [json.loads(l) for l in open(fn)]; D = np.array([r['D'] for r in R]); q = R[0]['q']
    print(f"q={q:2d} n={len(R)} P(D>0)={(D > 1e-7).mean():.3f} ED={D.mean():.4f} sd={D.std():.4f} ED/q^(1/3)={D.mean() / q ** (1 / 3):.4f} E(3q-T)/q^(1/3)={np.mean([r['defT'] for r in R]) / q ** (1 / 3):.3f}")
print('== interpolation (hybrid.py): number of pair-additive legs 0..3 ==')
for fn in sorted(glob.glob('hybrid_q*.json')):
    H = json.load(open(fn)); q = H['q']
    for k in '0123':
        D = np.array(H['D'][k]); print(f"q={q:2d} n={H['n']} legs={k}: P(D>0)={(D > 1e-7).mean():.3f} ED={D.mean():.4f} E(3q-T)/q^(1/3)={np.mean(H['defT'][k]) / q ** (1 / 3):.3f}")
