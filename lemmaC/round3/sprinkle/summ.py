import json, glob, numpy as np
for fn in sorted(glob.glob('regret_q*.jsonl')):
    R = [json.loads(l) for l in open(fn)]
    q = R[0]['q']; n = len(R)
    D = np.array([r['D'] for r in R]); T = np.array([r['T'] for r in R])
    pos = D > 1e-7
    # per-instance mean regret over sampled M pairs, scaled to the sum over all 3q M-pairs
    sIP = np.array([3*q*np.mean([x['rIP'] for x in r['del_']]) for r in R])
    sLPM = np.array([3*q*np.mean([x['rLP'] for x in r['del_']]) for r in R])
    dD = np.array([3*q*np.mean([x['D2'] for x in r['del_']]) for r in R]) - 3*q*D
    se = lambda a: a.std()/len(a)**.5 if len(a) > 1 else float('nan')
    print(f'q={q} n={n} P(D>0)={pos.mean():.3f} ED={D.mean():.4f}')
    print(f'   E sum_M rho_IP = {sIP.mean():.3f}+-{se(sIP):.3f}   E(3q-T) = {(3*q-T).mean():.3f}+-{se(3*q-T):.3f}   paired diff {(sIP-(3*q-T)).mean():.3f}+-{se(sIP-(3*q-T)):.3f}')
    print(f'   E sum_M rho_LP = {sLPM.mean():.3f}+-{se(sLPM):.3f}  (M pairs only)   E sum_M (D^-i - D) = {dD.mean():.3f}+-{se(dD):.3f}')
    z = [x for r in R if r['D'] <= 1e-7 for x in r['del_']]
    p = [x for r, pp in zip(R, pos) if pp for x in r['del_']]
    if z:
        d2 = np.array([x['D2'] for x in z]); print(f'   on D=0 ({(~pos).sum()} inst, {len(z)} deletions): P(D^-i>0) = {(d2>1e-7).mean():.3f}, E D^-i = {d2.mean():.4f}, E sum_M D^-i = {3*q*d2.mean():.3f}')
    if p:
        d2 = np.array([x['D2'] for x in p]); Dp = np.array([r['D'] for r, pp in zip(R, pos) if pp for x in r['del_']])
        print(f'   on D>0 ({pos.sum()} inst, {len(p)} deletions): P(D^-i=0) = {(d2<=1e-7).mean():.3f}, P(D^-i>D) = {(d2>Dp+1e-9).mean():.3f}, P(D^-i<D) = {(d2<Dp-1e-9).mean():.3f}, E(D^-i-D) = {(d2-Dp).mean():.4f}, x3q = {3*q*(d2-Dp).mean():.3f}')
