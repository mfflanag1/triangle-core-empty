import json, glob, numpy as np
EPS = 1e-7
se = lambda a: np.std(a) / len(a) ** .5
print('q | n | P(D>0) | E[N_kill;D>0] | E[N_create;D=0] | N_kill given D>0 | N_create/3q given D=0 | P(D^-i>0) | sum rho_IP vs 3q-T | sum rho_LP vs 3q-LP | drift at D=0')
for fn in sorted(glob.glob('balance_q*.jsonl')):
    R = [json.loads(l) for l in open(fn)]; q = R[0]['q']; n = len(R)
    D = np.array([r['D'] for r in R]); pos = D > EPS
    kill = np.array([sum(1 for (m, a, b) in r['res'] if (a - b) + r['D'] <= EPS) if r['D'] > EPS else 0 for r in R], float)
    crea = np.array([sum(1 for (m, a, b) in r['res'] if m and (a - b) > EPS) if r['D'] <= EPS else 0 for r in R], float)
    sIP = np.array([sum(a for (m, a, b) in r['res']) for r in R]); sLP = np.array([sum(b for (m, a, b) in r['res']) for r in R])
    T = np.array([r['T'] for r in R]); LP = np.array([r['LP'] for r in R])
    anyD = np.array([np.mean([(a - b) + r['D'] > EPS for (m, a, b) in r['res'] if m]) for r in R])
    drift = sIP - sLP
    z = ~pos
    print(f'{q} | {n} | {pos.mean():.3f}+-{se(pos):.3f} | {kill.mean():.2f}+-{se(kill):.2f} | {crea.mean():.2f}+-{se(crea):.2f} | {kill[pos].mean():.1f} | '
          f'{(crea[z].mean()/(3*q) if z.any() else float("nan")):.2f} (n={z.sum()}) | {anyD.mean():.3f}+-{se(anyD):.3f} | {sIP.mean():.2f}+-{se(sIP):.2f} vs {(3*q-T).mean():.2f} | '
          f'{sLP.mean():.2f}+-{se(sLP):.2f} vs {(3*q-LP).mean():.2f} | {(drift[z].mean() if z.any() else float("nan")):.2f}')
