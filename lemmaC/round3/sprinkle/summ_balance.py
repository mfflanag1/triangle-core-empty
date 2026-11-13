import json, glob, numpy as np
EPS = 1e-7
se = lambda a: np.std(a) / len(a) ** .5
for fn in sorted(glob.glob('balance_q*.jsonl')):
    R = [json.loads(l) for l in open(fn)]; q = R[0]['q']; n = len(R)
    D = np.array([r['D'] for r in R]); pos = D > EPS
    kill = np.array([sum(1 for (m, a, b) in r['res'] if (a - b) + r['D'] <= EPS) if r['D'] > EPS else 0 for r in R], float)
    crea = np.array([sum(1 for (m, a, b) in r['res'] if m and (a - b) > EPS) if r['D'] <= EPS else 0 for r in R], float)
    tot = np.array([sum(a - b for (m, a, b) in r['res']) for r in R])          # sum_i (D^{-i} - D)
    sIP = np.array([sum(a for (m, a, b) in r['res']) for r in R]); sLP = np.array([sum(b for (m, a, b) in r['res']) for r in R])
    T = np.array([r['T'] for r in R]); LP = np.array([r['LP'] for r in R])
    anyD = np.array([np.mean([(a - b) + r['D'] > EPS for (m, a, b) in r['res'] if m]) for r in R])   # P(D^{-i}>0), i uniform in M
    print(f'q={q} n={n}  P(D>0)={pos.mean():.3f}+-{se(pos):.3f}  E D={D.mean():.4f}+-{se(D):.4f}')
    print(f'   E[N_kill;D>0]={kill.mean():.3f}+-{se(kill):.3f}   E[N_create;D=0]={crea.mean():.3f}+-{se(crea):.3f}   paired diff {(kill-crea).mean():.3f}+-{se(kill-crea):.3f}')
    print(f'   E sum_i(D^-i - D)={tot.mean():.4f}+-{se(tot):.4f}  vs E D   | E sum rho_IP={sIP.mean():.3f}+-{se(sIP):.3f} vs E(3q-T)={(3*q-T).mean():.3f} | E sum rho_LP={sLP.mean():.3f}+-{se(sLP):.3f} vs E(3q-LP)={(3*q-LP).mean():.3f}')
    print(f'   P(D^-i>0), i unif in M: {anyD.mean():.3f}+-{se(anyD):.3f}  => bound P(D>0) >= {anyD.mean()/5:.3f};  mean #pairs in M u supp: {np.mean([r["nU"] for r in R]):.1f} (max {max(r["nU"] for r in R)}, 12q-6={12*q-6})')
    print(f'   on D>0: mean N_kill={kill[pos].mean():.2f}; on D=0: mean N_create={crea[~pos].mean():.2f} of {3*q}')
