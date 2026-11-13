import json, glob, numpy as np
print('q n | mean mu  var mu  4E[Q]  ratio | var T  4E[R_T]  E(3q-T) | E[1/d]  max c1  max cinf | P(mu>0)  4E[Q]/mean^2')
for fn in sorted(glob.glob('out/poin_q*.jsonl')):
    R = [json.loads(l) for l in open(fn)]
    if len(R) < 15: continue
    g = lambda k: np.array([r[k] for r in R])
    mu, Q, T, RT, d = g('mu'), g('Q'), g('T'), g('R_T'), g('d')
    print(f"{R[0]['q']:2d} {len(R):3d} | {mu.mean():+.4f} {mu.var(ddof=1):.2e} {4*Q.mean():.2e} {4*Q.mean()/mu.var(ddof=1):6.1f} | "
          f"{T.var(ddof=1):.3f} {4*RT.mean():.3f} {g('def_T').mean():.2f} | {np.mean(1/d):.3f} {g('c1').max():.3f} {g('cinf').max():.3f} | "
          f"{(mu>0).mean():.2f} {4*Q.mean()/mu.mean()**2:.2f}")
