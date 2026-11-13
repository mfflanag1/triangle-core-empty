"""Tail of the rescaled team cost in the optimal assignment M and in the LP optimum x*.

Floating point: evidence only. Uses ../../scaling/cg.py (column-generation LP, reduced-cost-fixing MILP).
Reports, per q: mean cost per team (1/q) sum xi_e w_e, and the share carried by teams with xi_e > K.
"""
import os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scaling'))
import cg
from scipy.optimize import milp, LinearConstraint, Bounds

KS = (2, 4, 6, 8, 10)


def inst(q, idx):
    A, B, C = cg.draw(q, idx, 31)
    V = A[:, :, None] + B[None, :, :] + C[:, None, :]
    vf = V.reshape(-1)
    xi = (3 - vf) * q ** (2 / 3)
    LPs, p, r, cols, x, it = cg.solve_lp(V, q)
    lp_tail = [float((xi[cols] * x * (xi[cols] > K)).sum() / q) for K in KS]
    lp_mean = float((xi[cols] * x).sum() / q)
    lp_max = float(xi[cols][x > 1e-9].max())
    delta = 0.05
    while True:
        cand = np.nonzero(r >= -delta)[0]
        mi = milp(-vf[cand], constraints=LinearConstraint(cg.amat(cand, q), np.ones(3 * q), np.ones(3 * q)),
                  integrality=np.ones(len(cand)), bounds=Bounds(0, 1), options=dict(time_limit=600, mip_rel_gap=1e-12))
        if mi.status == 0 and -mi.fun >= LPs - delta + 1e-9:
            break
        delta *= 2
    m = xi[cand][mi.x > 0.5]
    assert len(m) == q
    return lp_mean, lp_max, lp_tail, float(m.mean()), float(m.max()), [float((m * (m > K)).sum() / q) for K in KS], float((m ** 2).mean())


def main():
    print("q reps | IP: mean cost/team, mean max xi, mean xi^2, tail share E[xi;xi>K]/q for K=2,4,6,8,10 | LP: mean, mean max, tails")
    for q, reps in ((16, 60), (24, 60), (32, 40), (48, 16)):
        t0 = time.time(); rows = [inst(q, i) for i in range(reps)]
        ipm = np.mean([r[3] for r in rows]); ipx = np.mean([r[4] for r in rows]); ip2 = np.mean([r[6] for r in rows])
        ipt = np.mean([r[5] for r in rows], axis=0)
        lpm = np.mean([r[0] for r in rows]); lpx = np.mean([r[1] for r in rows]); lpt = np.mean([r[2] for r in rows], axis=0)
        print(f"{q} {reps} | IP {ipm:.3f} {ipx:.2f} {ip2:.2f} " + " ".join(f"{v:.4f}" for v in ipt) +
              f" | LP {lpm:.3f} {lpx:.2f} " + " ".join(f"{v:.4f}" for v in lpt) + f" | {time.time() - t0:.0f}s")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
