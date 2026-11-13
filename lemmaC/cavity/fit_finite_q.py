"""Fit the finite-q means (STRUCTURE.md section 1, q = 12..48) to a q^(1/3) + b + c q^(-1/3),
weighted by standard errors sd/sqrt(n); also with a fixed to the cavity values."""
import numpy as np
q = np.array([12, 16, 24, 32, 40, 48.]); n = np.array([200, 200, 200, 200, 200, 120.])
T = np.array([5.89, 6.44, 7.48, 8.21, 8.90, 9.35]); sdT = np.array([.58, .52, .48, .46, .46, .43])
L = np.array([5.87, 6.40, 7.41, 8.10, 8.77, 9.19])          # sd(3q-LP) ~ sd(3q-T)
G = np.array([.022, .039, .069, .108, .127, .156]); sdG = np.array([.036, .046, .053, .054, .053, .054])
def wfit(cols, y, se):
    A = np.vstack(cols).T / se[:, None]; b = y / se
    c, *_ = np.linalg.lstsq(A, b, rcond=None); cov = np.linalg.inv(A.T @ A)
    chi2 = np.sum((A @ c - b)**2)
    return c, np.sqrt(np.diag(cov)), chi2
x = q**(1/3)
for name, y, se in [('3q-T', T, sdT / np.sqrt(n)), ('3q-LP*', L, sdT / np.sqrt(n)), ('LP*-T', G, sdG / np.sqrt(n))]:
    for lab, cols in [('a q^1/3 + b', (x, np.ones_like(q))), ('a q^1/3 + b + c q^-1/3', (x, np.ones_like(q), 1 / x))]:
        c, s, chi2 = wfit(cols, y, se)
        print(f"{name:7s} {lab:24s} coef={np.round(c,3)} se={np.round(s,3)} chi2={chi2:.1f} (dof {len(q)-len(c)})")
for name, y, a in [('3q-T', T, 2.2495), ('3q-LP*', L, 2.197), ('LP*-T', G, 0.052)]:
    se = (sdG if name == 'LP*-T' else sdT) / np.sqrt(n)
    c, s, chi2 = wfit((np.ones_like(q), 1 / x), y - a * x, se)
    print(f"{name:7s} a fixed={a}: b={c[0]:.3f}+-{s[0]:.3f} c={c[1]:.3f}+-{s[1]:.3f} chi2={chi2:.1f} (dof 4)")
