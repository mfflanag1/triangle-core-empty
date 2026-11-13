"""Weighted linear fits of the brief's E D table to D = a*g(q) - b. se = 0.05/sqrt(n). Evidence only."""
import numpy as np
q = np.array([16, 24, 32, 40, 48, 56, 64.]); n = np.array([200, 200, 200, 200, 120, 70, 23.])
D = np.array([0.039, 0.069, 0.108, 0.127, 0.156, 0.188, 0.212]); se = 0.05 / np.sqrt(n)
G = {'q^(1/2)': lambda x: x ** .5, 'q^(1/3)': lambda x: x ** (1 / 3), 'q^(1/3)/log q': lambda x: x ** (1 / 3) / np.log(x),
     'q^(1/4)': lambda x: x ** .25, 'q^(1/6)': lambda x: x ** (1 / 6), 'log q': np.log, 'q^(1/3)/log^2 q': lambda x: x ** (1 / 3) / np.log(x) ** 2}
for nm, g in G.items():
    X = np.stack([g(q), -np.ones(7)], 1) / se[:, None]
    p, res, *_ = np.linalg.lstsq(X, D / se, rcond=None)
    chi2 = float((((X @ p) - D / se) ** 2).sum())
    pr = [round(float(p[0] * g(np.array([float(Q)]))[0] - p[1]), 3) for Q in (96, 128, 200, 400)]
    print(f'D = a*{nm:16s} - b: a={p[0]:8.4f} b={p[1]:7.4f} chi2={chi2:6.2f} (dof 5)  pred D(96,128,200,400)={pr}')
