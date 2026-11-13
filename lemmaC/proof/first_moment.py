"""Exact finite-q constants of Theorem 1 (first-moment lower bound on E(3q - T)) and Theorem 2
(union-bound dual lower bound on E(3q - LP*)), against the measured means in structure/STRUCTURE.md."""
from math import lgamma, exp, log
data = {8: (4.89, 4.88), 12: (5.89, 5.87), 16: (6.44, 6.40), 24: (7.48, 7.41), 32: (8.21, 8.10), 40: (8.90, 8.77), 48: (9.35, 9.19)}
print('q   t*_q  EIPlb  EIPlb/q^(1/3)  ELPlb/q^(1/3)  data(3q-T)/q^(1/3)  data(3q-LP*)/q^(1/3)')
for q in list(data) + [10**3, 10**5, 10**7]:
    tstar = exp((lgamma(3*q+1) - 2*lgamma(q+1)) / (3*q))     # (q!)^2 t^{3q}/(3q)! = 1 at t = t*
    ip = tstar * 3*q / (3*q + 1)                              # = int_0^t* (1 - (t/t*)^{3q}) dt
    lp = q * 0.75 * (6 / q**2) ** (1/3)                       # q * int_0^inf (1 - min(1, q^2 t^3/6)) dt
    s = q ** (1/3); d = data.get(q)
    print(q, round(tstar, 3), round(ip, 3), round(ip/s, 4), round(lp/s, 4), *(round(x/s, 3) for x in d) if d else '')
print('limit beta = 3 e^(-1/3) =', 3*exp(-1/3), '  LP constant (3/4) 6^(1/3) =', 0.75*6**(1/3))
