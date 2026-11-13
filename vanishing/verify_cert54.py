"""Exact re-check (Fractions) of the q=3 certificate idx=54 used in Proposition 3."""
from fractions import Fraction as F
import itertools
A=[["33/100","3/20","17/25"],["1/25","43/50","12/25"],["7/20","3/50","31/100"]]
B=[["9/100","19/20","87/100"],["1/25","3/20","33/50"],["83/100","12/25","22/25"]]
C=[["7/50","77/100","59/100"],["13/25","67/100","9/100"],["71/100","1/2","93/100"]]
A,B,C=[[[F(x) for x in r] for r in M] for M in (A,B,C)]
v=lambda a,b,c: A[a][b]+B[a][c]+C[b][c]
lam={(0,0,1):F(1,2),(0,2,2):F(1,2),(1,1,0):F(1,2),(1,1,1):F(1,2),(2,0,2):F(1,2),(2,2,0):F(1,2)}
for s in range(3):
    for i in range(3): assert sum(w for t,w in lam.items() if t[s]==i)==1
LP=sum(w*v(*t) for t,w in lam.items())
T=max(sum(v(i,s[i],t[i]) for i in range(3)) for s in itertools.permutations(range(3)) for t in itertools.permutations(range(3)))
print('fractional value', LP, ' T', T, ' gap', LP-T)
assert LP-T==F(3,50)
