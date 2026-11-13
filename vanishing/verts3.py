"""Exact vertex enumeration of the 3x3x3 axial assignment polytope (pycddlib, exact rationals),
and of its image under the pair-marginal map (the polytope that matters for pair-additive values)."""
import itertools, json, cdd
from fractions import Fraction as F
q = 3
teams = list(itertools.product(range(q), repeat=3))
rows = []
# equalities: load of each agent = 1  ->  1 - sum = 0  as [b, a...] with b + a.x >= 0 / = 0
for s in range(3):
    for i in range(q):
        rows.append([1] + [-1 if t[s] == i else 0 for t in teams])
neq = len(rows)
for j in range(len(teams)):
    rows.append([0] + [1 if k == j else 0 for k in range(len(teams))])
mat = cdd.matrix_from_array(rows, rep_type=cdd.RepType.INEQUALITY, lin_set=set(range(neq)))
poly = cdd.polyhedron_from_matrix(mat)
gen = cdd.copy_generators(poly)
V = [[F(x) for x in g[1:]] for g in gen.array if g[0] == 1]
print('vertices of axial polytope:', len(V))
def marg(lam):
    m = {}
    for t, w in zip(teams, lam):
        for key in [('xy', t[0], t[1]), ('xz', t[0], t[2]), ('yz', t[1], t[2])]:
            m[key] = m.get(key, 0) + w
    keys = sorted(m)
    return tuple(m[k] for k in keys)
imgs = {}
for lam in V:
    imgs.setdefault(marg(lam), lam)
print('distinct marginal images:', len(imgs))
# which images are extreme points of the image polytope: check via cdd on V-rep
pts = list(imgs)
gm = cdd.matrix_from_array([[1] + [x for x in p] for p in pts], rep_type=cdd.RepType.GENERATOR)
red = cdd.matrix_canonicalize(gm)
keep = [i for i in range(len(pts)) if i not in red[0] and i not in red[1]] if False else None
gm2 = cdd.matrix_from_array([[1] + [x for x in p] for p in pts], rep_type=cdd.RepType.GENERATOR)
r = cdd.redundant_rows(gm2)
ext = [i for i in range(len(pts)) if i not in r]
print('extreme image points:', len(ext))
out = []
for i in ext:
    lam = imgs[pts[i]]
    sup = {str(t): str(w) for t, w in zip(teams, lam) if w != 0}
    out.append(dict(integral=all(w in (0, 1) for w in lam), support=sup, marg=[str(x) for x in pts[i]]))
print('fractional extreme image points:', sum(not o['integral'] for o in out))
from collections import Counter
print(Counter(tuple(sorted(Counter(o['support'].values()).items())) for o in out if not o['integral']))
json.dump(out, open('out/verts3.json', 'w'), indent=1)
