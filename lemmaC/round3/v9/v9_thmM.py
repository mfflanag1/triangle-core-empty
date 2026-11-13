"""Adversarial test of Theorem M (and Lemma 3) on random IFR pairs, incl. atoms at the top, bounded supports,
log lower tails, nearly degenerate laws, band perturbations. Also a non-IFR control (power check).
Q'(theta) = sum_i a_i theta^(-b_i) 1{theta < c_i}, 0<=b_i<=1, cutoffs c_i on panel edges (so quadrature stays exact-ish)."""
import numpy as np, sys
from v9_core import *
G = Grid(panel=0.25, k=int(sys.argv[2]) if len(sys.argv) > 2 else 6)
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 2026)
lev = (G.th > 1e-9) & (G.th < 35)
edges_th = np.exp(G.edges)
def F(b, s): return np.log(s) if abs(b - 1) < 1e-12 else (s ** (1 - b) - 1) / (1 - b)
def randterms(nmax=3, atom=None):
    out = []
    for _ in range(rng.integers(1, nmax + 1)):
        a = np.exp(rng.uniform(-3, 2)); b = rng.choice([0.0, 1.0, rng.uniform(0, 1), 2 / 3])
        c = np.inf if rng.random() < 0.5 else edges_th[rng.integers(60, len(edges_th) - 4)]
        out.append((a, b, c))
    if atom is not None: out = [(a, b, min(c, atom)) for a, b, c in out]
    return out
def Qof(terms): return sum(a * F(b, np.minimum(G.th, c)) for a, b, c in terms)
def slopes(q): t, s, P = T(q, G); return t, s, P.last_cond < 1e7
def compare(q0, q1, tag, stats):
    t0, s0, ok0 = slopes(q0); t1, s1, ok1 = slopes(q1)
    ok = lev & ok0 & ok1 & (s0 > 0)
    r = np.where(ok, (s1 - s0) / np.where(s0 > 0, s0, 1), -9.0)[lev]; i = np.argmax(r)
    stats.setdefault('levels tested', [0])[0] += int(ok.sum())
    g0 = mean(t0, G) + 2 * mean(q0, G); g1 = mean(t1, G) + 2 * mean(q1, G)
    st = stats.setdefault(tag, [0, -9.0, None, -9.0]); st[0] += 1
    if r[i] > st[1]: st[1] = r[i]; st[2] = G.th[lev][i]
    st[3] = max(st[3], (g1 - g0) / abs(g0))            # Lemma 3: should be <= 0
stats = {}
for trial in range(int(sys.argv[3]) if len(sys.argv) > 3 else 40):
    # A: Q1 = Q0 + lam*Qb
    ta, tb = randterms(), randterms(); q0 = Qof(ta); compare(q0, q0 + np.exp(rng.uniform(-4, 2)) * Qof(tb), 'A generic', stats)
    # B: Q0 = Q1 truncated at an atom
    c = edges_th[rng.integers(100, len(edges_th) - 6)]; q1 = Qof(ta); compare(Qof([(a, b, min(cc, c)) for a, b, cc in ta]), q1, 'B truncation->atom', stats)
    # C: both with atoms at the top
    c0 = edges_th[rng.integers(110, len(edges_th) - 6)]; tt = randterms(atom=c0); tu = randterms(atom=edges_th[rng.integers(110, len(edges_th) - 6)])
    q0 = Qof(tt); compare(q0, q0 + np.exp(rng.uniform(-3, 1)) * Qof(tu), 'C atoms both', stats)
    # D: nearly degenerate Q0 (eps*Qa) or exactly delta_0, small or large increments
    eps = 10.0 ** rng.uniform(-6, -1); q0 = eps * Qof(ta) * (rng.random() < 0.7)
    compare(q0, q0 + 10.0 ** rng.uniform(-6, 0) * Qof(tb), 'D near-degenerate', stats)
    # E: band increment: Delta' = kap*1{c1<theta<c2} (needs Q0' to jump down by >= kap at c1)
    i1 = rng.integers(90, len(edges_th) - 10); c1 = edges_th[i1]; c2 = edges_th[i1 + rng.integers(1, 6)]; kap = np.exp(rng.uniform(-3, 2))
    base = Qof(ta); compare(base + kap * np.minimum(G.th, c1), base + kap * np.minimum(G.th, c2), 'E band', stats)
    # F: infinitesimal increments (sign of d s/d tau)
    q0 = Qof(ta); compare(q0, q0 + 1e-4 * Qof(tb), 'F infinitesimal', stats)
print('levels tested in total:', stats.pop('levels tested')[0])
for tag, (n, worst, th, g) in stats.items():
    print('%-22s n=%3d  max (s1-s0)/s0 = %+.3e at theta=%.2e   max rel (Gamma1-Gamma0) = %+.3e' % (tag, n, worst, th, g))
# control: NON-IFR inputs (convex pieces in Q: gap in the support / two-point laws); Delta nondecreasing
viol = 0; w = -9
for trial in range(int(sys.argv[3]) if len(sys.argv) > 3 else 40):
    c = edges_th[rng.integers(120, len(edges_th) - 12)]; J = np.exp(rng.uniform(0, 2.5))
    q0 = Qof(randterms()) * 0.2 + J * (G.th > c)                 # jump => not concave
    q1 = q0 + np.exp(rng.uniform(-2, 1)) * np.minimum(G.th, c)   # dispersive increase below the jump
    a0 = slopes(q0); a1 = slopes(q1); ok = lev & a0[2] & a1[2] & (a0[1] > 0)
    r = np.where(ok, (a1[1] - a0[1]) / np.where(a0[1] > 0, a0[1], 1), -9.0); w = max(w, r.max()); viol += r.max() > 1e-6
print('non-IFR control: %d/40 pairs violate, max (s1-s0)/s0 = %+.3e' % (viol, w))
