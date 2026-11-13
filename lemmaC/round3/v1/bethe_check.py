"""Referee check of Theorem B (outward Bethe construction), exact arithmetic.

Finite connected piece S of a hypertree (3 agents per team, incidence graph a tree).
x_e on S with sum over e in S at v <= 1; the remainder r_v stands for v's teams outside S
(so interior behaviour of an infinite / infinite-degree hypertree is modelled exactly: the
law of m restricted to S depends on x on S only). For every root agent the law of m|_S is
enumerated and compared (a) across roots, (b) with marginals x_e, (c) with the closed form
prod_v x_{m(v)} * prod_{e in m} x_e^-2 * prod_{f not in m} (1-x_f)^-2.
Edge cases forced with positive probability: x_f = 1, x_g = 0, r_v = 0.
"""
import random, sys
from fractions import Fraction as F
from itertools import product

def build(rng, max_teams):
    teams = []            # triples of agents
    frontier = [0]; nag = 1
    while frontier and len(teams) < max_teams:
        w = frontier.pop(rng.randrange(len(frontier)))
        for _ in range(rng.choice([1, 1, 2, 2, 3])):
            if len(teams) >= max_teams: break
            teams.append((w, nag, nag + 1)); frontier += [nag, nag + 1]; nag += 2
    return teams, nag

def rand_x(rng, teams, nag):
    inc = {v: [] for v in range(nag)}
    for i, t in enumerate(teams):
        for v in t: inc[v].append(i)
    x = {}; r = {}
    def split(total, k, complete):
        # split 'total' among k teams (+ remainder unless complete), allow zeros
        parts = [F(rng.choice([0, 0, 1, 1, 2, 3, 5])) for _ in range(k + (0 if complete else 1))]
        s = sum(parts)
        if s == 0: parts[-1] = F(1); s = F(1)
        parts = [total * p / s for p in parts]
        return parts[:k], (F(0) if complete else parts[k])
    def go(w, f):
        ch = [g for g in inc[w] if g != f]
        total = F(1) if f is None else 1 - x[f]
        complete = bool(ch) and rng.random() < 0.5
        vals, rem = split(total, len(ch), complete)
        r[w] = rem
        for g, val in zip(ch, vals):
            x[g] = val
        for g in ch:
            for u in teams[g]:
                if u != w: go(u, g)
    go(0, None)
    for v in range(nag):
        assert sum(x[g] for g in inc[v]) + r[v] == 1
    return inc, x, r

def conv(d1, d2):
    out = {}
    for (a, p), (b, q) in product(d1.items(), d2.items()):
        if p * q: out[a | b] = out.get(a | b, 0) + p * q
    return out

def law(teams, inc, x, r, root):
    def below(g, w, state):          # agents of g other than w, each reached through g
        d = {frozenset(): F(1)}
        for u in teams[g]:
            if u != w: d = conv(d, sub(u, g, state))
        return d
    def sub(w, f, state):
        ch = [g for g in inc[w] if g != f]
        if state == 'matched':
            d = {frozenset(): F(1)}
            for g in ch: d = conv(d, below(g, w, 'unmatched'))
            return d
        den = F(1) if f is None else 1 - x[f]
        if den == 0:                 # null branch; convention: pick outside S
            opts = [(None, F(1))]
        else:
            opts = [(g, x[g] / den) for g in ch] + [(None, r[w] / den)]
        out = {}
        for pick, p in opts:
            if p == 0: continue
            d = {frozenset(): F(1)}
            for g in ch:
                if g == pick:
                    d = conv(d, conv({frozenset([g]): F(1)}, below(g, w, 'matched')))
                else:
                    d = conv(d, below(g, w, 'unmatched'))
            for k, q in d.items(): out[k] = out.get(k, 0) + p * q
        return out
    return sub(root, None, 'root')

def closed_form(teams, inc, x, r, nag, m):
    covered = {}
    for g in m:
        for v in teams[g]:
            if v in covered: return F(0)
            covered[v] = g
    p = F(1)
    for v in range(nag):
        p *= x[covered[v]] if v in covered else r[v]
    for g in range(len(teams)):
        if g in m:
            if x[g] == 0: return F(0)
            p /= x[g] ** 2
        else:
            if x[g] == 1: return F(0) if p == 0 else None   # 0/0: must already be null
            p /= (1 - x[g]) ** 2
    return p

def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    trials = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    bad_root = bad_marg = bad_cf = edge1 = edge0 = 0; tot_roots = 0
    for t in range(trials):
        teams, nag = build(rng, rng.randint(2, 7))
        inc, x, r = rand_x(rng, teams, nag)
        edge1 += any(v == 1 for v in x.values()); edge0 += any(v == 0 for v in x.values())
        laws = [law(teams, inc, x, r, o) for o in range(nag)]
        tot_roots += nag
        ref = laws[0]
        assert sum(ref.values()) == 1
        for o in range(1, nag):
            if laws[o] != ref:
                bad_root += 1
                if bad_root <= 3: print('ROOT DEPENDENCE', teams, x, r, 'root', o)
        for g in range(len(teams)):
            if sum(p for m, p in ref.items() if g in m) != x[g]: bad_marg += 1
        # closed form over all subsets
        for mask in range(1 << len(teams)):
            m = frozenset(g for g in range(len(teams)) if mask >> g & 1)
            cf = closed_form(teams, inc, x, r, nag, m)
            if cf is None or cf != ref.get(m, F(0)):
                bad_cf += 1
                if bad_cf <= 3: print('CLOSED FORM MISMATCH', teams, x, r, sorted(m), cf, ref.get(m, 0))
    print(f'trials={trials} roots_compared={tot_roots} with_x=1:{edge1} with_x=0:{edge0}')
    print(f'root-dependence failures={bad_root} marginal failures={bad_marg} closed-form failures={bad_cf}')

main()
