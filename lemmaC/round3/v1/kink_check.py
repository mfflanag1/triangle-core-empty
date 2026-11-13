"""Referee check of Theorem K (i), (ii), (iii) in exact arithmetic on finite trees.

Recursion (tree.md section 4): cells g^k_m = xi_k + s^{j_k}_m + s^{l_k}_m, new
s_m = -c_(n+1-m), c_(1)<=...<=c_(n) the n smallest cells over all child teams.
Finite branching K per agent stands in for the Poisson family (the algebra of the claims is
pointwise). Integer / rational data, so ties are exact.
"""
import random, sys
from fractions import Fraction as F

def step(children, n):
    """children: list of (xi, sj, sl). Returns new slopes, kstar, Delta, cells."""
    cells = [[xi + sj[m] + sl[m] for m in range(n)] for xi, sj, sl in children]
    firsts = [c[0] for c in cells]
    ks = min(range(len(cells)), key=lambda k: firsts[k])
    others = [firsts[k] for k in range(len(cells)) if k != ks]
    delta = min(others) - firsts[ks] if others else None
    c = sorted(v for cl in cells for v in cl)[:n]
    return [-c[n - 1 - m] for m in range(n)], ks, delta, cells

def test_i(rng, trials):
    stats = {}
    for _ in range(trials):
        n = rng.choice([2, 3, 4])
        K = rng.choice([2, 3])
        ch = []
        for _ in range(K):
            sj = sorted(rng.randint(-3, 3) for _ in range(n))
            sl = sorted(rng.randint(-3, 3) for _ in range(n))
            ch.append((rng.randint(0, 6), sj, sl))
        s, ks, delta, cells = step(ch, n)
        firsts = [c[0] for c in cells]
        if firsts.count(min(firsts)) > 1: continue          # non-unique argmin: skip
        kap = s[-1] - s[0]
        ksum = (ch[ks][1][-1] - ch[ks][1][0]) + (ch[ks][2][-1] - ch[ks][2][0])
        st = stats.setdefault(n, dict(cases=0, ineq_fail=0, if_fail=0, onlyif_fail=0, n2_fail=0, ex=None))
        st['cases'] += 1
        if kap > ksum: st['ineq_fail'] += 1
        if ksum <= delta and kap != ksum: st['if_fail'] += 1
        if kap == ksum and not ksum <= delta:
            st['onlyif_fail'] += 1
            top_tie = cells[ks][-1] == cells[ks][-2]
            if st['ex'] is None: st['ex'] = (ch, s, delta, ksum, 'top cells of k* tied' if top_tie else 'NO TIE')
            st.setdefault('onlyif_without_tie', 0)
            st['onlyif_without_tie'] += (not top_tie)
        if n == 2 and kap != min(ksum, delta): st['n2_fail'] += 1
    for n in sorted(stats):
        st = stats[n]; ex = st.pop('ex')
        print('K(i) n=%d' % n, st)
        if ex: print('   example (children (xi,sj,sl), new s, Delta, kappa_j*+kappa_l*, note):', ex)

def test_ii(rng, trials, g, n, K):
    ok = bad = notE = 0; sharp_fail = 0
    for _ in range(trials):
        a = sorted(F(rng.randint(-4, 4), 1) for _ in range(n))
        scale = F(1, rng.choice([20, 60, 200, 1000]))
        a = [v * scale for v in a]
        A = max(abs(v) for v in a)
        def build(depth):
            if depth == g:
                X = F(rng.randint(-300, 300), 100)
                return dict(X=X, s=[-X + a[m] for m in range(n)], kids=None)
            kids = [(F(rng.randint(0, 600), 100), build(depth + 1), build(depth + 1)) for _ in range(K)]
            vals = [xi - j['X'] - l['X'] for xi, j, l in kids]
            X = min(vals); ks = vals.index(X)
            sv = sorted(vals); delta0 = sv[1] - sv[0]
            s, ksp, deltap, _ = step([(xi, j['s'], l['s']) for xi, j, l in kids], n)
            return dict(X=X, s=s, kids=kids, ks=ks, delta0=delta0, ksp=ksp, deltap=deltap, depth=depth)
        root = build(0)
        E = True
        stack = [root]
        while stack:
            v = stack.pop()
            if v['kids'] is None: continue
            r = g - v['depth']
            if not v['delta0'] > 2 ** (r + 1) * A: E = False
            xi, j, l = v['kids'][v['ks']]
            stack += [j, l]
        sig = lambda m, t: m if t % 2 == 0 else n - 1 - m
        pred = [-root['X'] + (-2) ** g * a[sig(m, g)] for m in range(n)]
        if E:
            if root['s'] == pred: ok += 1
            else:
                bad += 1
                if bad <= 3: print('K(ii) FAIL on E_g', root['s'], pred)
        else:
            notE += 1
            sharp_fail += (root['s'] != pred)
    print(f'K(ii) g={g} n={n} K={K}: on E_g exact={ok} violations={bad}; off E_g: {notE} cases, formula fails in {sharp_fail}')

def test_iii(rng, trials):
    bad = 0
    for _ in range(trials):
        c = F(rng.randint(-50, 50), 7)
        kids = [(F(rng.randint(0, 600), 100), F(rng.randint(-300, 300), 100), F(rng.randint(-300, 300), 100)) for _ in range(5)]
        X = min(xi - a - b for xi, a, b in kids)
        Xc = min(xi - (a + c) - (b + c) for xi, a, b in kids)
        bad += (Xc != X - 2 * c)
    print('K(iii) shift identity violations:', bad, 'of', trials)

rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
test_i(rng, 200000)
for g, n, K in [(1, 2, 3), (2, 2, 3), (3, 2, 3), (2, 3, 3), (3, 4, 2), (4, 3, 2), (3, 5, 3)]:
    test_ii(rng, 400, g, n, K)
test_iii(rng, 2000)
