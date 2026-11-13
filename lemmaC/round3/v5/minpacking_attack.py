"""V5 attack on Lemma 2 accounting and on point 8 (exact enumeration, floating point sums).
Finite unimodular stand-in: a fixed tripartite 3-uniform hypergraph with a uniform random root,
iid discrete costs on teams, m = uniform random minimum-cost perfect matching (ties = extra
randomness). x = E[m | costs]. y^v_e = E[x_e | costs on the labelled r-ball of v], x^(r) = min.
Checks, exactly over all cost configurations:
  loads <= 1;  cost(x^(r)) <= cost(m);  mean slack <= 3 * mean_v E sum_{e ni v} (x_e - y^v_e)_+ ;
  and what happens to the integral rule 1{y > 1/2} (point 8)."""
import itertools, numpy as np

def build(q, extra, seed):
    rng = np.random.default_rng(seed)
    teams = {(i, i, i) for i in range(q)} | {(i, (i + 1) % q, (i + 2) % q) for i in range(q)}
    while len(teams) < 2 * q + extra:
        teams.add(tuple(int(t) for t in rng.integers(0, q, 3)))
    return sorted(teams)

def matchings(q, teams):
    out = []
    by_a = [[k for k, t in enumerate(teams) if t[0] == a] for a in range(q)]
    def rec(a, usedb, usedc, cur):
        if a == q: out.append(list(cur)); return
        for k in by_a[a]:
            _, b, c = teams[k]
            if b not in usedb and c not in usedc:
                rec(a + 1, usedb | {b}, usedc | {c}, cur + [k])
    rec(0, frozenset(), frozenset(), [])
    return out

def run(q=4, extra=3, levels=3, seed=0):
    teams = build(q, extra, seed); E = len(teams)
    PM = matchings(q, teams); P = len(PM)
    Mmat = np.zeros((P, E)); [Mmat[p].__setitem__(PM[p], 1) for p in range(P)]
    cfg = np.array(list(itertools.product(range(levels), repeat=E)), dtype=float) + 1.0  # costs 1..levels
    cost = cfg @ Mmat.T; opt = np.isclose(cost, cost.min(1, keepdims=True))
    x = (opt / opt.sum(1, keepdims=True)) @ Mmat                      # x_e = P(e in m | costs)
    agents = [(s, i) for s in range(3) for i in range(q)]
    inc = {v: [k for k, t in enumerate(teams) if t[v[0]] == v[1]] for v in agents}
    nbr = {v: {(s, teams[k][s]) for k in inc[v] for s in range(3)} for v in agents}
    c_m = (cfg * x).sum(1).mean() / q                                 # per-team cost of m
    print(f"q={q} teams={E} matchings={P} configs={len(cfg)}  P(tie)={np.mean(opt.sum(1) > 1):.3f}  c(m)={c_m:.4f}")
    for r in (1, 2, 3):
        y = {}
        for v in agents:
            ball = {v}
            for _ in range(r - 1): ball = set().union(*[nbr[a] for a in ball])
            bt = sorted({k for a in ball for k in inc[a]})
            key = (cfg[:, bt] - 1) @ (levels ** np.arange(len(bt)))
            key = key.astype(np.int64); cnt = np.bincount(key)
            for k in inc[v]:
                y[v, k] = (np.bincount(key, weights=x[:, k]) / np.maximum(cnt, 1))[key]
        xr = np.stack([np.minimum.reduce([y[(s, teams[k][s]), k] for s in range(3)]) for k in range(E)], 1)
        xi = np.stack([np.minimum.reduce([(y[(s, teams[k][s]), k] > 0.5).astype(float) for s in range(3)]) for k in range(E)], 1)
        load = np.stack([xr[:, inc[v]].sum(1) for v in agents], 1)
        loadi = np.stack([xi[:, inc[v]].sum(1) for v in agents], 1)
        slack = (1 - load).mean(); rhs = 3 * np.mean([sum(np.maximum(x[:, k] - y[v, k], 0) for k in inc[v]).mean() for v in agents])
        print(f"  r={r}: max load={load.max():.6f}  cost(x^r)={(cfg*xr).sum(1).mean()/q:.4f} <= c(m)  "
              f"slack={slack:.4f} <= 3E(x-y)+={rhs:.4f}: {slack <= rhs + 1e-12}   "
              f"integral rule 1{{y>1/2}}: max load={loadi.max():.0f} slack={(1-loadi).mean():.4f}")

run(4, 3, 3, 0)
run(4, 4, 3, 1)
run(5, 2, 2, 2)
