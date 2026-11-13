"""Referee check of the sprinkling chain in Theorem C at small q (floating point, evidence only).

Pathwise: (a) sum_{e in M} (d_e - kappa)_+ <= (3q - T) - Phi(delta)
          (b) 3q - T <= Phi(delta') - kappa*k + T''(S)
Means:    (c) E Phi(delta') <= E Phi(delta)/(1-eta);  E T''(S) <= E g(k)/eta
          (d) g(k)/k^(1/3) for small k (is C0 = 4.87 a uniform constant at small sizes?)
"""
import sys, itertools
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix

def amat(k):
    idx = np.arange(k ** 3); a, b, c = np.unravel_index(idx, (k, k, k))
    rows = np.concatenate([a, k + b, 2 * k + c]); cols = np.concatenate([idx, idx, idx])
    return csr_matrix((np.ones(3 * k ** 3), (rows, cols)), shape=(3 * k, k ** 3))

def cost(A, B, C):
    return (A[:, :, None] + B[None, :, :] + C[:, None, :]).reshape(-1)

def perfect(d, k, AM):
    r = milp(d, constraints=LinearConstraint(AM, np.ones(3 * k), np.ones(3 * k)), integrality=np.ones(k ** 3), bounds=Bounds(0, 1), options=dict(mip_rel_gap=0))
    assert r.status == 0
    return r.fun, r.x > 0.5

def phi(d, k, AM, kappa):
    r = milp(d - kappa, constraints=LinearConstraint(AM, -np.inf, np.ones(3 * k)), integrality=np.ones(k ** 3), bounds=Bounds(0, 1), options=dict(mip_rel_gap=0))
    assert r.status == 0
    return r.fun + kappa * k, r.x > 0.5

def main(q, K, eta, reps, seed):
    rng = np.random.default_rng(seed); AM = amat(q); kappa = K * q ** (-2 / 3)
    AMs = {k: amat(k) for k in range(1, q + 1)}
    acc = dict(tail=0, def_=0, phi=0, phip=0, k=0, k13=0, T2=0); viol_a = viol_b = 0
    gk = {}
    for _ in range(reps):
        dp = [1 - (1 - rng.random((q, q))) ** (1 / (1 - eta)) for _ in range(3)]
        dpp = [1 - (1 - rng.random((q, q))) ** (1 / eta) for _ in range(3)]
        d3 = [np.minimum(x, y) for x, y in zip(dp, dpp)]
        d = cost(*d3); defi, M = perfect(d, q, AM); ph, _ = phi(d, q, AM, kappa)
        tail = np.maximum(d[M] - kappa, 0).sum()
        if tail > defi - ph + 1e-9: viol_a += 1
        php, Np = phi(cost(*dp), q, AM, kappa)
        a, b, c = np.unravel_index(np.nonzero(Np)[0], (q, q, q))
        Sa = np.setdiff1d(np.arange(q), a); Sb = np.setdiff1d(np.arange(q), b); Sc = np.setdiff1d(np.arange(q), c)
        k = len(Sa); T2 = 0.0
        if k > 0:
            T2, _ = perfect(cost(dpp[0][np.ix_(Sa, Sb)], dpp[1][np.ix_(Sb, Sc)], dpp[2][np.ix_(Sa, Sc)]), k, AMs[k])
        if defi > php - kappa * k + T2 + 1e-9: viol_b += 1
        for key, v in (('tail', tail), ('def_', defi), ('phi', ph), ('phip', php), ('k', k), ('k13', k ** (1 / 3)), ('T2', T2)): acc[key] += v / reps
    print(f"q={q} K={K} eta={eta} reps={reps}: pathwise violations (a)={viol_a} (b)={viol_b}")
    print(f"   E(3q-T)={acc['def_']:.3f}  E Phi(delta)={acc['phi']:.3f}  E tail={acc['tail']:.3f} <= E(3q-T)-E Phi = {acc['def_']-acc['phi']:.3f}")
    print(f"   E Phi(delta')={acc['phip']:.3f} <= E Phi(delta)/(1-eta) = {acc['phi']/(1-eta):.3f};  E k={acc['k']:.2f}, kappa*E k={kappa*acc['k']:.3f}, E T''(S)={acc['T2']:.3f}")
    sys.stdout.flush()

def gtable(reps, seed):
    rng = np.random.default_rng(seed)
    print("g(k)/k^(1/3), uniform deficits:", end=" ")
    for k in (1, 2, 3, 4, 6, 8, 12, 16):
        AM = amat(k); v = np.mean([perfect(cost(rng.random((k, k)), rng.random((k, k)), rng.random((k, k))), k, AM)[0] for _ in range(reps if k <= 8 else reps // 4)])
        print(f"k={k}: {v / k ** (1/3):.3f}", end="  ")
    print(); sys.stdout.flush()

if __name__ == "__main__":
    gtable(200, 5)
    for q in (8, 12):
        for K, eta in ((2, 0.5), (4, 0.5), (4, 0.25)):
            main(q, K, eta, 150, 1000 + q)
