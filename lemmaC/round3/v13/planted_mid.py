"""V13: q = 12..24.  (a) void bound of Lemma 2 on the full family and on the family F avoiding pairs(N);
(b) class form of Lemma 3: for fixed N and fixed cheap d_N,
    P(X = x, |H^-| = n, s(H^-) = 0 | d_N) <= exp(Delta + K^3(n+q)/(2q) + K^6/(36q)) W^x Bin(x) Bin(n-x)
    (sum of Lemma 3 over h with x touching teams and no sharing inside h; then s(h;N)=x; P' side bounded by product of binomials).
Single bin, W = 36 q^(2/3)/K."""
import sys, math, numpy as np
from collections import Counter
from math import comb, exp, log
rng = np.random.default_rng(1313)
def binom_pmf(n, p, k): return comb(n, k) * p ** k * (1 - p) ** (n - k)
def run(q, K, nsamp, mode):
    t = K * q ** (-2 / 3); pK = t ** 3 / 6
    Delta = 3 * q ** 3 * (q - 1) * t ** 5 / 40
    W = 36 * q ** (2 / 3) / K
    idx = np.arange(q)
    diag = np.zeros((q, q, q), bool); diag[idx, idx, idx] = True
    a, b, c = np.meshgrid(idx, idx, idx, indexing="ij")
    touchmask = ((a == b) | (b == c) | (a == c)) & ~diag
    freemask = ~touchmask & ~diag
    nfree = int(freemask.sum()); ntouch = int(touchmask.sum()); assert ntouch == 3 * q * (q - 1)
    cnt = Counter(); void_full = 0; void_F = 0; done = 0; bs = max(1, 4000000 // q ** 3)
    void_full_u = 0; done_u = 0
    while done < nsamp:
        AB = rng.random((bs, q, q)); BC = rng.random((bs, q, q)); AC = rng.random((bs, q, q))
        # unconditioned market: void bound on the full family, and on F (teams with no pair of N)
        D = AB[:, :, :, None] + BC[:, None, :, :] + AC[:, :, None, :]
        near = D <= t
        void_full_u += int((~near.any(axis=(1, 2, 3))).sum())
        void_F += int((~(near & freemask).any(axis=(1, 2, 3))).sum()); done_u += bs
        v = {"zero": 0.0, "cheap": 0.75 * t / 3}[mode]
        AB[:, idx, idx] = v; BC[:, idx, idx] = v; AC[:, idx, idx] = v
        D = AB[:, :, :, None] + BC[:, None, :, :] + AC[:, :, None, :]
        near = (D <= t) & ~diag
        n = near.sum(axis=(1, 2, 3)); x = (near & touchmask).sum(axis=(1, 2, 3))
        # sharing inside H^-: count pair multiplicities
        mAB = near.sum(axis=3); mBC = near.sum(axis=1); mAC = near.sum(axis=2)
        s = (np.clip(mAB - 1, 0, None).sum(axis=(1, 2)) + np.clip(mBC - 1, 0, None).sum(axis=(1, 2)) + np.clip(mAC - 1, 0, None).sum(axis=(1, 2)))
        for nn, xx, ss in zip(n.tolist(), x.tolist(), s.tolist()):
            cnt[(nn, xx, ss)] += 1
        done += bs
    mu = q ** 3 * pK
    pv = void_full_u / done_u; pvF = void_F / done_u
    print("q=%d K=%.2f mode=%s samples=%d  mu=%.3f Delta=%.4f W=%.1f" % (q, K, mode, done, mu, Delta, W))
    print("   void, full family: emp=%.5f (+-%.5f)  IT=(1-p)^{q^3}=%.5f  Lemma2 bound exp(-mu+Delta)=%.5f"
          % (pv, math.sqrt(pv * (1 - pv) / done_u), (1 - pK) ** (q ** 3), exp(-mu + Delta)))
    print("   void, family F (no pair of N), |F|=%d: emp=%.5f  IT=%.5f  bound exp(-|F|p+Delta)=%.5f"
          % (nfree, pvF, (1 - pK) ** nfree, exp(-nfree * pK + Delta)))
    worst = 0; shown = 0
    for (nn, xx, ss), c in sorted(cnt.items(), key=lambda kv: -kv[1]):
        if ss != 0 or c < 100: continue
        emp = c / done; se = math.sqrt(c) / done
        Pp = binom_pmf(ntouch, pK, xx) * binom_pmf(nfree, pK, nn - xx)
        bound = exp(Delta + K ** 3 * (nn + q) / (2 * q) + K ** 6 / (36 * q)) * W ** xx * Pp
        worst = max(worst, (emp - 3 * se) / bound)
        if shown < 6:
            print("   n=%d x=%d s=0: emp=%.3e  P'(n,x)=%.3e  emp/P'=%.3e  (3q^(2/3)/K)^x=%.3e  W^x=%.3e  emp/bound=%.3e"
                  % (nn, xx, emp, Pp, emp / Pp, (3 * q ** (2 / 3) / K) ** xx, W ** xx, emp / bound)); shown += 1
    print("   max over classes (>=100 events) of (emp-3se)/bound = %.4e" % worst); sys.stdout.flush()
if __name__ == "__main__":
    for q, K, ns in ((12, 0.8, 150000), (16, 1.0, 100000), (24, 1.0, 50000), (24, 1.3, 50000)):
        for mode in ("cheap", "zero"):
            run(q, K, ns, mode)
