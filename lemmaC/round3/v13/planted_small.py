"""V13: pointwise test of TRANSFER Lemma 3 at tiny q (single bin [0,K]).
For fixed N = {(i,i,i)} and fixed d_N, estimate P(H^- = h | d_N) for the most frequent h by simulation,
and compare with rho(h,N) * P'(H'^- = h)   (P' exact).
rho = exp(Delta + K^3(|h|+q)/(2q) + K^6/(36q)) * W^{s(h;N)},  W = 36 q^(2/3)/K  (delta_b = K).
Also prints the sharper internal quantity emp / (W^s P') to compare with exp(Delta + ...).
"""
import sys, math, numpy as np
from collections import Counter
rng = np.random.default_rng(13)
def run(q, K, dN_mode, nsamp, top=14):
    t = K * q ** (-2 / 3); assert t <= 1
    pK = t ** 3 / 6
    Delta = 3 * q ** 3 * (q - 1) * t ** 5 / 40
    W = 36 * q ** (2 / 3) / K
    idx = np.arange(q)
    w = (np.uint64(1) << np.arange(q ** 3, dtype=np.uint64)).reshape(q, q, q)
    diag = np.zeros((q, q, q), bool); diag[idx, idx, idx] = True
    cnt = Counter(); done = 0; bs = 200000
    while done < nsamp:
        AB = rng.random((bs, q, q)); BC = rng.random((bs, q, q)); AC = rng.random((bs, q, q))
        if dN_mode == "zero": v = 0.0
        elif dN_mode == "one": v = 1.0
        elif dN_mode == "cheap": v = 0.75 * t / 3      # every team of N has xi = 0.75 K
        if dN_mode != "random":
            AB[:, idx, idx] = v; BC[:, idx, idx] = v; AC[:, idx, idx] = v
        D = AB[:, :, :, None] + BC[:, None, :, :] + AC[:, :, None, :]
        near = (D <= t) & ~diag
        code = (near * w).sum(axis=(1, 2, 3), dtype=np.uint64)
        u, c = np.unique(code, return_counts=True)
        for a, b in zip(u.tolist(), c.tolist()): cnt[a] += b
        done += bs
    worst = 0.0; rows = []
    for code, c in cnt.most_common(400):
        if c < 200: break
        teams = [(a, b, cc) for a in range(q) for b in range(q) for cc in range(q) if (code >> (a * q * q + b * q + cc)) & 1]
        n = len(teams)
        pairs = set()
        for a, b, cc in teams:
            pairs |= {("AB", a, b), ("BC", b, cc), ("AC", a, cc)}
        free = {p for p in pairs if p[1] != p[2]}
        s = 3 * n - len(free)
        touch = sum(1 for a, b, cc in teams if a == b or b == cc or a == cc)
        Pp = pK ** n * (1 - pK) ** (q ** 3 - q - n)
        emp = c / done; se = math.sqrt(c) / done
        rho = math.exp(Delta + K ** 3 * (n + q) / (2 * q) + K ** 6 / (36 * q)) * W ** s
        r = emp / (rho * Pp); worst = max(worst, r - 3 * se / (rho * Pp))
        rows.append((n, touch, s, emp, emp / Pp, rho, r))
    print("q=%d K=%.2f d_N=%s  samples=%d  Delta=%.4f W=%.1f pK=%.5f  configs tested=%d  max (emp-3se)/bound=%.4f"
          % (q, K, dN_mode, done, Delta, W, pK, len(rows), worst))
    seen = set()
    for n, touch, s, emp, ratio, rho, r in rows:
        key = (n, touch, s)
        if key in seen: continue
        seen.add(key)
        print("    |h|=%d touch=%d s=%d   emp=%.3e  emp/P'=%9.3f   rho=%10.3e   emp/(rho P')=%.4f" % (n, touch, s, emp, ratio, rho, r))
        if len(seen) >= top: break
    sys.stdout.flush()
if __name__ == "__main__":
    for q, K in ((3, 1.0), (3, 1.6), (4, 1.0), (4, 1.6)):
        for mode in ("zero", "cheap", "one", "random"):
            run(q, K, mode, 2000000 if q == 3 else 1200000)
