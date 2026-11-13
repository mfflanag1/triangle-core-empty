"""V12 independent recomputation of Psi_0 (K=5,B=50,theta=1.6712, weights w_K5_e10.json).
Method 1 (no Laplace-transform integrals): the law of S = sum_b D_b w_b is compound Poisson; put the atoms w_b on a
lattice (mean-preserving two-point split), get the law by FFT (exp of the Levy measure), then take E ln W,
E[w_b/(W+w_b)], E[w_inf/W] as plain sums.  Run at two lattice widths.
Method 2: mpmath tanh-sinh quadrature (22 digits) of the Frullani / Laplace forms.
kappa_b, J: closed form AND direct quadrature."""
import json, sys, time, re
import numpy as np, mpmath as mp
LANE = "/Users/maxf/projects/economics-research/game-theory/coopgpt/average_case_protection_2026-11-09/lemmaC/round3"
d = json.load(open(LANE + "/firstmoment/w_K5_e10.json"))
# weights as printed in THEOREM_I.md (section 5) -- compare with the json
md = open(LANE + "/firstmoment/THEOREM_I.md").read()
blk = md[md.index("u_1..u_50 ="):md.index("### 5.1")].replace("−", "-")
nums = [float(x) for x in re.findall(r"-?\d+\.\d+", blk)]
print("weights in THEOREM_I.md == weights in json:", nums[:50] == d["u"] and nums[50] == d["uinf"], len(nums))
B, eta, K, th = d["B"], d["eta"], d["K"], d["theta"]
u = np.array(d["u"]); uinf = d["uinf"]
b = np.arange(1, B + 1)
lam = eta ** 3 * (3 * b * b - 3 * b + 1) / 6
w = np.exp(u); winf = np.exp(uinf)
print("sum lam =", lam.sum(), " K^3/6 =", K ** 3 / 6, " EW =", winf + (lam * w).sum())

def kappa_closed(th):
    A = lambda x: -np.exp(-th * x) * (x * x / (2 * th) + x / th ** 2 + 1 / th ** 3)
    return A(eta * b) - A(eta * (b - 1)), -A(K)

def psi0_from(Phi, rho, eps, kap, J):
    C = np.exp(-3 * u) * kap / lam; Cinf = np.exp(-3 * uinf) * J
    return 3 * Phi + (lam * np.log((1 - rho) ** 3 + C * rho ** 3)).sum() - 3 * eps + Cinf * eps ** 3

kap, J = kappa_closed(th)
# direct quadrature of kappa_b and J (independent of the antiderivative)
mp.mp.dps = 25
kq = [mp.quad(lambda x: mp.exp(-mp.mpf(th) * x) * x * x / 2, [mp.mpf(eta) * (i - 1), mp.mpf(eta) * i]) for i in (1, 2, 17, 50)]
Jq = mp.quad(lambda x: mp.exp(-mp.mpf(th) * x) * x * x / 2, [K, 10, 40, mp.inf])
print("kappa closed vs quad (b=1,2,17,50):", [float(abs(kq[j] - kap[i - 1]) / kq[j]) for j, i in enumerate((1, 2, 17, 50))], " J rel diff:", float(abs(Jq - J) / Jq))

def fft_method(nlog, Smax=32.0):
    N = 1 << nlog; h = Smax / N
    nu = np.zeros(N)
    pos = w / h; k = np.floor(pos).astype(int); fr = pos - k
    np.add.at(nu, k, lam * (1 - fr)); np.add.at(nu, k + 1, lam * fr)
    p = np.fft.irfft(np.exp(np.fft.rfft(nu) - lam.sum()), n=N)
    x = winf + h * np.arange(N)
    Phi = float((p * np.log(x)).sum()); eps = float((p * winf / x).sum())
    rho = np.array([float((p * wb / (x + wb)).sum()) for wb in w])
    return Phi, rho, eps, float(p.sum()), float(p.min()), float(p[-N // 8:].sum())

res = {}
for nlog in (22, 23, 24):
    t = time.time(); Phi, rho, eps, tot, pmin, ptail = fft_method(nlog)
    res[nlog] = psi0_from(Phi, rho, eps, kap, J)
    print("FFT 2^%d: Phi=%.10f eps=%.11f rho1=%.10f sum lam rho+eps=%.12f mass=%.14f minp=%.1e tailmass=%.1e  Psi0=%.9f  -Psi0/theta=%.8f  [%.0fs]"
          % (nlog, Phi, eps, rho[0], (lam * rho).sum() + eps, tot, pmin, ptail, res[nlog], -res[nlog] / th, time.time() - t), flush=True)
print("FFT Richardson (h^2): Psi0 = %.9f" % (res[24] + (res[24] - res[23]) / 3))

# Method 2: tanh-sinh
mp.mp.dps = 22
lam_m = [mp.mpf(3 * i * i - 3 * i + 1) / 6000 for i in range(1, B + 1)]
w_m = [mp.exp(mp.mpf(str(v))) for v in d["u"]]; winf_m = mp.exp(mp.mpf(str(uinf))); th_m = mp.mpf(str(th))
def lpsi(t): return -t * winf_m - mp.fsum(l * (-mp.expm1(-t * x)) for l, x in zip(lam_m, w_m))
pts = [0, mp.mpf(10) ** -6, mp.mpf(10) ** -3, mp.mpf(1) / 10, 1, 4, 16, 64, 256, 1024, 4096, mp.inf]
t = time.time()
Phi_m = mp.quad(lambda t: (mp.expm1(-t) - mp.expm1(lpsi(t))) / t, pts)
eps_m = winf_m * mp.quad(lambda t: mp.exp(lpsi(t)), pts)
rho_m = [x * mp.quad(lambda t: mp.exp(lpsi(t) - t * x), pts) for x in w_m]
def Am(x): return -mp.exp(-th_m * x) * (x * x / (2 * th_m) + x / th_m ** 2 + 1 / th_m ** 3)
e_m = mp.mpf(1) / 10
kap_m = [Am(e_m * i) - Am(e_m * (i - 1)) for i in range(1, B + 1)]; J_m = -Am(mp.mpf(5))
Psi_m = 3 * Phi_m + mp.fsum(l * mp.log((1 - r) ** 3 + mp.exp(-3 * mp.mpf(str(uu))) * kk / l * r ** 3) for l, r, kk, uu in zip(lam_m, rho_m, kap_m, d["u"])) \
        - 3 * eps_m + mp.exp(-3 * mp.mpf(str(uinf))) * J_m * eps_m ** 3
print("tanh-sinh: Phi=%s eps=%s rho1=%s  sum lam rho + eps = %s  [%.0fs]" % (mp.nstr(Phi_m, 14), mp.nstr(eps_m, 13), mp.nstr(rho_m[0], 13), mp.nstr(mp.fsum(l * r for l, r in zip(lam_m, rho_m)) + eps_m, 16), time.time() - t))
print("tanh-sinh: Psi0 = %s   -Psi0/theta = %s   certified upper bound -3.84751328; bound - value = %s" % (mp.nstr(Psi_m, 13), mp.nstr(-Psi_m / th_m, 12), mp.nstr(mp.mpf("-3.84751328") - Psi_m, 5)))
for beta in ("2.30", "2.302", "2.3022", "2.30224", "2.30225"):
    print("  beta=%s  theta*beta+Psi0 = %s" % (beta, mp.nstr(th_m * mp.mpf(beta) + Psi_m, 8)))
