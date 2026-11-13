"""V13: recompute the perturbed exponent of THEOREM_I_PA.md in exact rational arithmetic.
Inputs taken from the text only: K=5, B=50, theta=1.6712, Psi0 <= -3.84751328, and the rows of the table.
Transcendentals are replaced by rational upper bounds proved inside the script:
  ln 2      <= sum_{k<=n} 1/(k 2^k) + 1/((n+1) 2^n)
  exp(-x)   <= 1 / (partial Taylor sum of exp(x)),  x > 0
  ell_max   <= theta*K  (|ell_b| <= theta*K since kappa_b/lambda_b is a mean of exp(-theta x), x in [0,K];
                         |ln J| <= theta*K iff 1 <= S <= exp(2*theta*K)... checked below as 1 <= S <= exp(theta K))
"""
from fractions import Fraction as Fr
def ln2_up(n=120):
    return sum(Fr(1, k * 2**k) for k in range(1, n + 1)) + Fr(1, (n + 1) * 2**n)
def exp_low(x, n=200):          # lower bound of exp(x), x>0 rational
    s, term = Fr(0), Fr(1)
    for j in range(n):
        s += term; term = term * x / (j + 1)
    return s
K, B, th, PSI0 = Fr(5), 50, Fr("1.6712"), Fr("-3.84751328")
eta = K / B
lam_max = eta**3 * (3 * B * B - 3 * B + 1) / 6
S = K * K / (2 * th) + K / th**2 + 1 / th**3
assert 1 <= S <= exp_low(th * K), "need -theta K <= ln J <= 0"
ell_max = th * K
LN2 = ln2_up()
print("lam_max =", lam_max, "=", float(lam_max), " S =", float(S), " ln2_up =", float(LN2))
rows = [("2.30", "14", "1/100000", "1/5000", "1/200000"),
        ("2.302", "17", "1/1000000", "1/50000", "1/4000000"),
        ("2.30224", "21", "1/200000000", "1/2000000", "1/2000000000")]   # third row: mine, not in the text
for beta, c1, d1, d2, tau in rows:
    beta, c1, d1, d2, tau = map(Fr, (beta, c1, d1, d2, tau))
    A = (8 * c1 + 4 + 2 * ell_max) / 3
    d1p = d1 + 2 * tau
    d2p = d2 + 2 * tau * (A + LN2)
    M = lam_max + d1p
    # hypothesis of Lemma E: c1 >= 1 + ln^+ M ; ln^+ M <= max(M-1,0)
    assert c1 >= 1 + max(M - 1, 0)
    tail = (6 * B + 2) / exp_low(1 + c1)
    Err = 3 * d2p + 4 * c1 * d1p + tail
    m = -(th * beta + PSI0)
    tot = th * beta + PSI0 + Err
    # coefficient of ||mu-mu'||_1 in Lemma Q
    coef = 4 * c1 + 2 + ell_max - Fr(3, 2) * A
    print("beta=%s c'=%s: A=%.6f  m=%.10f  Err<=%.10f  exponent<=%.10f  coef=%s  -> %s"
          % (beta, c1, A, m, Err, tot, coef, "NEGATIVE" if tot < 0 else "FAIL"))
    print("    pieces: 3*d2'=%.3e  4c'*d1'=%.3e  tail=%.3e   cost of enlargement per tau = %.2f"
          % (3 * d2p, 4 * c1 * d1p, tail, float(2 * (3 * (A + LN2) + 4 * c1))))
# plain first moment constant: beta^3 e / 27 <= exp(0.22) at beta = 2.31
e_up = 1 / (1 / exp_low(Fr(1)))  # lower bound of e actually; use upper: e <= 2.7182819
e_up = Fr("2.7182819")
assert Fr("2.31")**3 * e_up / 27 <= exp_low(Fr("0.22"))
print("beta^3 e/27 at 2.31 =", float(Fr("2.31")**3 * e_up / 27), "<= exp(0.22) lower bound", float(exp_low(Fr("0.22"))))
# transfer constants
C1 = 36 * K**5 / Fr(1, 10)**6
print("C1 =", float(C1), " K^3 C1/2 =", float(K**3 * C1 / 2))
