# Theorem I (step I-a): integer lower bound in the Poissonised model

Round 2, line `firstmoment`. Written 2026-11-11. Nothing here has been read by a second
agent or a different model family. Every part is labelled **complete** or **gap**.

Result in one line: in the model with independent team costs, with probability tending
to 1 no perfect assignment has total deficit below 2.302·q^(1/3). The constant is certified
by interval arithmetic (best certified value with the parameters used: 2.30224).

| Part | Content | Label |
|---|---|---|
| 1 | Statement, parameters, profile space | complete |
| 2 | Exact conditional first moment, simplicity conditioning, high option, teams above K | complete |
| 3 | Choice-function bound, number of profiles, typical stars, passage to the limit functional | complete (the limit step is a compactness argument, not quantitative) |
| 4 | Concavity and the tangent bound | complete |
| 5 | Certified evaluation | complete as an interval-arithmetic computation (`certify.py`, mpmath.iv); the script has one author and no independent re-implementation in interval arithmetic |

---

## 1. Statement

### 1.1 Model

Three sectors V_1, V_2, V_3 of q agents. A team is a triple e in V_1 × V_2 × V_3. The q^3
team deficits d_e are independent, each distributed as U + U' + U'' with U, U', U''
independent U[0,1]. The density f of d_e satisfies

  (F1) f(t) = t^2/2 on [0,1],  (F2) f(t) ≤ t^2/2 for all t ≥ 0

(on [1,2], 2f = t^2 − 3(t−1)^2; on [2,3], 2f = (3−t)^2). Rescaled cost: xi_e = q^(2/3)·d_e.
For a perfect assignment N (q pairwise disjoint teams) put X(N) = sum of xi_e over e in N.
The total deficit of N is q^(−2/3)·X(N), so "total deficit below beta·q^(1/3)" is
"X(N) < beta·q".

Why this version: it is the pair-additive deficit law with the dependence between teams
removed, so the transfer step (I-b) compares like with like. The proof uses only (F1), (F2)
and independence; it is identical for every iid deficit law whose density equals t^2/2 on
some interval [0, t_*] and is at most t^2/2 everywhere. By monotone coupling the theorem
then holds for every iid law whose distribution function is pointwise below one of these;
this covers the "first point of a Poisson process of intensity x^2/(2q^2) per team" version,
whose distribution function 1 − exp(−t^3/6) is below min(t^3/6, 1).

Limit object: with Lambda(dx) = (x^2/2) dx, the number of teams at a given agent with
xi_e in a set I ⊂ [0, q^(2/3)] is Binomial(q^2, Lambda(I)/q^2), which tends to
Poisson(Lambda(I)).

### 1.2 Parameters

- Truncation level K > 0; number of bins B; bin width eta = K/B; bins
  I_b = ((b−1)·eta, b·eta], b = 1..B. A team is **low** with bin b if xi_e in I_b, and
  **high** if xi_e > K.
- lambda_b = Lambda(I_b) = eta^3·(3b^2 − 3b + 1)/6. Sum of lambda_b = K^3/6.
- Tilt theta > 0. kappa_b = ∫ over I_b of exp(−theta·x)·x^2/2 dx, and
  J = ∫ from K to infinity of exp(−theta·x)·x^2/2 dx
  = exp(−theta·K)·(K^2/(2 theta) + K/theta^2 + 1/theta^3).
- Weights u = (u_1, ..., u_B, u_inf) in R^(B+1); w_b = exp(u_b), w_inf = exp(u_inf).
- Reference star: D = (D_1, ..., D_B) independent, D_b ~ Poisson(lambda_b).
  W = w_inf + sum_b D_b·w_b.

**Definition 1.1.**
  Phi = E ln W,  rho_b = E[w_b/(W + w_b)],  eps_0 = E[w_inf/W],
  C_b = exp(−3u_b)·kappa_b/lambda_b,  C_inf = exp(−3u_inf)·J,

  Psi_0(K, B, theta, u) = 3·Phi + sum_b lambda_b·ln((1 − rho_b)^3 + C_b·rho_b^3)
                          − 3·eps_0 + C_inf·eps_0^3.

- **Stars.** The binned star of agent v is the vector (D_{v,1}, ..., D_{v,B}), D_{v,b} =
  number of low teams of bin b containing v.
- **Profile space.** A profile is m = (m_1, ..., m_B) in N^B with |m| = sum m_b ≤ q: an
  assignment has profile m if it uses m_b teams of bin b and r = q − |m| high teams. There
  are at most (q+1)^B profiles. Normalised: mu_b = m_b/q, eps = r/q; the limit domain is
  P = {mu in R^(B+1): 0 ≤ mu_b ≤ lambda_b, mu_inf = eps ≥ 0, sum_b mu_b + eps = 1}.

### 1.3 Theorem

**Theorem I (general form).** If theta·beta + Psi_0(K, B, theta, u) < 0 for some K, B,
theta, u, then P(some perfect assignment N has X(N) < beta·q) → 0 as q → infinity.

**Theorem I (certified constant).** For K = 5, B = 50 (eta = 1/10), theta = 1.6712 and the
weights of `w_K5_e10.json` (section 5), Psi_0 ≤ −3.847513. Hence the conclusion holds for
every beta ≤ 2.3022; in particular with probability tending to 1 every perfect assignment
has total deficit at least 2.302·q^(1/3).

The proof gives P ≤ P(stars not typical) + exp(−c·q); the first term tends to 0 at a
polynomial rate (not tracked).

---

## 2. The exact conditional first moment  (complete)

Assume q^(2/3) ≥ K. Put p_b = lambda_b/q^2 = P(xi_e in I_b) (by F1) and
p = sum p_b = K^3/(6q^2).

### 2.1 Tilting and the split into low and high teams

Let G be the labelled low hypergraph: the set of pairs (e, b) with xi_e in I_b. Given G the
costs are independent; a low team of bin b has density x^2/(2·lambda_b) on I_b, a high team
has density g(x)/(1 − p) on (K, infinity) with g(x) ≤ x^2/(2q^2) by (F2). So

  E[exp(−theta·xi_e) | bin b] = kappa_b/lambda_b,
  E[exp(−theta·xi_e) | high] ≤ J/(q^2·(1 − p)).

Let Z = number of perfect assignments with X(N) < beta·q. Since
1{X < beta·q} ≤ exp(theta·(beta·q − X)),

  E[Z | G] ≤ exp(theta·beta·q) · sum over N of prod_{e in N∩G} (kappa_{b(e)}/lambda_{b(e)})
             · (J/(q^2(1−p)))^{|N \ G|}.

Group the assignments by L = N ∩ G (a matching inside G) and U = N \ G. Given L with
|L| = q − r, the set U is a perfect assignment of the remaining r + r + r agents; there are
(r!)^2 of those (the condition that U avoids G is dropped; this is the union bound over
teams above K). With Y_m(G) = number of matchings L ⊂ G with m_b teams of bin b:

  **(2.1)** E[Z | G] ≤ exp(theta·beta·q) · sum over profiles m of
            Y_m(G) · prod_b (kappa_b/lambda_b)^{m_b} · (r!)^2 · (J/(q^2(1−p)))^r.

An agent not covered by L "takes the high option". This is the checker's fix: without it an
agent whose star is empty below K has no choice and the count is wrong.

### 2.2 Configuration model and the simplicity conditioning

A degree array is D = (D^s_{v,b}) (s = 1,2,3; v in V_s; b ≤ B) with the same column sums
n_b = sum_v D^s_{v,b} in each sector. Let n = sum n_b and let 𝒢(D) be the set of labelled
low hypergraphs (distinct teams, one bin each) with degree array D.

(a) Every g in 𝒢(D) has P(G = g) = prod_b p_b^{n_b}·(1 − p)^{q^3 − n}. So, given the binned
stars, G is uniform on 𝒢(D).

(b) Configuration model. For each b and s take half-edges (v, b, i), v in V_s,
1 ≤ i ≤ D^s_{v,b}; there are n_b of them. A configuration is a family of bijections
(sigma_b, tau_b)_b from the sector-1 half-edges of bin b to those of sectors 2 and 3; there
are prod_b (n_b!)^2 configurations. Each configuration is a multiset of labelled triples.
Every g in 𝒢(D) is the image of exactly prod_{s,v,b} D^s_{v,b}! configurations (at each
agent and bin, the teams of g through it are assigned bijectively to its half-edges; teams
of g are distinct, so different assignments give different configurations).

(c) For a configuration c let Ytilde_m(c) be the number of sets of configuration triples,
m_b from bin b, whose agents are all distinct. If c maps to g in 𝒢(D), Ytilde_m(c) = Y_m(g).

**Lemma 2.1 (exact count).** Summing over all configurations with degree array D,

  sum_c Ytilde_m(c) = prod_{s=1}^3 N_s(m; D^s) · prod_b (m_b!·(n_b − m_b)!)^2,

where N_s(m; D^s) is the number of **choice functions** of sector s: ways to select m_b
half-edges of bin b in sector s, for all b, all at distinct agents. Equivalently
N_s(m; D^s) is the coefficient of z^m·z_inf^{q−|m|} in prod_{v in V_s}(z_inf + sum_b D^s_{v,b}·z_b).
In other words E_config[Ytilde_m] = prod_s N_s(m; D^s) · prod_b C(n_b, m_b)^(−2).

Proof. A set of half-edge triples as in (c) is a triple of choice functions together with,
for each b, two bijections between the selected half-edges: prod_s N_s · prod_b (m_b!)^2
sets. Each lies in exactly prod_b ((n_b − m_b)!)^2 configurations. QED.

(Checked by exact enumeration of all 20,736 configurations of a 3+3+3 instance with two
bins, all 10 profiles including partial ones: `check_partial.py`. The perfect-matching case
was enumerated by the checker in `../v2/bf.py`.)

**Lemma 2.2 (simplicity conditioning).** Put
Q(D) = prod_b p_b^{n_b}·(1 − p)^{q^3 − n} · prod_b (n_b!)^2 / prod_{s,v,b} D^s_{v,b}!.
Then for every degree array D and profile m

  E[Y_m(G); D(G) = D] ≤ Q(D) · prod_s N_s(m; D^s) · prod_b C(n_b, m_b)^(−2),

and the sum of Q(D) over all degree arrays is at most exp(K^6/(36·q·(1 − p))) = 1 + o(1).

Proof. By (a)–(c), E[Y_m; D(G) = D] = prod p_b^{n_b}(1−p)^{q^3−n} · sum_{g in 𝒢(D)} Y_m(g)
= prod p_b^{n_b}(1−p)^{q^3−n} · (prod D!)^(−1) · sum over simple c of Ytilde_m(c), and the sum
over simple configurations is at most the sum over all of them; apply Lemma 2.1.
For the total mass, the multinomial theorem gives, for fixed (n_b), sum over arrays of
1/prod D! = prod_b (q^{n_b}/n_b!)^3, so the sum of Q over arrays with column sums (n_b) is
(1 − p)^{q^3 − n}·prod_b (q^3 p_b)^{n_b}/n_b!. Summing over (n_b) in N^B gives at most
(1 − p)^{q^3}·exp(q^3 p/(1 − p)) ≤ exp(q^3 p^2/(1 − p)), and q^3 p^2 = K^6/(36q). QED.

Remark. Q(D) = P(D(G) = D)/P_config(simple | D). Lemma 2.2 says the conditional first moment
given the stars is at most the configuration-model value divided by the probability of
simplicity, and that dividing costs nothing after averaging over D. No lower bound on the
probability of simplicity is needed.

### 2.3 The conditional first-moment bound

Combining (2.1) and Lemma 2.2: for any set 𝒯 of degree arrays,

  **(2.2)** E[Z; D(G) in 𝒯] ≤ (1 + o(1)) · exp(theta·beta·q) · sum over profiles m of
     max over D in 𝒯 of  prod_s N_s(m; D^s) · prod_b C(n_b, m_b)^(−2)
       · prod_b (kappa_b/lambda_b)^{m_b} · (r!)^2 · (J/(q^2(1 − p)))^r.

---

## 3. Exponential bounds  (complete)

### 3.1 Choice functions

For u in R^(B+1) let Phihat_s(u) = (1/q)·sum_{v in V_s} ln(exp(u_inf) + sum_b D^s_{v,b}·exp(u_b)).
The polynomial in Lemma 2.1 has nonnegative coefficients, so for every u

  **(3.1)** N_s(m; D^s) ≤ exp(q·[Phihat_s(u) − sum_b mu_b·u_b − eps·u_inf]).

This is the bound "exp(q[E ln W − ∫ rho ln w dLambda])" of the brief, in binned form, with
the empirical star distribution in place of the Poisson one and the high option included.
It holds for every weight vector and every profile.

### 3.2 Binomials, factorials, number of profiles

C(n, m) ≥ exp(n·h(m/n))/(n + 1) with h the binary entropy (natural logarithm);
r! ≤ e·(r + 1)·(r/e)^r; (1 − p)^(−r) ≤ (1 − p)^(−q) = 1 + o(1). Write nu_b = n_b/q and
g(x, y) = y·h(x/y) for 0 ≤ x ≤ y (g(0,0) = 0; g is continuous). Then each term of (2.2) is
at most poly(q)·exp(q·Xi), where

  Xi(mu; D) = theta·beta + sum_s inf_u [Phihat_s(u) − <mu, u>] − 2·sum_b g(mu_b, nu_b)
              + sum_b mu_b·ln(kappa_b/lambda_b) + eps·ln J + 2·eps·ln eps − 2·eps,

<mu, u> = sum_b mu_b·u_b + eps·u_inf. The polynomial factor is at most
e^2·(q+1)^2·prod_b (n_b + 1)^2, and there are at most (q+1)^B profiles: exp(o(q)) in total
as long as n_b = O(q).

### 3.3 Typical stars

Let nu be the law of the reference star D (product of Poisson(lambda_b)) and nuhat_s the
empirical law of the q binned stars of sector s. For delta > 0 let 𝒯_delta be the set of
degree arrays with

  sum over d in N^B of |nuhat_s(d) − nu(d)|·(1 + |d|) ≤ delta  for s = 1, 2, 3.

**Lemma 3.1.** P(D(G) in 𝒯_delta) → 1 for every delta > 0.

Proof. Within a sector the stars are independent (agents of one sector lie in disjoint sets
of teams) with law nu_q = Multinomial(q^2; p_1, ..., p_B, 1 − p), and the total variation
distance between nu_q and nu is at most q^2·p^2 = K^6/(36q^2) (Le Cam coupling, team by
team). For R ≥ 1: the expected contribution of {|d| ≤ R} is at most
(1 + R)·[sum_{|d| ≤ R} sqrt(nu_q(d)/q) + 2·TV] → 0, and the expected contribution of
{|d| > R} is at most E_{nu_q}[(1+|D|); |D| > R] + E_nu[(1+|D|); |D| > R] ≤ 2C/(1 + R) with
C = (1 + K^3/6)^2 + K^3/6, uniformly in q. Let q → infinity, then R → infinity, and use
Markov's inequality and a union bound over the three sectors. QED.

On 𝒯_delta: (i) |nu_b − lambda_b| ≤ delta; (ii) for A > 0 and u in the box
[−A, 0]^(B+1), exp(−A) ≤ exp(u_inf) + sum_b d_b·exp(u_b) ≤ 1 + |d|, so the logarithm is at
most A + |d| in absolute value and |Phihat_s(u) − Phi(u)| ≤ (A + 1)·delta, where
Phi(u) = E ln(exp(u_inf) + sum_b D_b·exp(u_b)) is the function of Definition 1.1.

### 3.4 The limit functional

For mu with sum_b mu_b + eps = 1 define

  H(mu) = inf over u in R^(B+1) of [Phi(u) − <mu, u>],
  H_A(mu) = the same infimum over u in [−A, 0]^(B+1),
  S_0(mu) = sum_b lambda_b·h(mu_b/lambda_b) − eps·ln eps + eps,
  ell_b = ln(kappa_b/lambda_b), ell_inf = ln J,
  L(mu) = theta·beta + 3·H(mu) − 2·S_0(mu) + <ell, mu>   on P.

Since Phi(u + c·1) = Phi(u) + c and sum mu = 1, the bracket is invariant under shifts of u,
so one may normalise max u = 0 and H_A decreases to H as A → infinity.

**Lemma 3.2.** limsup (1/q)·ln E[Z; D(G) in 𝒯_delta] ≤ L* + o_delta(1), where
L* = sup of L over P and o_delta(1) → 0 as delta → 0. In particular if L* < 0 then
E[Z; D(G) in 𝒯_delta] → 0 exponentially fast for delta small.

Proof. Let 𝒦_delta = {(mu, nu): mu a probability vector on {1..B, inf}, 0 ≤ mu_b ≤ nu_b,
|nu_b − lambda_b| ≤ delta} (compact) and

  Lam_A(mu, nu) = theta·beta + 3·H_A(mu) − 2·sum_b g(mu_b, nu_b) + <ell, mu> + 2 eps ln eps − 2 eps.

H_A is a finite concave function on R^(B+1), hence continuous, so Lam_A is continuous on
𝒦_delta. By 3.1–3.3, for D in 𝒯_delta every profile has Xi(mu; D) ≤ Lam_A(mu, nu) +
3(A + 1)·delta, with (mu, nu) in 𝒦_delta. Hence by (2.2)

  limsup (1/q)·ln E[Z; 𝒯_delta] ≤ sup over 𝒦_delta of Lam_A + 3(A + 1)·delta.

As delta decreases to 0 the compact sets 𝒦_delta decrease to 𝒦_0 = P × {lambda}, so the sup
of the continuous function Lam_A over 𝒦_delta decreases to its sup over P, where
Lam_A = theta·beta + 3H_A − 2S_0 + <ell, ·> (note 2 eps ln eps − 2 eps is the eps-part of
−2·S_0). As A → infinity these continuous functions decrease pointwise to L on the compact
set P, so their suprema decrease to L* (if all suprema were ≥ c, the nested nonempty compact
sets {Lam_A ≥ c} would have a common point, where L ≥ c). Choose A first, then delta. QED.

The step is not quantitative (no explicit A, delta). That is enough for the theorem.

---

## 4. Concavity and the tangent bound  (complete)

Everything in this section concerns the reference star D (independent Poisson(lambda_b))
and is exact for finite B.

Notation. For u in R^(B+1): w = exp(u), W = w_inf + sum_b D_b w_b, and the choice
probabilities pi_b = D_b·w_b/W (b ≤ B), pi_inf = w_inf/W; sum pi = 1.
grad Phi(u) = E pi =: mu(u), and Hess Phi(u) = M(u) := E[diag(pi) − pi·pi^T].
T = {x: sum x = 0}. M·1 = 0, and x^T M x = E[Var_pi(x)] > 0 for x not a multiple of 1
(with positive probability pi charges every coordinate), so M is positive definite on T.

**Lemma 4.1 (Poisson identities).** mu_b(u) = lambda_b·rho_b with rho_b = E[w_b/(W + w_b)]
in (0,1); mu_inf(u) = eps_0 = E[w_inf/W]; and

  E[pi·pi^T] = ∫_0^infinity t·E[exp(−tW)]·v_t·v_t^T dt + diag(lambda_b·E[(w_b/(W + w_b))^2], 0),

where (v_t)_b = lambda_b·w_b·exp(−t·w_b) and (v_t)_inf = w_inf.

Proof. Mecke's formula for Poisson variables: E[D_b·f(D)] = lambda_b·E f(D + e_b),
E[D_b·D_c·f(D)] = lambda_b·lambda_c·E f(D + e_b + e_c) for b ≠ c, and
E[D_b^2·f(D)] = lambda_b^2·E f(D + 2e_b) + lambda_b·E f(D + e_b). Apply them to
f = w_b/W, w_b·w_c/W^2, w_b·w_inf/W^2, and write 1/x^2 = ∫ t·exp(−tx) dt. Entry by entry:
(b,c): lambda_b lambda_c w_b w_c E[(W + w_b + w_c)^(−2)]; (b,b): lambda_b^2 w_b^2
E[(W + 2w_b)^(−2)] + lambda_b E[w_b^2 (W + w_b)^(−2)]; (b,inf): lambda_b w_b w_inf
E[(W + w_b)^(−2)]; (inf,inf): w_inf^2 E[W^(−2)]. QED.

**Lemma 4.2 (covariance domination; the checker's inequality, with the high option).**
M(u) ⪯ Rt(u) := diag(lambda_b·rho_b·(1 − rho_b) (b ≤ B), eps_0) in the positive
semidefinite order.

Proof. The integral in Lemma 4.1 is positive semidefinite, and E[(w_b/(W + w_b))^2] ≥ rho_b^2
by Jensen. So E[pi pi^T] ⪰ diag(lambda_b rho_b^2, 0), and
M = diag(lambda_b rho_b, eps_0) − E[pi pi^T] ⪯ diag(lambda_b rho_b − lambda_b rho_b^2, eps_0). QED.

(The checker's code tests "lambda_max(R^(−1/2) M R^(−1/2)) ≤ 3/2" and finds exactly 1.0000;
Lemma 4.2 is the proof that it is at most 1, and the high-option coordinate is new here.
Floating-point confirmation on 60 random models with 1–3 bins: `check_concave.py`,
maximum 0.999999.)

**Lemma 4.3 (Legendre facts).** Let 𝒟 = grad Phi(R^(B+1)). (i) 𝒟 is relatively open in the
hyperplane {sum = 1}, contained in {0 < mu_b < lambda_b, eps > 0}; on 𝒟, H is C^2,
H(mu(u)) = Phi(u) − <mu(u), u>, its derivative along x in T is −<u, x>, and its second
derivative is −x^T (M|_T)^(−1) x. (ii) If mu_0 in 𝒟 and H(mu) > −infinity, then
mu_t = mu_0 + t(mu − mu_0) lies in 𝒟 for all t in [0,1).

Proof. (i) Phi restricted to T is smooth and strictly convex, so grad Phi is injective on T
with invertible derivative M|_T; apply the inverse function theorem and differentiate
H(mu(u)) = Phi(u) − <mu(u), u>. The inclusion follows from Lemma 4.1. (ii) For u in T,
Phi(u) − <mu_t, u> ≥ (1 − t)·[Phi(u) − <mu_0, u>] + t·H(mu). The bracket is a strictly convex
function on T with a minimiser, so it tends to infinity along every ray; hence so does the
left side, which is continuous and therefore attains its minimum at some u_t, and
grad Phi(u_t) = mu_t. QED.

**Proposition 4.4 (concavity along segments; tangent inequality).** Let u_0 in R^(B+1),
mu_0 = mu(u_0), with rho_b, eps_0 as in Lemma 4.1. For every mu in P,

  **(4.1)** H(mu) − S_0(mu) ≤ Phi(u_0) + sum_b lambda_b·ln(1 − rho_b) − eps_0
            − sum_b mu_b·[u_{0,b} + ln((1 − rho_b)/rho_b)] − eps·[u_{0,inf} − ln eps_0].

Proof. If H(mu) = −infinity there is nothing to prove. Otherwise put x = mu − mu_0 in T and
phi(t) = H(mu_t), a concave function on [0,1], finite at both ends. By Lemma 4.3, on [0,1)
both phi and t ↦ S_0(mu_t) are C^2, with

  phi''(t) = −x^T (M|_T)^(−1) x,   (d/dt)^2 S_0(mu_t) = −x^T Rt^(−1) x,

M, Rt evaluated at the point u_t with mu(u_t) = mu_t (the Hessian of S_0 at a point of 𝒟 is
−diag(1/(lambda_b rho_b (1 − rho_b)), 1/eps) with rho_b = mu_b/lambda_b). Since x ⊥ 1 and
M·1 = 0,

  x^T (M|_T)^(−1) x = sup_y [2<x, y> − y^T M y] ≥ sup_y [2<x, y> − y^T Rt y] = x^T Rt^(−1) x

by Lemma 4.2. So t ↦ phi(t) − S_0(mu_t) is concave on [0,1), hence lies below its tangent
at t = 0:

  phi(t) − S_0(mu_t) ≤ H(mu_0) − S_0(mu_0) − t·<u_0 + grad S_0(mu_0), x>,

grad S_0(mu_0) = (ln((1 − rho_b)/rho_b) for b ≤ B; −ln eps_0). Let t increase to 1: S_0 is
continuous on P, and phi(1) ≤ liminf phi(t) because phi is concave and finite on [0,1].
Finally H(mu_0) = Phi(u_0) − <mu_0, u_0> and
S_0(mu_0) − <grad S_0(mu_0), mu_0> = −sum_b lambda_b·ln(1 − rho_b) + eps_0
(use h(r) − r·h'(r) = −ln(1 − r)). QED.

Consequence (the checker's concavity statement): H − S_0 is concave on 𝒟, so
Sigma = 3H − 2S_0 = 3(H − S_0) + S_0 is concave. Only (4.1) is used below.

(Floating-point confirmation of (4.1): 360 random pairs (u_0, mu) in models with 1–3 bins,
largest value of left side minus right side −1.3e−11: `check_concave.py`.)

**Proposition 4.5 (tangent bound).** For every u_0 in R^(B+1),

  L* = sup over P of L ≤ theta·beta + Psi_0(K, B, theta, u_0).

Proof. L − theta·beta = 3(H − S_0) + S_0 + <ell, mu>. Insert (4.1):

  L(mu) − theta·beta ≤ 3·[Phi(u_0) + sum_b lambda_b ln(1 − rho_b) − eps_0] + S_0(mu) + <a, mu>,
  a_b = ell_b − 3u_{0,b} − 3 ln((1 − rho_b)/rho_b),  a_inf = ell_inf − 3u_{0,inf} + 3 ln eps_0.

Dropping the constraint sum mu = 1 and maximising coordinate by coordinate,
sup_{0 ≤ r ≤ 1} [h(r) + a·r] = ln(1 + e^a) and sup_{eps ≥ 0} [−eps ln eps + eps + a·eps] = e^a:

  S_0(mu) + <a, mu> ≤ sum_b lambda_b·ln(1 + exp(a_b)) + exp(a_inf).

Now (1 − rho_b)^3·(1 + exp(a_b)) = (1 − rho_b)^3 + C_b·rho_b^3 and exp(a_inf) = C_inf·eps_0^3.
QED.

Remarks. (1) A Lagrange multiplier tau for sum mu = 1 would replace a by a − tau and add
tau; this equals the shift u_0 → u_0 + (tau/3)·1, under which Psi_0 is invariant up to that
replacement, so nothing is lost by taking tau = 0. (2) At a stationary point
(3u_b = −theta·c_b + 2 ln(rho_b/(1 − rho_b)) with exp(−theta·c_b) = kappa_b/lambda_b; this
is tree.md's stationarity equation with multiplier 3·beta_tree = theta) the bound is tight,
so minimising Psi_0 over u_0 loses nothing. (3) The tilt uses the exact conditional Laplace
transform kappa_b/lambda_b of the cost within a bin, so binning loses nothing in cost
accounting; the only loss is that the conditioning is on coarser stars. Any B gives a valid
theorem. Floating-point values of the optimised bound −Psi_0/theta at K = 5: 2.30114
(B = 20), 2.30225 (B = 50); the checker's continuum value at K = 5 is 2.30246, and 2.30248
without truncation.

---

## 5. Certified evaluation  (complete as an interval computation)

Script: `certify.py` (mpmath 1.3.0 interval context `iv`, 80-bit, outward rounding).
Inputs are exact rationals: K = 5, eta = 1/10, B = 50, theta = 1.6712, and u from
`w_K5_e10.json` (five-decimal numbers produced by the floating-point search `opt_float.py`;
how they were found is irrelevant, Proposition 4.5 holds for every u_0):

u_1..u_50 = 0.71406, 0.58839, 0.44812, 0.30536, 0.16226, 0.01936, −0.12321, −0.26544,
−0.40739, −0.54916, −0.69086, −0.8326, −0.97452, −1.11672, −1.25934, −1.4025, −1.54629,
−1.69084, −1.83624, −1.98257, −2.12991, −2.27832, −2.42785, −2.57854, −2.73041, −2.88345,
−3.03766, −3.19302, −3.34949, −3.50701, −3.66554, −3.82501, −3.98535, −4.14649, −4.30837,
−4.47092, −4.63405, −4.7977, −4.96183, −5.12636, −5.29125, −5.45647, −5.62197, −5.7877,
−5.95361, −6.11975, −6.28597, −6.45224, −6.6188, −6.78546;  u_inf = −4.61708.

### 5.1 Formulas enclosed

Let psi(t) = E exp(−tW) = exp(−t·w_inf − sum_b lambda_b·(1 − exp(−t·w_b))). Then

- Phi = E ln W = ∫_0^infinity (exp(−t) − psi(t)) dt/t (Frullani)
      = −gamma + ∫_0^1 (1 − psi) dt/t − ∫_1^infinity psi dt/t, gamma = Euler's constant;
- rho_b = ∫_0^infinity w_b·exp(−t·w_b)·psi(t) dt;  eps_0 = ∫_0^infinity w_inf·psi(t) dt;
- lambda_b, kappa_b, J in closed form (the antiderivative of exp(−theta x)·x^2/2 is
  −exp(−theta x)·(x^2/(2 theta) + x/theta^2 + 1/theta^3)).

No truncation of the star size is needed: the Poisson expectation is carried exactly by psi.

### 5.2 Error bounds (every one is rigorous)

psi(t) = E exp(−tW) and exp(−t·w_b)·psi(t) = E exp(−t(W + w_b)) are convex and decreasing in
t. Mesh: t_0 = 2^(−27), t_i = t_{i−1}·(1 + 0.005) in floating point (the nodes are exact
binary numbers, adjacent cells share endpoints), 5,109 cells, up to t_1 ≈ 868.

- On a cell [a, b] with midpoint c: chord ≥ function ≥ tangent at c. For ∫ psi dt/t this
  gives the enclosure [(psi(c) − psi'(c)·c)·ln(b/a) + psi'(c)·(b − a),
  psi(a)·ln(b/a) + ((psi(b) − psi(a))/(b − a))·((b − a) − a·ln(b/a))], with
  psi' = −psi·(w_inf + sum_b lambda_b w_b exp(−t w_b)). For the dt-integrals it gives
  [midpoint rule, trapezoid rule]. The cell [0, t_0] is included for rho_b and eps_0.
- Head of Phi: 0 ≤ ∫_0^{t_0} (1 − psi) dt/t ≤ t_0·E W, E W = w_inf + sum lambda_b w_b ≤ 0.7232.
- Tails: for t ≥ t_1, psi(t) ≤ psi(t_1)·exp(−(t − t_1)·w_inf). So
  ∫_{t_1}^infinity psi dt/t ≤ psi(t_1)/(t_1·w_inf), and the tails of rho_b and eps_0 are at
  most psi(t_1). The mesh stops when both bounds are below 1e−12.
- rho-dependence of Psi_0: r ↦ (1 − r)^3 + C·r^3 and e ↦ −3e + C·e^3 are convex on [0,1]
  and [0, infinity), so over an enclosing interval their maximum is at an endpoint; the
  script takes the larger endpoint value. (This avoids the dependency problem: near the
  optimum these functions are at their minimum in r, so the bound is insensitive to the
  width of the rho enclosures.)

### 5.3 Output (`out_certify_K5_e10.txt`)

  Phi in [−0.5222263134, −0.5222231985]   (width 3.1e−6)
  eps_0 in [0.02081077389, 0.02081090331];  rho_1 in [0.7564957391, 0.7565004437];
  largest relative width of a rho_b enclosure: 6.2e−6
  check: sum_b lambda_b·rho_b + eps_0 in [0.999997927, 1.000004146]  (contains 1)
  **Psi_0 ≤ −3.84751328**
  **−Psi_0/theta ≥ 2.30224586**
  theta·beta + Psi_0 ≤ −0.00375 at beta = 2.30;  ≤ −0.00041 at beta = 2.302;
  ≤ −0.00007 at beta = 2.3022.

Second, smaller certificate (`w_K5_e25.json`, `out_certify_K5_e25.txt`): K = 5, B = 20
(eta = 1/4), theta = 1.6671, 21 weights: Psi_0 ≤ −3.83621828, −Psi_0/theta ≥ 2.30113.

Cross-check (floating point, different quadrature: Frullani integral on a logarithmic grid,
torch): `opt_float.py` gives 2.302248 and 2.301135 for the two parameter sets before the
weights are rounded.

What is rigorous and what is not: the enclosure of Psi_0 is rigorous given the correctness
of mpmath's interval arithmetic and of the 70-line script. The search for the weights, the
cross-check, and the two sanity scripts are floating point and carry no weight in the proof.

---

## 6. Proof of Theorem I  (complete)

Fix K, B, theta, u_0 with theta·beta + Psi_0 < 0. By Proposition 4.5, L* < 0. By Lemma 3.2
there is delta > 0 with E[Z; D(G) in 𝒯_delta] ≤ exp(q·L*/4) for large q. By Lemma 3.1,
P(D(G) not in 𝒯_delta) → 0. So
P(Z ≥ 1) ≤ P(D(G) not in 𝒯_delta) + E[Z; D(G) in 𝒯_delta] → 0.
With the certified parameters, theta·2.3022 + Psi_0 ≤ 3.84743664 − 3.84751328 < 0. QED.

---

## 7. Status

- **Statement proved:** in the Poissonised (independent team cost) model, whp every perfect
  assignment has total deficit at least 2.302·q^(1/3). Certified constant 2.30224; the
  continuum value of the variational problem is 2.30248 (floating point), so the method
  cannot give more than that.
- **Gaps:** none known inside this theorem. Not done: quantitative rate (Lemma 3.2 is a
  compactness argument); independent interval re-implementation of section 5; any read by a
  second agent or model family. The transfer to pair-additive costs (I-b) is a separate line.
- **Differences from tree.md §5 and the checker:** (1) the simplicity conditioning is handled
  by Lemma 2.2 with no lower bound on the probability of simplicity; (2) costs are handled
  by an exponential tilt with exact within-bin Laplace transforms, so there is no binning
  error in cost and coarse bins are valid; (3) concavity is proved in the binned model with
  the high option (Lemmas 4.1–4.2) and used only through the tangent inequality (4.1), which
  avoids boundary and regularity questions for Sigma; (4) the tangent bound holds at any
  weight vector, so the certificate does not need an exact stationary point.

## Files

- `THEOREM_I.md` — this file.
- `certify.py`, `out_certify_K5_e10.txt`, `out_certify_K5_e25.txt` — interval certificate.
- `w_K5_e10.json`, `w_K5_e25.json` — parameters and weights (exact decimals).
- `opt_float.py` — floating-point search for weights (evidence only). `w_K6_e10.json` is
  its output for K = 6, B = 60 (floating-point value 2.302264; not certified, not used).
- Floating-point re-evaluation at the rounded weights (40,001-point logarithmic grid):
  Psi_0 = −3.8475164 (B = 50) and −3.8362214 (B = 20), each 3e−6 below the certified upper
  bound; Monte Carlo over 400,000 Poisson stars gives Phi = −0.52226 ± 0.00103 against the
  enclosure [−0.5222263, −0.5222232].
- `check_partial.py` — exact enumeration check of Lemma 2.1.
- `check_concave.py` — floating-point check of Lemma 4.2 and inequality (4.1).
