# Proof notes (post-gate extension, 2026-11-11)

Written by a Claude Code lane agent. Not reviewed; not reproduced by a different model.
Setting: three sectors X, Y, Z of q agents each; a full team is e = (x, y, z), one agent per
sector; v(x, y, z) = a(x, y) + b(x, z) + c(y, z) with pair values a, b, c >= 0 (pair-additive
triangle market). T = largest total value of a set of disjoint full teams.
Core (P1, P7): C = {p in R^N : p >= 0, p(N) = T, p(e) >= v(e) for every full team e}.
H_i = max {p_i : p in C}. rho = max {t : p in C, p_i >= t H_i for every i with H_i > 0}.

## Theorem 1 (q = 2: the core is never empty). Status: proof complete, unreviewed.

**LP formulation.** Let LP* = max { sum_e lam_e v(e) : lam >= 0, sum_{e ni i} lam_e <= 1 for every
player i }. By LP duality LP* = min { y(N) : y >= 0, y(e) >= v(e) for all e }. Every core point is
feasible for the dual with y(N) = T, and T <= LP* always; so C is nonempty iff LP* = T (then any
dual optimum is a core point).

**Step 1 (padding).** Because v >= 0 we may restrict to lam with every load equal to 1. Indeed,
for a feasible lam the three sectors carry the same total load sum_e lam_e, so the deficits
d_i = 1 - load_i have equal sums over X, Y and Z. While some deficit is positive, every sector
has a player with positive deficit; add mass min(d_x, d_y, d_z) on the team (x, y, z) formed by
one such player per sector. This zeroes at least one deficit, keeps feasibility and does not
lower the objective (v >= 0). After at most 3q steps all loads are 1.

**Step 2 (coordinates).** Write X = {x1, x2} and so on, lam_{ijk} = lam on (x_i, y_j, z_k).
The objective depends on lam only through its three two-dimensional marginals:
sum_e lam_e v(e) = sum_{x,y} a(x,y) lam_XY(x,y) + sum_{x,z} b(x,z) lam_XZ(x,z) + sum_{y,z} c(y,z) lam_YZ(y,z).
With all loads 1, each marginal is a 2 x 2 doubly stochastic matrix, so it is fixed by one
number: s = lam_XY(x1, y1), t = lam_XZ(x1, z1), u = lam_YZ(y1, z1). Hence the objective is an
affine function K + alpha s + beta t + gamma u of (s, t, u).

**Step 3 (the image of the polytope).** Put w = lam_111. The six load equations (rank 4) give
lam_112 = s - w, lam_121 = t - w, lam_211 = u - w, lam_122 = 1 - s - t + w,
lam_212 = 1 - s - u + w, lam_221 = 1 - t - u + w, lam_222 = s + t + u - 1 - w
(checked symbolically: every load is identically 1 and the marginals are s, t, u).
Nonnegativity asks for a w with
max(0, s+t-1, s+u-1, t+u-1) <= w <= min(s, t, u, s+t+u-1).
Such w exists iff each lower bound is at most each upper bound, i.e. iff
s, t, u in [0, 1], s + t + u >= 1, s + t - u <= 1, s - t + u <= 1, -s + t + u <= 1.
The last four inequalities define the tetrahedron D with vertices (1,1,1), (1,0,0), (0,1,0),
(0,0,1) (each inequality is tight on exactly three of the four vertices), and D lies inside the
unit cube, so the image of the fractional polytope is exactly D.

**Step 4.** The four vertices of D are the images of the four perfect matchings
(sigma = X-to-Y bijection, tau = X-to-Z bijection: s = [sigma = id], t = [tau = id],
u = [sigma = tau]). An affine function on D attains its maximum at a vertex, so LP* is attained
by a perfect matching: LP* = T, and the core is nonempty. QED.

Remark. Pair-additivity is used only in Step 2. Without it the 2 x 2 x 2 axial polytope has the
fractional vertex lam = 1/2 on {111, 122, 212, 221}, whose image is (1/2, 1/2, 1/2), the centre
of D: for pair-additive values it is dominated by the average of the four matchings, but a
general three-sided value can reward it (the Lucas 1995 / Crowley-Graham empty cores at q = 2).
Consistent with gate 1 (2,000 of 2,000 q = 2 instances nonempty).

## Theorem 2 (rho <= 1/3 whenever every H_i > 0). Status: proof complete, unreviewed.

The literal claim "rho <= 1/3 on every triangle market" is FALSE. Exact counterexample (q = 2,
in fact a path market, c = 0): a = [[0,0],[1,1]] (rows x1, x2; columns y1, y2),
b = [[1,0],[1,0]] (rows x1, x2; columns z1, z2), c = 0. Then T = 2. The teams (x2, y, z1) have
value 2 for both y, so p(N) = 2 forces p_x1 = p_y1 = p_y2 = p_z2 = 0; (x2, y, z2) has value 1,
so p_x2 >= 1; (x1, y, z1) has value 1, so p_z1 >= 1. The core is the single point
p = (0, 1 | 0, 0 | 1, 0), H = p, and rho = 1 (only x2 and z1 have H_i > 0). Small random
perturbations of this instance give rho = 1/3 (float, 90 perturbations), so the failure is
degenerate. The correct statement:

**Theorem 2.** In any three-sector market with q agents per sector (any team values v >= 0,
not only pair-additive, and any sector graph) with nonempty core, if H_i > 0 for every player
then rho <= 1/3. More generally, with K sectors of equal size, rho <= 1/K.

**Shift lemma.** If p is in C and alpha_X + alpha_Y + alpha_Z = 0, then
p' = p + alpha_X 1_X + alpha_Y 1_Y + alpha_Z 1_Z is in C provided p' >= 0. Proof: every full
team has exactly one member per sector, so p'(e) = p(e) + alpha_X + alpha_Y + alpha_Z = p(e);
the sectors have equal size q, so p'(N) = p(N) + q (alpha_X + alpha_Y + alpha_Z) = T.

**Proof of Theorem 2.** Let (p, t) be optimal for rho, and let m_S = min_{i in S} p_i be
attained at i_S, for S = X, Y, Z. Put M = m_X + m_Y + m_Z. Shifting by
(m_Y + m_Z, -m_Y, -m_Z) keeps p >= 0, so by the lemma H_{i_X} >= p_{i_X} + m_Y + m_Z = M; in the
same way H_{i_S} >= M for each S. If M > 0, all three H_{i_S} are positive, so
m_S = p_{i_S} >= t H_{i_S} >= t M for each S; summing gives M >= 3 t M, so t <= 1/3. If M = 0,
some m_S = 0, i.e. p_{i_S} = 0 with H_{i_S} > 0 by hypothesis, so t <= 0. Either way
rho <= 1/3. The K-sector version is word for word the same. QED.

Corollaries and links (not separately checked):
- On trees P1 gives rho >= 1/K, so rho = 1/K exactly on trees whenever all H_i > 0; this matches
  the gate 1 path control (rho = 1/3 in all 5,000 instances).
- Gate 1 recorded rho never above 1/3 + 2e-9 in any instance and the P7 link "s = 0 iff
  rho >= 1/3"; with Theorem 2 that link reads "s = 0 iff rho = 1/3" on markets with all H_i > 0.
- The hypothesis H_i > 0 for all i holds whenever some core point has M > 0: for any i in X,
  the same shift gives H_i >= p_i + m_Y + m_Z >= M (likewise in Y and Z). Whether it holds almost surely for continuous
  iid values was not checked as a proof; gate 1 saw smallest H_i = 0.034 over all instances.
