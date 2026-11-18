# Independent audit of the T1 fractional limit and numerical separation

17 November 2026. Theorem 9.12 passes the focused audit of its two cost comparisons. Independent Arb implementations certify the integer lower bound 2.302 and a fractional upper bound below 2.250535. For the gap $D_q$ this gives $\liminf_{q\to\infty} q^{-1/3}\mathbb E D_q\ge 0.051465$. The exact fractional limit and the optimality of the particular recursive matching remain open.

## Reproduce the calculations

Run from the repository root with Python 3 and python-flint 0.9.0 (FLINT 3.6.0). The finite-cylinder check uses the Python standard library.

~~~sh
python3 -B lemmaC/independent/cert_integer.py --bits 160 --panels 384
python3 -B lemmaC/independent/cert_fractional.py
python3 -B lemmaC/independent/cert_fractional.py --grid-bits 11 --precision 224
python3 -B lemmaC/independent/cert_fractional_checks.py
python3 -B lemmaC/independent/check_cylinders.py
~~~

| Check | Independently certified result | Method and output |
|---|---|---|
| Integer, $\beta=2.302$ | Corrected market exponent in $[-0.0001680096,-0.0001679879]$ | [Quadrature and error proof](cert_integer.md) · [Two precision settings](cert_integer_output.txt) |
| Fractional, starting at $\delta_0$ | Stage-ten upper endpoint $10135505296350745/4503599627370496 < 2.250535$ | [Lattice and cost enclosures](cert_fractional.md) · [Fine run](cert_fractional_run.txt) · [Coarse run](cert_fractional_coarse.txt) |

The integer program reads the paper's published decimal witness. It computes the Poisson-star integrals by Frullani's identity and composite Simpson quadrature with rigorous derivative and tail estimates, then includes the full correction for shared pair costs. A separate 192-bit, 512-panel run put the corrected exponent in $[-0.0001680037,-0.0001679938]$. The 160-bit run also enclosed the *unperturbed* exponent at $2.30225$ in $[0.0000037952,0.0000038169]$. Thus this specific witness and tilt fail the negativity test at that value.

The fractional program constructs ten rounded laws on a fresh uniform mesh. Arb encloses each cell integral. Exact integer polynomial convolution and integer moment sums enclose the hazard and its derivative in the directions required by Lemmas 7.13–7.15. Integration of the resulting lower supporting lines, including the final infinite tail, bounds the cost. An independent replay gave the identical fine endpoint. A coarser grid with higher precision also certifies the published $2.2519$ threshold. An additional direct test integrated the first moments on all 1,085 cells of a three-slope example and checked the cell and chord inequalities. Both implementations use the analytic inequalities of the paper. Each calculates its own enclosures without reading the deposited endpoints, mpmath routines or saved stage arrays.

## Theorem 6.7 gives exact equality of expected costs

Condition on the complete weighted tree. Let $x_e$ be the fractional weight on team $e$. At a chosen root the chain selects $e$ with probability $x_e$. A child whose parent chose their shared team must choose that team. Otherwise the child chooses an alternative $g$ with conditional probability $x_g/(1-x_e)$. Induction gives $\Pr(t(v)=g\mid T,x)=x_g$ at every incidence. A parent team of weight one is selected almost surely, so the exceptional transition has probability zero. Since the agents form a countable set, the chosen teams cover all agents simultaneously.

There is also a direct check that the law is independent of the starting root. Take a finite connected piece of the incidence tree containing all three agents of each included team. Specify a state $t(v)$ at every agent in the piece, allowing a named outside team at the boundary. A compatible pattern selects an included team at all three of its agents or at none. Set $p_e=x_e$ for a selected team and $p_e=1-x_e$ for an unused team. For patterns with positive state and status probabilities, the probability is

$$
\frac{\prod_{v}x_{t(v)}}{\prod_{e\ \mathrm{included}}p_e^2}.
$$

The denominator has two factors for the two child transitions contributed by each team. When the team is selected, each child transition is $x_e/p_e=1$. A zero-probability status or state gives the cylinder probability zero directly. The expression is independent of the choice of root. A finite cylinder and any two candidate roots fit in a finite connected piece of this form. Summing out the unspecified states proves equality on all cylinders. The [exact finite-cylinder check](check_cylinders.py) enumerates 7,776 patterns at each of seven roots, including zero and unit weights. It checks both normalization and the one-agent marginals.

The rank coding makes the resulting rule a measurable probability kernel. Its independence of root and invariance under relabeling pass the adjacent-root mass-transport identity to the matching. Finally, because costs are nonnegative and each conditional edge marginal is $x_e$, Tonelli gives $c(m)=\mathbb E\sum_{e\ni o}\xi_e x_e=c(x)$, including infinite expected cost. The equality concerns expected root cost under the two laws. [Construction and proof](../../manuscript/sections/06_tree.tex)

## Theorem 9.1 gives the upper bound in the finite market

1. Conditional probabilities given the entire cost-marked tree define common weights on each team. The signed adjacent-root identity in Lemma 9.2 forces the expected sum of the absolute endpoint disagreements to be zero. Its integrability follows because a matched team has two other agents.
2. Condition these probabilities on the truncated radius-$r$ ball and use the smallest of each team's three proposals. This makes the local assignment feasible as a packing. Its expected root cost is at most $c(m)$ because cost up to the cutoff is visible on the same ball. The expected root slack tends to zero by martingale convergence for each team rank and domination by the summable full-tree probabilities.
3. At fixed $r$, total-variation convergence of the cost-marked balls transfers both bounded observables to the market. The local rule only needs measurability. Each market team contributes to three of the $3q$ potential roots, so expected root cost on the $q^{2/3}$ cost scale equals expected total market deficit divided by $q^{1/3}$.
4. Two independent families of uniform pair costs support the initial packing and its repair. A single uniform phase rounds the slack in every sector, leaving the same random number $N$ of agents in each. Conditional on these sets and the packing, the unused cost family gives a fresh market of size $N$. Averaging its optimal assignment over the phase fills *every* slack exactly. Its expected cost is at most $C_0(\mathbb E N)^{1/3}$ before the cost-splitting factor. The required expectation estimate for completed assignments is Theorem 2(b) in the authors' uniform-cost version. [Frieze, Pegden and Tkocz](https://www.math.cmu.edu/~ttkocz/mypapers/mathematics/ASS.pdf)
5. First let $q$ grow with $r$ and the splitting parameter $\eta$ fixed. Next let $r$ grow so the normalized patch cost vanishes. Finally let $\eta$ decrease to zero to remove the packing's factor $1/(1-\eta)$. This gives $\limsup_q q^{-1/3}\mathbb E\Phi_q^{\mathrm{LP}}\le c(m)$ for every finite-cost invariant matching. The explicit finite-market estimate also holds when $c(m)=\infty$ because it uses bounded truncated costs. We clarified this last distinction in §9.1. [Transfer proof](../../manuscript/sections/09_lp.tex)

## Closing Theorem 9.12

Write $a_q=q^{-1/3}\mathbb E\Phi_q^{\mathrm{LP}}$. Take a subsequence achieving its liminf and mark each finite team by an optimal LP weight before choosing a uniform root. Compact marks on fixed finite shapes give a subsequential decorated limit. Atom-free Poisson costs make the countably many cutoff restrictions consistent. Counting outgoing triples (root, shared team, neighbor) preserves finite multiplicities in the adjacent-root identity. At fixed cutoff $K$, the finite root degree converges in law to a Poisson variable and its mean equals the limiting $K^3/6$ for all sufficiently large $q$. This provides the uniform integrability needed to pass the identity to the limit.

The expected uncovered root weight after discarding teams of cost above $K$ is at most $a_q/K$. Consequently the limiting weights cover the root fully, and mass transport propagates coverage to all agents. Truncating the cost gives $c(x)\le\liminf a_q$. Theorem 6.7 now supplies a matching $m_*$, and Theorem 9.1 applies after forgetting its auxiliary weights:

$$
\inf_m c(m)\ \le\ c(m_*)=c(x)\ \le\ \liminf_q a_q\ \le\ \limsup_q a_q\ \le\ \inf_m c(m).
$$

This proves both attainment and convergence of the means. The fractional variance bound $\operatorname{Var}(\Phi_q^{\mathrm{LP}})\le4\mathbb E\Phi_q^{\mathrm{LP}}=O(q^{1/3})$ then gives convergence in mean square after normalization. The audit also checked that Theorem 8.4 supplies fixed-cutoff total variation and that its use here requires no uniformity in $r$ or $K$. It takes the paper's earlier counting and fixed-point machinery as dependencies. [Variational limit](../../manuscript/sections/09_lp.tex)

## Manuscript checked

The result, the methods of reproduction and the improved gap are recorded in §11 and Appendix B of the local manuscript. Its 116-page audited PDF has SHA-256 206a7708360d5a9892c9fc4670243ea349f2d1a1b60ed05a0d6c965c33b828e4. A separate build produced the same extracted text (SHA-256 44fffe8bcd3f82604cc91181c95c9375499ec1bc8a080bfe51c3faf8ef5e85fc with `pdftotext -layout`). The proof and certificate pages were inspected after rendering. The preceding 115-page version and its scoped review are identified in the research handoff.
