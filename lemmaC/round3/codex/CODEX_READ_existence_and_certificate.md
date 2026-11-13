1. **GAP — `cert/cert.py:230–233`; `cert/CERTIFICATE.md` §0: the exact decimal `2.251887` is not certified by the printed output.**  
   `f_hi` rounds the internal bound upward, but `%.6f` subsequently rounds to nearest. A value slightly greater than `2.251887` can therefore print as `2.251887`. This does **not** establish that the claimed inequality is false; it leaves its last digit unsupported. Verify the unformatted `c_hi` against the exact rational `2251887/1000000`, or publish a decimal rounded upward. The margin below `2.30` is unaffected.

2. **GAP — `rdeB/ROUTES.md` §1 and Corollary A: “a solution … hence \(Q^*\in\mathcal Q\).”**  
   \(\mathcal Q\) requires all absolute moments, whereas the cited `rde/RDE.md` §1.2 explicitly establishes only all positive moments and the first negative moment. The extension is straightforward but missing: for integers \(k\ge1\),
   \[
   \mathbb E X_-^k
   \le \frac{\mathbb E Y_+^{k+3}}{(k+1)(k+2)(k+3)}<\infty.
   \]
   Adding this calculation closes the hypothesis needed for Theorem M and Lemma 3.

3. **COSMETIC — `rdeB/ROUTES.md` Lemma 4(i): “\(Y\) … has all exponential moments.”**  
   This is false if it includes negative parameters: \(R(E)=\log E\) near zero, so \(\mathbb E e^{-\lambda R(E)}=\infty\) for \(\lambda\ge1\). What the argument needs—and has—is \(\mathbb E e^{\lambda Y}<\infty\) for every **positive** \(\lambda\). Together with IFR, that forces the upper-tail hazard to diverge. The IFR convolution assertion itself is valid on the real line; see Theorem 3.2 and its following qualification in [Barlow–Marshall–Proschan](https://kylewoodward.com/blog-data/pdfs/references/barlow%2Bmarshall%2Bproschan-the-annals-of-mathematical-statistics-1963A.pdf).

4. **COSMETIC — `rdeB/ROUTES.md` §4: “explicit super-solution.”**  
   \(R\) is explicit, but its multiplier \(a\) is defined through an unspecified threshold and supremum. This proves that a suitable multiplier exists; it does not supply a numerical super-solution. No numerical multiplier is needed for the stated existence limit.

5. **GAP — `cert/cert.py:26,146,222`: the arithmetic justification needs an explicit parameter scope.**  
   For runs A and B, \(h=2^{-10}\) or \(2^{-12}\), so lattice coordinates and mesh coordinates share the \(2^{-20}\) grid; their asserted size bounds make addition exact. The unrestricted CLI accepts other `HEXP` values, while line 146 does not check that `y` shares that grid. Thus the exact-argument justification is established for the advertised runs, not every accepted invocation. Restrict the supported parameters or verify the common-grid condition directly.

I found **no reversed inequality** in the principal mathematical or numerical steps. In particular, I checked:

- **Prékopa and atoms:** \(M_0\), extended by zero where appropriate, remains log-concave for a top atom and for δ₀. Its primitives give the required \(M_1^2\ge M_0M_2\); the likelihood-ratio comparison has the correct direction.
- **Cost comparison:** centering removes the mean of \(\Delta\); the convex supporting-line inequality and nonnegative covariance establish \(g_1\ge g_0\), hence the claimed inequality for \(\Gamma\).
- **Existence limit:** normalization at level 1, domination by \(a|R|\), local convergence of \(g,g'\), and recovery of a fixed point by shifting are consistent. Derivative statements can be read almost everywhere.
- **Certificate:** R1 makes the rounded law **less dispersed** and concave. R3’s tangent terms, R4’s chord weights—including the top atom—and both mean formulas have the claimed directions. The final exponential tail is included.
- **Arithmetic:** for the advertised runs and stated IEEE/BLAS assumptions, the nonnegative-sum error allowance is ample. `rexp` has a valid remainder bound. Floating-point mesh choices become exact law parameters; the other floating-point calculations on the rigorous path are protected by enclosures or conservative error margins.

| Claim | Verdict | Confidence |
|---|---|---|
| **(1) Invariant matching and cost identities** | Sound conditional on a proper fixed point; no substantive defect found. | 95% |
| **(2) Theorem M, cost comparison, existence** | Sound after the minor moment clarification; super-solution is existential rather than numerical. | 94% |
| **(3) Certificate** | Enclosure scheme appears sound for runs A/B. Exact `≤ 2.251887` remains unsupported by six-decimal output alone. | 90% |

**Single weakest step:** converting the internal numerical enclosure into the published exact decimal.

Limitations: this was a source audit, not a reproduction of runs A/B or acceptance checks. Their referenced logs and saved stages were not supplied; the unformatted final bound remains unresolved.