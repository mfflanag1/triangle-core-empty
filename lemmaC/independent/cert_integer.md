# Independent integer certificate at β = 2.302

From the repository root, run

```sh
python3 -B lemmaC/independent/cert_integer.py
python3 -B lemmaC/independent/cert_integer.py --bits 160 --panels 384
```

The program requires Python 3 and python-flint. It reads the paper’s exact decimal witness from `lemmaC/round3/firstmoment/w_K5_e10.json`. It derives every functional and transfer term from those inputs and the displayed formulas in §§4–5. All transcendental operations and interval arithmetic use FLINT’s Arb real balls. The two commands used for the deposited output ran with python-flint 0.9.0 and FLINT 3.6.0.

At 128 bits and 256 Simpson panels per octave the corrected market exponent lies in \([-0.0001680444,-0.0001679531]\). At 160 bits and 384 panels it lies in \([-0.0001680096,-0.0001679879]\). Both certify the negative exponent required for the lower bound \(2.302\). The fine computation also puts \(\theta(2.30225)+\Psi_0\) in \([0.0000037952,0.0000038169]\) for this specific witness. The sign at that value concerns the unperturbed exponent with the stipulated weights. Other witnesses and the actual asymptotic integer optimum are outside that calculation.

## Integrals and truncation

Set \(w_0=e^{u_\infty}\), \(W=w_0+\sum_bD_bw_b\), and \(\psi(t)=\exp(-tw_0+\sum_b\lambda_b(e^{-tw_b}-1))\). Use Frullani’s identity directly to obtain
\[
\mathbb E\log W=\int_0^\infty [e^{-t}-\psi(t)]\,dt/t.
\]
The other integrands are \(w_b e^{-tw_b}\psi(t)\) for \(\rho_b\), and \(w_0\psi(t)\) for \(\varepsilon_0\). For \(s=2^{-32}\), the absolute error from omitting the head is at most \(s(1+\mathbb EW)\) for the log integral, and \(s\) times the prefactor for the other integrals. The log bound follows from \(|e^{-t}-e^{-tW}|\le t|1-W|\). For \(T=2^{13}\), \(\psi(t)\le e^{-w_0t}\) bounds the discarded log tail in absolute value by \((e^{-T}+e^{-w_0T}/w_0)/T\). The remaining tails are at most \(w_b e^{-(w_0+w_b)T}/(w_0+w_b)\) and \(e^{-w_0T}\), respectively. These head and tail bounds also justify integration of the expected Frullani difference.

## Quadrature error

Each of the four positive integrands \(e^{-t}/t,\ \psi(t)/t,\ w_b e^{-tw_b}\psi(t),\ w_0\psi(t)\) is a Laplace transform \(f\) of a positive measure. For the two integrands involving \(1/t\), multiply the Laplace transform by that of Lebesgue measure, or use \(e^{-tW}/t=\int_W^\infty e^{-ty}\,dy\).

On any slab \([a,2a]\), differentiating four times yields
\[
0\le f^{(4)}(t)\le 81\,a^{-4}f(a/2).
\]
Indeed \(y^4e^{-ay/2}\le (8/e)^4a^{-4}<81a^{-4}\). For an even number \(n\) of uniform subintervals on the slab, the absolute error of composite Simpson integration is at most
\[
\frac{a}{180}(a/n)^4\sup_{[a,2a]}|f^{(4)}|
 \le \frac{81 a}{180n^4} f(a/2).
\]
The program sums these errors for the two components of the log difference and bounds the other terms separately. All quadrature nodes are computed from integer ratios with ball arithmetic. The outer head, tail and discretization errors are added on both sides of the computed balls. The star-mass identity \(\varepsilon_0+\sum_b\lambda_b\rho_b=1\) provides a normalization check.

## Functional and market correction

The antiderivative tail \(F(x)=e^{-\theta x}(1+\theta x+(\theta x)^2/2)/\theta^3\) gives \(\kappa_b=F((b-1)/10)-F(b/10)\) and \(J=F(5)\). Substitute the star enclosures into the expression for \(\Psi_0\) in §4. For §5 take \(c'=17,\ \delta_1=10^{-6},\ \delta_2=1/50000,\ \tau=1/4000000,\ A=(8c'+4+2\theta K)/3\). The program checks each tilt magnitude against \(\theta K\) and checks the required lower bound on \(c'\). It then adds \(3(\delta_2+2\tau(A+\log 2))+4c'(\delta_1+2\tau)+(6B+2)e^{-1-c'}\) to the unperturbed exponent. Success requires the resulting whole ball to be strictly below zero. The output saves the exact upper dyadic endpoint and decimal enclosures rounded outwards.
