"""Was the monotonisation of the validated slopes ever active?  (run A, saved parameters)
For each stage: rebuild the lattice bounds of the saved law, re-evaluate the raw upper bounds
g'_hi(t_i) on the saved mesh, and compare with the saved shat (= running maximum of the raw values).
Also compares with the suffix minimum (the referee's variant).  Both are valid: s(theta) is
nondecreasing, so min_{k>=j} raw_k >= s(theta_j), and max_{k<=j} raw_k >= raw_j >= s(theta_j)."""
import sys
import numpy as np
sys.argv = ['cert.py', '10', '1500', '2500', '8']
import cert
law = (0, np.array([1.0]), np.array([0.0]), np.array([1.0]))
for n in range(1, 9):
    d = np.load('stages/A_stage%d.npz' % n)
    Kmin, p_lo, D_lo, nu_hi = law
    P2lo = cert.lo(cert.conv_nonneg(p_lo, p_lo)); CDlo = cert.lo(cert.conv_nonneg(D_lo, p_lo)); NU2hi = cert.hi(cert.conv_nonneg(nu_hi, nu_hi))
    y = (2 * Kmin + np.arange(len(P2lo))) * cert.hf
    glo, gplo, gphi = cert.evalsums(d['t'], y, P2lo, CDlo, NU2hi)
    pmax = np.maximum.accumulate(gphi); smin = np.minimum.accumulate(gphi[::-1])[::-1]
    print("stage %d: %d levels; theta reproduced: %s; raw g'_hi strictly increasing: %s; running max changed %d entries; "
          "suffix min changed %d entries; saved shat == raw: %s"
          % (n, len(gphi), np.array_equal(glo, d['theta']), bool(np.all(np.diff(gphi) > 0)), int((pmax != gphi).sum()),
             int((smin != gphi).sum()), np.array_equal(d['shat'], gphi)), flush=True)
    law = cert.law_to_lattice(d['theta'], d['shat'])[:4]
