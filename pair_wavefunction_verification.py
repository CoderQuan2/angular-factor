"""
F3 computation: pair-wavefunction overlaps for the 0+ and 1- Fierz channels
of g_F(chi chi-bar)^2, massless chi, pair (k,-k) on the shared Fermi surface.

Conventions (identical for both channels -- this is what makes ratios meaningful):
  Dirac rep; C = i g2 g0; channel operators phi_J = chi^T C Gamma^J chi with
  Gamma = i g5 (0+) and Gamma^i = i g5 g^i (1-), exactly as in the P1a Fierz fit.
  Spinors u(k,h) = (chi_h(khat), h chi_h(khat))/sqrt2, u^dagger u = 1,
  chi_h = helicity eigenspinor, sigma.khat chi_h = h chi_h.
Pair overlap: psi_J(khat) = u^T(k,h1) C Gamma^J u(-k,h2), summed/selected
over the helicity pairs the operator actually creates (nonzero entries).
"""
import numpy as np

rng = np.random.default_rng(20261013)

# Pauli matrices
s1 = np.array([[0,1],[1,0]], complex); s2 = np.array([[0,-1j],[1j,0]], complex)
s3 = np.array([[1,0],[0,-1]], complex); I2 = np.eye(2)
pauli = [s1,s2,s3]

def chi_hel(theta, phi, h):
    """Helicity eigenspinor of sigma.khat, eigenvalue h=+1/-1."""
    if h == +1:
        return np.array([np.cos(theta/2), np.exp(1j*phi)*np.sin(theta/2)])
    else:
        return np.array([-np.exp(-1j*phi)*np.sin(theta/2), np.cos(theta/2)])

def u_spinor(theta, phi, h):
    c = chi_hel(theta, phi, h)
    return np.concatenate([c, h*c])/np.sqrt(2)

# C Gamma matrices (4x4), built from blocks computed by hand:
# C i g5        = diag(s2, s2)
# C i g5 g^i    = [[0, s2 s_i],[-s2 s_i, 0]]
M_scalar = np.block([[s2, 0*I2],[0*I2, s2]])
M_vec = [np.block([[0*I2, s2@pauli[i]],[-s2@pauli[i], 0*I2]]) for i in range(3)]

def overlaps(theta, phi):
    """Return (scalar overlaps dict, vector overlaps dict) keyed by (h1,h2)."""
    th2, ph2 = np.pi-theta, phi+np.pi   # direction of -k
    sc, vec = {}, {}
    for h1 in (+1,-1):
        for h2 in (+1,-1):
            a = u_spinor(theta, phi, h1); b = u_spinor(th2, ph2, h2)
            sc[(h1,h2)] = a @ M_scalar @ b
            vec[(h1,h2)] = np.array([a @ M_vec[i] @ b for i in range(3)])
    return sc, vec

# --- Check 1: structure of nonzero overlaps, |psi|^2 patterns at random khat
print("=== random Fermi-surface directions ===")
maxdev_sc, maxdev_vec_sum, maxdev_trans = 0.0, 0.0, 0.0
for _ in range(200):
    theta = np.arccos(rng.uniform(-1,1)); phi = rng.uniform(0, 2*np.pi)
    khat = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
    sc, vec = overlaps(theta, phi)
    # scalar: nonzero only for h1==h2, |psi|^2 = 1
    for (h1,h2), v in sc.items():
        if h1==h2: maxdev_sc = max(maxdev_sc, abs(abs(v)**2 - 1.0))
        else: assert abs(v) < 1e-12, ("scalar wrong-pair nonzero", v)
    # vector: nonzero only for h1==-h2; weight checked PER PAIR (v1.1 --
    # the deposited v1.0 check summed the two pairs before comparing, which
    # a compensating per-pair error would have passed)
    for (h1,h2), v in vec.items():
        if h1==-h2:
            maxdev_vec_sum = max(maxdev_vec_sum, abs(np.sum(np.abs(v)**2) - 2.0))
            maxdev_trans = max(maxdev_trans, abs(np.dot(v, khat)))  # transversality psi.khat = 0
        else: assert np.all(np.abs(v) < 1e-12), ("vector wrong-pair nonzero", v)
print(f"scalar |psi|^2 = 1 (allowed pairs): max dev over 200 dirs = {maxdev_sc:.2e}")
print(f"vector sum_i|psi_i|^2 = 2, checked per pair: max dev = {maxdev_vec_sum:.2e}")
print(f"vector transversality |psi . khat| = 0: max = {maxdev_trans:.2e}")

# --- Check 2: susceptibility tensor Pi_ij = sum_pairs int dOmega/4pi psi_i psi_j*
# Monte Carlo over sphere (fine grid via random sampling, seed fixed)
N = 400000
cth = rng.uniform(-1,1,N); phi = rng.uniform(0,2*np.pi,N); theta = np.arccos(cth)
Pi = np.zeros((3,3), complex); Pi_sc = 0.0
for t,p in zip(theta,phi):
    sc, vec = overlaps(t,p)
    for (h1,h2), v in vec.items():
        if h1==-h2: Pi += np.outer(v, np.conj(v))
    for (h1,h2), v in sc.items():
        if h1==h2: Pi_sc += abs(v)**2
Pi /= N; Pi_sc /= N
print("\n=== susceptibility (units of radial factor K) ===")
print("Pi_ij (vector channel):"); print(np.round(Pi.real,4))
print(f"trace Pi = {np.trace(Pi).real:.4f}   (expect 4)")
print(f"Pi_ii avg = {np.trace(Pi).real/3:.4f}   (expect 4/3 = {4/3:.4f})")
print(f"Pi_scalar total = {Pi_sc:.4f}   (expect 2)")
F = (np.trace(Pi).real/3)/Pi_sc
print(f"F = Pi_ii / Pi_scalar = {F:.6f}   (expect 2/3 = {2/3:.6f})")

# --- Check 3: what the factor does to the gap (the note's own arithmetic,
# in neutral parameters -- v1.0 of this script carried application targets
# from the parent programme that the note never states; removed in v1.1).
# With A = arsinh(Lambda/Delta) fixed by the scalar channel, a channel whose
# effective coupling is F * lambda_BCS has gap exponent A/F, i.e. the
# suppression exponent rescales as S(F) = A * (1/F - 1).
import math
A = 10.0  # arsinh(Lambda/Delta) for the scalar channel, chosen round
print("\n=== gap-exponent rescaling at fixed scalar anchor A = arsinh(Lambda/Delta) = 10 ===")
for Fv, name in [(1.0, "F = 1   (transplant)"), (2/3, "F = 2/3 (this note)"),
                 (1/2, "F = 1/2"), (1/3, "F = 1/3 (longitudinal average)")]:
    print(f"{name}: S = A(1/F - 1) = {A*(1/Fv - 1):8.4f},  gap ratio exp(-S) = {math.exp(-A*(1/Fv - 1)):.3e}")

# --- Check 4: the inversion identity both channels share (definitional):
# Delta = Lambda / sinh(1/lambda) returns the input gap exactly.
lam = 1/A
Delta_check = 1.0/math.sinh(1/lam)   # in units of Lambda
print(f"\ninversion identity: Lambda/sinh(1/lambda) with lambda = 1/A -> "
      f"Delta/Lambda = {Delta_check:.6e} (input was exp(-A) scale, definitional)")
