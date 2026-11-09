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
    # vector: nonzero only for h1==-h2
    tot_pairs_sum = 0.0
    for (h1,h2), v in vec.items():
        if h1==-h2:
            tot_pairs_sum += np.sum(np.abs(v)**2)   # per-pair sum_i |psi_i|^2
            maxdev_trans = max(maxdev_trans, abs(np.dot(v, khat)))  # transversality psi.khat = 0
        else: assert np.all(np.abs(v) < 1e-12), ("vector wrong-pair nonzero", v)
    maxdev_vec_sum = max(maxdev_vec_sum, abs(tot_pairs_sum/2 - 2.0))  # per pair: 2
print(f"scalar |psi|^2 = 1 (allowed pairs): max dev over 200 dirs = {maxdev_sc:.2e}")
print(f"vector sum_i|psi_i|^2 = 2 per allowed pair: max dev = {maxdev_vec_sum:.2e}")
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

# --- Check 3: anchors and S(F)
from mpmath import mp, asinh, sinh
mp.dps = 50
Lam = mp.mpf('9.7e26'); Dhad = mp.mpf('4.10e8')  # eV
inv_lam_had = asinh(Lam/Dhad)
print(f"\n1/lambda_had^eff = asinh(Lambda/Delta_had) = {inv_lam_had}  (record: 43.00082)")
for Fv, name in [(mp.mpf(1),'F=1 (R0 null)'), (mp.mpf(2)/3,'F=2/3 (F3 derived)'),
                 (mp.mpf(1)/2,'F=1/2 (R6 asserted)'), (mp.mpf(1)/3,'F=1/3 (R6 asserted)')]:
    S = inv_lam_had*(1/Fv - 1)
    lamG = Fv/inv_lam_had
    print(f"{name}: S = {S:.5f}, lambda_G = {lamG:.6f}")
print("required band: S in [47.4, 49.8], F in [0.46337, 0.47567], lambda_G in [0.010776, 0.011062]")

# --- Check 4: J=0 acceptance -- gap-equation inversion identity (both channels share it)
lam = 1/inv_lam_had
print(f"\nJ=0 inversion: Delta = Lambda/sinh(1/lambda) with lambda=lambda_had^eff -> "
      f"{Lam/sinh(inv_lam_had):.6e} eV (target Delta_had = 4.10e8 eV, definitional)")
lam_eff = (mp.mpf(2)/3)*lam
Dg = Lam/sinh(1/lam_eff)
print(f"J=1 with lambda_eff = (2/3)lambda: Delta_G = {Dg:.6e} eV vs target 3.7e-13 eV "
      f"-> log10 ratio = {float(mp.log10(Dg/mp.mpf('3.7e-13'))):.2f} dex high")
