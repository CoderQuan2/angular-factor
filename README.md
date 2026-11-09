# The Angular Factor Between the Scalar and Vector Pairing Channels of a Single Contact Interaction

### *A five-page note: the vector channel's pair susceptibility is exactly 2/3 of the scalar channel's — not 1, not 1/3*

<p align="center">
  <img src="https://img.shields.io/badge/Result-F%20%3D%202%2F3%20exactly-blueviolet?style=for-the-badge" alt="F = 2/3 exactly"/>
  <img src="https://img.shields.io/badge/Record-Zenodo%20DOI-orange?style=for-the-badge&logo=zenodo" alt="Zenodo record"/>
  <img src="https://img.shields.io/badge/License-CC%20BY%204.0-yellow?style=for-the-badge" alt="CC BY 4.0"/>
  <img src="https://img.shields.io/badge/Verification-Script%20Included-brightgreen?style=for-the-badge&logo=python" alt="Verification script included"/>
</p>

**Canonical record:** [doi:10.5281/zenodo.24085269](https://doi.org/10.5281/zenodo.24085269) — the Zenodo deposit is the version of record. This repository is the bookkeeping mirror: same files, same version, public history.

---

## 🧭 Overview

A single attractive contact interaction for a massless Dirac fermion pairs in both a scalar (0⁺) and a vector (1⁻) channel, and its Fierz decomposition assigns the two channels the same vertex coefficient. The standard move is therefore to write the vector channel's gap equation by transplanting the scalar one — same equation, coupling relabelled. This note computes the quantity the transplant assumes away: the ratio **F** of the two channels' pair susceptibilities at a common Fermi surface.

Three values suggest themselves before any computation:

| Candidate | Source |
| :--- | :--- |
| F = 1 | The transplant: same vertex coefficient, so same gap equation |
| F = 1/3 | The longitudinal angular average ⟨cos²θ⟩ |
| **F = 2/3** | **The computation** |

The mechanism, in one paragraph: helicity selection pairs the scalar channel with equal helicities and the vector channel with opposite helicities only. The vector pair wavefunction is then transverse to the relative momentum at *every* direction — ψ(k̂) ⊥ k̂ identically, by rotational covariance from a single evaluated direction — and its total weight per pair is 2 everywhere. Rotational invariance forces the vector susceptibility tensor to be a multiple of the identity, and its trace fixes that multiple: each of the three components carries 4/3, against the scalar channel's 2. The 1/3 candidate is the weight of a fixed polarization along the momentum — the longitudinal weight — and transversality makes that weight exactly zero.

$$F \;=\; \frac{4/3}{2} \;=\; \frac{2}{3} \qquad\Longrightarrow\qquad \lambda_{\mathrm{eff}}^{(1^-)} \;=\; \tfrac{2}{3}\,\lambda_{\mathrm{BCS}}$$

The factor is a property of the pair states, not of the interaction vertex — which is why no Fierz analysis of the vertex can contain it.

---

## 📦 Contents

| File | What it is |
| :--- | :--- |
| `angular-factor-abram.pdf` | The note itself — 5 pages, self-contained: premises, derivation, verification table, scope |
| `angular-factor-abram.tex` | The LaTeX source of the note |
| `pair_wavefunction_verification.py` | The numerical checks: spinor-overlap transversality at random directions, and the Monte Carlo evaluation of the susceptibility tensor |

## ✅ Verify it yourself

```bash
python3 pair_wavefunction_verification.py   # overlap checks + Monte Carlo susceptibility tensor
pdflatex angular-factor-abram.tex           # build the note from source (run twice for references)
pdflatex angular-factor-abram.tex
```

Expected results, as printed in the note's verification table: the scalar channel returns F ≡ 1 identically through the same machinery (the control); the susceptibility tensor's diagonal entries are 1.3334, 1.3340, 1.3326 with off-diagonal entries at or below 1.3×10⁻³ and trace exactly 4 (400,000 sampled directions); the Monte Carlo ratio is F = 0.666667; spinor-overlap deviations from the analytic values stay at or below 2×10⁻¹⁵, with wrong-helicity components below 10⁻¹².

## 🎯 Scope — what the note claims and what it does not

- **Exact at leading logarithm.** The angular factor is exact; the full nonlinear gap equation carries an anisotropic remainder of order unity *inside* the logarithm, which moves prefactors, not the factor.
- **Independent of the coupling, the cutoff, and the radial kernel** — all three are channel-common and cancel in the ratio.
- **Mass corrections** enter at order m²c⁴/ε_F²; the massless statement is the leading term, not an approximation smuggled past the reader.
- **No new theory is claimed.** The machinery is standard — pair susceptibilities, Fierz decomposition, helicity spinors. The claim is narrow: a standard shortcut commits a specific, quantified error, and the correct factor is 2/3. The note arose inside a larger modelling project's pairing analysis; it stands or falls entirely on its own contents.

## 📖 Citation

```bibtex
@misc{abram2026angularfactor,
  author       = {Abram, T.},
  title        = {The Angular Factor Between the Scalar and Vector Pairing Channels
                  of a Single Contact Interaction},
  year         = {2026},
  publisher    = {Zenodo},
  version      = {1.0},
  doi          = {10.5281/zenodo.24085269},
  url          = {https://doi.org/10.5281/zenodo.24085269}
}
```

## Status

- **v1.0** — published on Zenodo (version of record); mirrored here as GitHub release v1.0 with identical files.
- Corrections, if any ever land, will be issued as a new version on both records — never as a silent edit.

## License

Creative Commons Attribution 4.0 International (CC BY 4.0) — reuse freely with credit.
