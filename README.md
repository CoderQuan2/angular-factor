# The Angular Factor Between the Scalar and Vector Pairing Channels of a Single Contact Interaction

**Author:** T. Abram
**Version:** 1.1 (2026-11-09)
**Cite:** https://doi.org/10.5281/zenodo.24085268 (all versions) — version 1.0: https://doi.org/10.5281/zenodo.24085269
**Mirror:** https://github.com/CoderQuan2/angular-factor

## What this is

A short, self-contained note computing the ratio of the vector (1⁻) to scalar (0⁺) pair susceptibilities of a single attractive contact interaction at a common Fermi surface. The transplant of the scalar gap equation into the vector channel assumes the ratio is 1; the longitudinal angular average suggests 1/3; the computation gives **F = 2/3 exactly** (a ratio of angular weights; the effective-coupling identification is exact at leading logarithm). The vector pair wavefunction is transverse to the relative momentum at every direction, each created pair carries total weight 2, and the trace over the channel's three components gives 4/3 per component against the scalar's 2.

## What changed in version 1.1

Version 1.1 carries the full repair list of two independent adversarial reviews of v1.0. No derived number changed — F = 2/3 stands exactly as published. The repairs are to the note's defenses and its verification artifact:

- A new Discussion paragraph ("One derivation chain, counted once") argues explicitly against the double-counting objection: the Fierz coefficient and the bubble weight are a coefficient and a matrix element of the same operator, and the eigenvalue is invariant under operator rescaling.
- The leading-logarithm qualifier now appears in the abstract and inside the theorem, not only in the scope discussion.
- The Fierz pairing line is displayed in Section 1, with its premise status stated.
- The density-of-states step of the scalar control is displayed, with the pair-summed convention stated.
- One sentence added on condensate-state freedom (unitary vs non-unitary spin-1 states); F is state-independent at the linearised level.
- The abstract's arithmetic now shows the two-pairs step explicitly.
- A stated numerical bound in Section 5 corrected to the verification script's own output (6.4 × 10⁻¹⁶).
- The verification script now checks the vector weight per pair (v1.0 summed the two pairs first), runs to completion (v1.0 terminated early on a formatting error after its main checks), and no longer carries application targets from the parent programme that the note never states.

## Files

- `angular-factor-abram.pdf` — the note (this version).
- `angular-factor-abram.tex` — LaTeX source.
- `pair_wavefunction_verification.py` — the verification script.
- `README.md` — this file.

## Verify it yourself

Run: `python3 pair_wavefunction_verification.py` (requires numpy only).
Expected headline output: scalar weight deviation ≤ 9×10⁻¹⁶; vector per-pair weight deviation ≤ 2×10⁻¹⁵; transversality ≤ 6.4×10⁻¹⁶; susceptibility tensor diagonal 1.3334 / 1.3340 / 1.3326, trace 4.0000, scalar total 2.0000, **F = 0.666667**.

## Scope

Exact at leading logarithm under the stated premises: one contact interaction with equal scalar and vector Fierz coefficients, massless fermions, a common isotropic Fermi surface. Mass corrects the overlaps at order m²c⁴/ε_F². The machinery is standard; the claim is the factor.
