# [111]-aligned DAC low-power reconstruction

## What is reproduced

The optical calculation uses the reconstructed absorption spectra of Ho et al.
under hydrostatic pressure. At fixed incident optical power, the low-power
objective is

\[
A(\lambda,P)=\lambda\,\sigma_{\rm abs}(\lambda,P).
\]

At 120 GPa its interpolant has a numerical maximum at 440.65 nm, reported as
approximately 441 nm because the digitised blue-side nodes are several
nanometres apart. This is the validated optical
baseline, not an orientation-resolved [111] result.

Source: K. O. Ho et al., *Optical Stability and Photophysics of NV Centers in
Diamond up to 120 GPa* (2026), arXiv:2606.02399.
https://arxiv.org/abs/2606.02399

## What [111] literature establishes

Huang et al. calculate ISC rates under general stress and report that uniaxial
[111] stress preserves the NV symmetry and produces the largest contrast; they
also compare the calculation with (111)-cut DAC data at different
hydrostaticities. This constrains the non-optical contrast response, but the
paper does not provide a 120-GPa orientation-resolved absorption spectrum for
reconstructing a new excitation optimum.

Source: B. Huang et al., *Elucidating the Inter-system Crossing of the
Nitrogen-Vacancy Center up to Megabar Pressures* (2026), arXiv:2511.20750.
https://arxiv.org/abs/2511.20750

Wang et al. experimentally retain ODMR operation in a (111)-cut DAC to megabar
pressure and report about 30% contrast at 130 GPa. This is evidence for strong
[111] ODMR contrast, not a measurement of the excitation-action spectrum.

Source: M. Wang et al., *Nature Communications* **15**, 8843 (2024).
https://doi.org/10.1038/s41467-024-52272-y

Davies and Hamer identify the 1.945-eV band as an A1-to-E transition of a
trigonal center through uniaxial-stress spectroscopy. Their result supports the
symmetry treatment but does not supply a megabar [111] absorption kernel.

Source: G. Davies and M. F. Hamer, *Proceedings of the Royal Society A* **348**,
285-298 (1976). https://doi.org/10.1098/rspa.1976.0039

## Claim boundary

Changing a wavelength-independent ODMR contrast rescales sensitivity but does
not move the low-power optimum. Therefore the literature-supported calculation
retains approximately 441 nm at a hydrostatic optical-kernel pressure of
120 GPa. This pressure coordinate is not identified with either the axial or
mean stress of the illustrative tensor. The required [111] optical coupling coefficient
has not been measured in the cited work, so the code provides no adjustable
stress-induced spectral shift and makes no orientation-resolved wavelength
claim beyond this reproduced baseline.

## Conditional finite-power splitting

The existing mediated-response theorem can be applied without inventing an
orientation-dependent optical shift. Define `Ic` as the incident power at which
the maximum-absorption wavelength reaches the response-optimal pumping rate.
For `I <= Ic`, the optimum is the single low-power point. For `I > Ic`, every
solution of

\[
 A(\lambda)/A_{\max}=I_c/I
\]

is an exactly degenerate optimum under the mediation assumption. At 120 GPa,
for example, `I/Ic=1.01` gives 437.49 and 446.36 nm, while `I/Ic=1.10` gives
426.61 and 457.69 nm. The multimodal Ho kernel adds further members near the
ZPL above `I/Ic=1.441`. The absolute value of `Ic` is not predicted for [111];
it requires an experimental power sweep that determines saturation, contrast,
and linewidth scales in the actual DAC.
