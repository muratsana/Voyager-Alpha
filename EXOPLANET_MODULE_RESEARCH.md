# Exoplanet Module Research Notes

Reviewed on 2026-09-01. The repositories were inspected as scientific references; no external UI
or application code was copied.

## Repository Assessment

### exoplanet-dev/exoplanet

- MIT-licensed probabilistic modeling toolkit built around PyMC and PyTensor.
- Strong fit for posterior inference, joint transit/RV models, Gaussian processes and complex priors.
- Not embedded in the desktop screening path because the full PyMC stack materially increases EXE
  size, startup time and model runtime. It also requires carefully chosen priors and convergence
  diagnostics that cannot be inferred safely from a short amateur sequence.

### exoplanet-dev/exoplanet-core

- MIT-licensed compiled numerical backend with a lightweight NumPy interface.
- Provides tested Kepler solving and quadratic limb-darkened light-curve calculations.
- Integrated as `exoplanet-core==0.4.0` for deterministic physical transit fitting.

### exoplanet-dev/theano-exoplanet

- Historical Theano compatibility fork, not a transit-analysis application.
- The modern exoplanet stack uses PyTensor/PyMC; this fork is not used by Voyager Alpha.

### exoplanet-dev/transitory

- MIT repository described as quick transit fits, but the reviewed main branch contains an empty
  README and package metadata without a fitting implementation.
- No code or dependency was taken from it.

## Implemented Changes

- Per-frame local stellar centroid refinement.
- Aperture-minus-annulus measurement with empirical uncertainty and flags.
- Stability-weighted comparison-star ensemble instead of a raw sum.
- Constant, linear or quadratic robust baseline detrending.
- Single-event contiguous transit search with pre/post baseline coverage.
- Quadratic limb-darkened physical fit through the exoplanet-core NumPy backend.
- Depth uncertainty, SNR, duration, midpoint, radius ratio, impact parameter, reduced chi-square
  and delta-BIC reporting.
- Raw/detrended light-curve views, error bars, fitted model, comparison weights and quality flags.
- CSV export including uncertainty, validity mask and physical model samples.

## Remaining Scientific Work

- Convert exposure midpoints from JD UTC to BJD TDB using target coordinates and observatory location.
- Measure and detrend airmass, sky background, FWHM, centroid drift and meridian-flip covariates.
- Add multi-night period search and phase folding.
- Add optional advanced PyMC posterior inference with explicit priors and convergence diagnostics.
- Validate against real ETD/AAVSO-style amateur transit datasets and injected-transit recovery tests.
