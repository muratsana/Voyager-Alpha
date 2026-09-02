# Voyager Alpha - Scientific Architecture

## Product Boundary

Voyager Alpha is a two-module desktop workstation for amateur astronomers:

1. **Asteroid Hunter** finds and reviews moving-object candidates in time-ordered FITS sequences.
2. **Exoplanet Inspection** produces a differential light curve from a target and comparison stars.

Neither module makes an automatic discovery claim. A positive result is a candidate that requires
human review, independent observations and the relevant MPC or exoplanet follow-up process.

## Asteroid Pipeline

1. **Sequence validation**
   - Read `.fit`, `.fits`, `.fts` and `.fits.gz` without memory mapping.
   - Sort by exposure midpoint and reject invalid time or inconsistent image shape.
2. **Calibration**
   - Optional master bias subtraction.
   - Exposure-scaled master dark subtraction.
   - Median-normalized master flat division with invalid-pixel protection.
3. **Astrometry**
   - Use a valid celestial WCS already present in the header.
   - Otherwise solve the reference frame through ASTAP on a temporary copy.
   - Never modify the observer's source FITS file.
4. **Registration**
   - Estimate a subpixel translation by FFT phase correlation.
   - Match stellar centroids and fit a robust affine transform when enough stars exist.
   - Preserve registration RMS, matched-star count and fallback status per frame.
5. **Static sky model**
   - Build a temporal median from a bounded set of aligned sequence samples.
   - Photometrically match each aligned frame before subtraction.
6. **Residual detection**
   - Extract signed residual sources with SEP and robust local noise.
   - Record SNR, flux, area, shape, WCS position and artifact flags.
   - Suppress persistent detector-coordinate residuals and stationary sources.
7. **Time-aware tracklets**
   - Seed velocity hypotheses from frame pairs.
   - Predict positions using actual exposure-midpoint intervals, not frame number.
   - Robustly fit a constant-velocity trajectory with missed-frame tolerance.
8. **Known-object recovery**
   - Query SkyBoT once per field/time and cache the result.
   - Project ephemeris positions into image coordinates and test for local image signal.
   - Match tracklets against known predictions separately from unknown discovery.
9. **Review and confirmation**
   - Inspect original/difference blink at sequence-locked STF.
   - Display only the selected tracklet path to avoid overlapping overlays.
   - Shift-and-stack a selected trajectory with `Synthetic Track` for faint-signal confirmation.
   - Store Accept, Reject or Needs follow-up plus reviewer notes.
10. **Export**
    - CSV and HTML include measurement evidence and review state.
    - ADES is a draft and is limited to accepted candidates with WCS coordinates.

## Exoplanet Pipeline

1. Validate a time-ordered sequence with at least five frames.
2. Apply the same optional bias/dark/flat calibration.
3. Register every image to the selected reference with the asteroid registration engine.
4. Recenter every selected star locally and measure aperture-minus-annulus flux plus uncertainty.
5. Normalize each comparison star independently and build an inverse-scatter weighted ensemble.
6. Apply a robust constant, linear or quadratic baseline while protecting downward transit points.
7. Search contiguous single-transit windows and require baseline coverage, depth significance and SNR.
8. Fit a quadratic limb-darkened physical profile with the MIT-licensed `exoplanet-core` NumPy backend.
9. Report depth uncertainty, duration, midpoint, radius ratio, impact parameter, reduced chi-square and delta BIC.

The current transit classifier is a screening metric. It downweights unstable comparison stars and
models limb darkening, but it does not yet model airmass, meridian flips, color extinction, weather
covariates, a period search across multiple nights or BJD_TDB conversion.

## UI Contract

- Only `Asteroid Hunter` and `Exoplanet Inspection` appear as top-level modules.
- Guided Flow runs dependencies in scientific order; Manual Flow exposes validation, single-frame
  plate solve, known-object recovery and unknown discovery commands.
- Every long-running operation runs outside the UI thread and reports stage plus per-file progress.
- Paused inspection materializes full-resolution data; blink uses high-quality cached previews with
  one stretch solution for the complete sequence.
- Double-clicking the image auto-fits it to the viewport.
- INFO, WARN and ERROR entries remain single-line and color coded.

## Method Provenance

The workflow uses generally accepted astronomical image-processing methods and concepts reviewed in
the Citizen Astronomy asteroid/comet guide. That repository is licensed `CC BY-NC-ND`; no source code,
text, or interface implementation has been copied or adapted. Voyager Alpha's implementation and UI
are independent.
