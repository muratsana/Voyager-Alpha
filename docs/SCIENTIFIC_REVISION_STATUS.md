# Voyager Alpha Scientific Revision Status

Date: 2026-09-02

This document records the implementation status after reviewing
`VoyagerAlpha_bilimsel_paket.zip`. The package is treated as a technical review and
requirements source, not as executable instructions. Every change below was checked
against the existing implementation and covered by focused tests where practical.

## Completed in this revision

- FITS exposure time resolution now honors midpoint keywords first, handles declared
  `TIMESYS`, and records the exposure midpoint instead of silently using `DATE-OBS`.
- FITS observatory coordinates and target coordinates are read when present. Barycentric
  correction, airmass, target altitude, and Sun altitude are calculated with Astropy.
- Exoplanet light-curve times use `BJD_TDB` when WCS and observatory position are
  available; otherwise the result remains explicitly marked `JD_UTC`.
- Exoplanet aperture flux is measured on the calibrated, unwarped detector frame.
  Registration is used only to map selected reference coordinates back to the sensor.
- Transit detrending uses a preliminary transit window and fits the polynomial to
  out-of-transit samples, reducing depth suppression by the trend model.
- Additional FITS camera metadata is recorded: read noise, dark current, saturation
  level, and full-well capacity when the header provides them.
- ASTAP invocation includes field-of-view guidance when known, downsampling, WCS
  output, FITS update, and logging. Machine-readable `PLTSOLVD` and warning values are
  parsed when ASTAP produces an INI file.
- Tracklet linking rejects missing or non-monotonic times, uses physical angular-rate
  limits when plate scale is known, performs uncertainty-weighted motion fits, records
  reduced chi-square, and derives sky position angle when detections have WCS.
- Candidate astrometry uses each original frame's WCS at detector coordinates when
  available. Propagated reference astrometry is explicitly flagged otherwise.
- ADES draft rows use exposure midpoint UTC and include available position uncertainty,
  exposure, log-SNR, and seeing values. The file remains marked as not MPC-ready.
- Product cache paths and generated tracklet identifiers now use Voyager Alpha naming.

## Partial implementations

- ASTAP solution status is parsed, but Gaia residual RMS and matched-star acceptance
  thresholds are not yet measured independently.
- Per-frame WCS is used when present, but every frame is not yet independently solved.
  Reference-WCS propagation remains a flagged fallback.
- The transit model remains a fast deterministic candidate-screening fit. It is not a
  joint physical detrending and posterior-sampling pipeline.
- Position uncertainty combines centroid and registration error where applicable, but
  it does not yet include independently measured plate-solve catalog residuals.

## Required before scientific submission claims

- Master bias/dark/flat construction, calibration compatibility checks, and a persistent
  bad-pixel map.
- Independent per-frame astrometry with catalog residual validation and ADES schema
  validation against the MPC toolchain.
- Topocentric known-object matching over the full track time span with motion and
  position-angle agreement checks.
- Physical transit priors, airmass-aware joint detrending, red-noise beta analysis,
  uncertainty/posterior sampling, and publication-format exports.
- Real FITS regression datasets from multiple cameras, meridian-flip sequences, known
  asteroids, and known transit observations with accepted reference results.

Voyager Alpha remains a candidate detection and review workstation. It must not label a
candidate as a confirmed asteroid discovery or confirmed exoplanet transit without the
independent checks listed above.
