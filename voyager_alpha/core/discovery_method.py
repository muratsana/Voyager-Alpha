from dataclasses import dataclass


@dataclass(frozen=True)
class DiscoveryMethod:
    """Documented moving-object workflow defaults used by the independent engine."""

    detector_mode: str = "hybrid"
    detection_sigma: float = 5.0
    detection_fwhm_px: float = 3.0
    max_residuals_per_frame: int = 24
    edge_margin_px: int = 6
    streak_min_area_px: int = 6
    streak_min_elongation: float = 1.8
    min_linked_frames: int = 3
    min_seed_displacement_px: float = 1.5
    match_tolerance_px: float = 2.8
    potential_discovery_rms_px: float = 0.9
    borderline_review_rms_px: float = 1.8
    known_visible_snr: float = 4.0
    known_visible_offset_px: float = 8.0
    gaia_bin_width_mag: float = 0.5
    gaia_samples_per_bin: int = 6
    gaia_required_recoveries: int = 3
    synthetic_max_motion_px_hour: float = 12.0
    synthetic_motion_step_px_hour: float = 1.0
    synthetic_angle_step_deg: float = 30.0
    synthetic_min_stacked_snr: float = 6.0


DOCUMENTED_DISCOVERY_METHOD = DiscoveryMethod()
