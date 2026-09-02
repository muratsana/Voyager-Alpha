from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np


@dataclass
class CameraMetadata:
    instrument: Optional[str] = None
    detector: Optional[str] = None
    pixel_size_x_um: Optional[float] = None
    pixel_size_y_um: Optional[float] = None
    binning_x: Optional[int] = None
    binning_y: Optional[int] = None
    gain: Optional[float] = None
    gain_keyword: Optional[str] = None
    offset: Optional[float] = None
    sensor_temperature_c: Optional[float] = None
    set_temperature_c: Optional[float] = None
    filter_name: Optional[str] = None
    focal_length_mm: Optional[float] = None
    aperture_mm: Optional[float] = None
    image_scale_arcsec_px: Optional[float] = None
    image_scale_source: Optional[str] = None
    readout_mode: Optional[str] = None
    bayer_pattern: Optional[str] = None


@dataclass
class FrameRecord:
    index: int
    file_path: str
    date_obs: Optional[str]
    midpoint_jd: Optional[float]
    exposure_seconds: Optional[float]
    shape: tuple[int, int]
    has_wcs: bool
    quality_flags: list[str] = field(default_factory=list)
    registration_dx: float = 0.0
    registration_dy: float = 0.0
    registration_peak: float = 1.0
    registration_rms_px: Optional[float] = None
    matched_stars: int = 0
    background_median: Optional[float] = None
    background_rms: Optional[float] = None
    fwhm_px: Optional[float] = None
    calibration_state: str = "uncalibrated"
    camera: CameraMetadata = field(default_factory=CameraMetadata)


@dataclass
class Detection:
    frame_index: int
    x: float
    y: float
    ra: Optional[float]
    dec: Optional[float]
    flux: float
    snr: float
    fwhm_px: float
    eccentricity: float
    area_px: int = 0
    peak: float = 0.0
    local_rms: float = 0.0
    flags: list[str] = field(default_factory=list)
    sensor_x: Optional[float] = None
    sensor_y: Optional[float] = None

    def to_overlay(self, is_known: bool = False, match: Optional[dict[str, Any]] = None) -> dict[str, Any]:
        return {
            "x": self.x,
            "y": self.y,
            "ra": self.ra,
            "dec": self.dec,
            "flux": self.flux,
            "snr": self.snr,
            "fwhm_px": self.fwhm_px,
            "eccentricity": self.eccentricity,
            "area_px": self.area_px,
            "peak": self.peak,
            "is_known": is_known,
            "match": match,
            "flags": self.flags,
        }


@dataclass
class Tracklet:
    tracklet_id: str
    detections: list[Detection]
    frames_detected: int
    motion_px_per_frame: float
    motion_arcsec_per_min: Optional[float]
    position_angle_deg: float
    fit_rms_px: float
    median_snr: float
    confidence: float
    velocity_x_px_min: float = 0.0
    velocity_y_px_min: float = 0.0
    arc_minutes: float = 0.0
    classification: str = "unknown_candidate"
    artifact_flags: list[str] = field(default_factory=list)
    known_match: Optional[dict[str, Any]] = None
    review_status: str = "unreviewed"
    reviewer_notes: str = ""

    @property
    def first_frame(self) -> int:
        return min(d.frame_index for d in self.detections)

    @property
    def last_frame(self) -> int:
        return max(d.frame_index for d in self.detections)

    def to_summary(self) -> dict[str, Any]:
        return {
            "id": self.tracklet_id,
            "frames": self.frames_detected,
            "first_frame": self.first_frame,
            "last_frame": self.last_frame,
            "motion_px_per_frame": self.motion_px_per_frame,
            "motion_arcsec_per_min": self.motion_arcsec_per_min,
            "position_angle_deg": self.position_angle_deg,
            "fit_rms_px": self.fit_rms_px,
            "median_snr": self.median_snr,
            "confidence": self.confidence,
            "velocity_x_px_min": self.velocity_x_px_min,
            "velocity_y_px_min": self.velocity_y_px_min,
            "arc_minutes": self.arc_minutes,
            "classification": self.classification,
            "artifact_flags": list(self.artifact_flags),
            "known_match": self.known_match,
            "review_status": self.review_status,
        }


@dataclass
class RegistrationSolution:
    """Affine transform mapping moving-frame x/y coordinates into the reference."""

    matrix_xy: np.ndarray
    offset_xy: np.ndarray
    rms_px: float
    matched_stars: int
    method: str
    phase_peak: float = 0.0

    @classmethod
    def identity(cls) -> "RegistrationSolution":
        return cls(
            matrix_xy=np.eye(2, dtype=np.float64),
            offset_xy=np.zeros(2, dtype=np.float64),
            rms_px=0.0,
            matched_stars=0,
            method="reference",
            phase_peak=1.0,
        )


@dataclass
class KnownObjectPrediction:
    name: str
    ra: float
    dec: float
    magnitude: Optional[float]
    object_type: str = "asteroid"
    motion_ra_arcsec_hour: Optional[float] = None
    motion_dec_arcsec_hour: Optional[float] = None
    x: Optional[float] = None
    y: Optional[float] = None
    local_snr: Optional[float] = None
    local_offset_px: Optional[float] = None
    expected_trail_px: Optional[float] = None
    search_radius_px: int = 5
    near_edge: bool = False
    visible: bool = False
    confidence: float = 0.0
    status: str = "predicted_in_field"
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass
class SyntheticTrackResult:
    image: np.ndarray
    used_frames: int
    skipped_frames: int
    velocity_x_px_min: float
    velocity_y_px_min: float
    snr: float
    peak_offset_px: float
    center_x: float
    center_y: float


@dataclass
class SequenceResult:
    frames: list[FrameRecord]
    reference_index: int
    reference_header: Any
    known_objects: list[KnownObjectPrediction] = field(default_factory=list)
    tracklets: list[Tracklet] = field(default_factory=list)
    limiting_magnitude: Optional[float] = None
    registration_solutions: dict[int, RegistrationSolution] = field(default_factory=dict)
    static_model: Optional[np.ndarray] = None
    reference_data: Optional[np.ndarray] = None
    calibration_paths: dict[str, str] = field(default_factory=dict)
    detections_by_frame: dict[int, list[Detection]] = field(default_factory=dict)
    aligned_cache: dict[int, np.ndarray] = field(default_factory=dict)
    residual_cache: dict[int, np.ndarray] = field(default_factory=dict)

    @property
    def potential_discoveries(self) -> list[Tracklet]:
        return [item for item in self.tracklets if item.classification == "unknown_candidate"]

    @property
    def borderline_review(self) -> list[Tracklet]:
        return [item for item in self.tracklets if item.classification == "review_candidate"]

    @property
    def known_recovered(self) -> list[KnownObjectPrediction]:
        return [item for item in self.known_objects if item.visible]

    @property
    def known_missed(self) -> list[KnownObjectPrediction]:
        return [item for item in self.known_objects if not item.visible]
