import math
from datetime import datetime, timedelta, timezone


ADES_DRAFT_FIELDS = [
    "permID",
    "provID",
    "trkSub",
    "mode",
    "stn",
    "obsTime",
    "ra",
    "dec",
    "rmsRA",
    "rmsDec",
    "astCat",
    "mag",
    "band",
    "photCat",
    "photAp",
    "logSNR",
    "seeing",
    "exp",
    "notes",
]


def render_ades_psv_draft(tracklets, sequence_frames, observatory_code: str = "") -> str:
    """Render an ADES-like PSV draft for review, not direct MPC submission."""
    lines = [
        "# version=2022",
        "# software",
        "! name Voyager Alpha",
        "# observatory",
        f"! mpcCode {observatory_code or 'UNKNOWN'}",
        "# submitter",
        "! name Voyager Alpha draft export",
        "# data",
        "|".join(ADES_DRAFT_FIELDS),
    ]

    for tracklet in tracklets:
        summary = tracklet.to_summary()
        for detection in tracklet.detections:
            frame = sequence_frames[detection.frame_index] if detection.frame_index < len(sequence_frames) else None
            obs_time = _frame_midpoint_utc(frame)
            uncertainty = detection.position_uncertainty_arcsec
            if uncertainty is None:
                uncertainty = tracklet.fit_rms_arcsec
            image_scale = frame.camera.image_scale_arcsec_px if frame else None
            seeing = detection.fwhm_px * image_scale if image_scale else None
            log_snr = math.log10(detection.snr) if detection.snr > 0 else None
            lines.append(
                "|".join(
                    [
                        "",
                        "",
                        summary["id"],
                        "CCD",
                        observatory_code,
                        obs_time,
                        "" if detection.ra is None else f"{detection.ra:.8f}",
                        "" if detection.dec is None else f"{detection.dec:.8f}",
                        "" if uncertainty is None else f"{uncertainty:.3f}",
                        "" if uncertainty is None else f"{uncertainty:.3f}",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "" if log_snr is None else f"{log_snr:.3f}",
                        "" if seeing is None else f"{seeing:.3f}",
                        "" if frame is None or frame.exposure_seconds is None else f"{frame.exposure_seconds:.3f}",
                        "DRAFT_NOT_MPC_READY",
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _frame_midpoint_utc(frame) -> str:
    if frame is None:
        return ""
    if frame.midpoint_utc:
        return frame.midpoint_utc
    date_obs = frame.date_obs
    if not date_obs:
        return ""
    try:
        value = datetime.fromisoformat(date_obs.replace("Z", "+00:00"))
        if value.tzinfo is None:
            value = value.replace(tzinfo=timezone.utc)
        value += timedelta(seconds=max(float(frame.exposure_seconds or 0.0), 0.0) / 2.0)
        return value.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")
    except ValueError:
        return date_obs
