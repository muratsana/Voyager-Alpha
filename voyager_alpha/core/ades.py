from datetime import datetime, timezone


ADES_DRAFT_FIELDS = [
    "permID",
    "provID",
    "trkSub",
    "obsTime",
    "ra",
    "dec",
    "rmsRA",
    "rmsDec",
    "mag",
    "band",
    "stn",
    "notes",
]


def render_ades_psv_draft(tracklets, sequence_frames, observatory_code: str = "") -> str:
    """Render an ADES-like PSV draft for review, not direct MPC submission."""
    lines = [
        "# version=2022",
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
            obs_time = _date_obs_to_utc(frame.date_obs if frame else None)
            lines.append(
                "|".join(
                    [
                        "",
                        "",
                        summary["id"],
                        obs_time,
                        "" if detection.ra is None else f"{detection.ra:.8f}",
                        "" if detection.dec is None else f"{detection.dec:.8f}",
                        "",
                        "",
                        "",
                        "",
                        observatory_code,
                        "DRAFT_NOT_MPC_READY",
                    ]
                )
            )
    return "\n".join(lines) + "\n"


def _date_obs_to_utc(date_obs: str | None) -> str:
    if not date_obs:
        return ""
    try:
        value = datetime.fromisoformat(date_obs.replace("Z", "+00:00"))
        if value.tzinfo is None:
            value = value.replace(tzinfo=timezone.utc)
        return value.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")
    except ValueError:
        return date_obs
