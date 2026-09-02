from html import escape


def render_tracklet_html_report(sequence_frames, tracklets) -> str:
    rows = []
    for tracklet in tracklets:
        summary = tracklet.to_summary()
        known_match = summary["known_match"]["name"] if summary["known_match"] else "No match"
        motion = summary["motion_arcsec_per_min"]
        motion_text = f"{motion:.3f} arcsec/min" if motion is not None else f"{summary['motion_px_per_frame']:.3f} px/frame"
        rows.append(
            "<tr>"
            f"<td>{escape(summary['id'])}</td>"
            f"<td>{summary['frames']}</td>"
            f"<td>{escape(motion_text)}</td>"
            f"<td>{summary['position_angle_deg']:.2f}</td>"
            f"<td>{summary['fit_rms_px']:.3f}</td>"
            f"<td>{summary['median_snr']:.2f}</td>"
            f"<td>{summary['confidence']:.3f}</td>"
            f"<td>{escape(known_match)}</td>"
            f"<td>{escape(summary['review_status'])}</td>"
            "</tr>"
        )

    flags = sorted({flag for frame in sequence_frames for flag in frame.quality_flags})
    wcs_ok = sum(1 for frame in sequence_frames if frame.has_wcs)
    time_ok = sum(1 for frame in sequence_frames if frame.midpoint_jd is not None)
    max_shift = max(
        [abs(frame.registration_dx) + abs(frame.registration_dy) for frame in sequence_frames],
        default=0,
    )

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Voyager Alpha Asteroid Report</title>
  <style>
    body {{
      margin: 32px;
      background: #111417;
      color: #e8edf2;
      font-family: Arial, sans-serif;
      line-height: 1.45;
    }}
    h1, h2 {{ margin: 0 0 12px; }}
    .summary {{
      display: grid;
      grid-template-columns: repeat(4, minmax(140px, 1fr));
      gap: 10px;
      margin: 20px 0;
    }}
    .metric {{
      border: 1px solid #2d343b;
      border-radius: 6px;
      padding: 12px;
      background: #181d22;
    }}
    .label {{ color: #9aa6b2; font-size: 12px; }}
    .value {{ color: #ffffff; font-size: 20px; margin-top: 4px; }}
    table {{
      width: 100%;
      border-collapse: collapse;
      background: #151a1f;
    }}
    th, td {{
      border-bottom: 1px solid #2d343b;
      padding: 8px 10px;
      text-align: left;
      font-size: 13px;
    }}
    th {{ color: #9edff2; background: #1b2228; }}
    .warning {{ color: #f2c94c; }}
    .note {{ color: #aeb8c2; max-width: 980px; }}
  </style>
</head>
<body>
  <h1>Voyager Alpha Asteroid Report</h1>
  <p class="note">This report lists multi-frame tracklet candidates. Unknown status is not an MPC discovery claim; it requires human review and astrometric validation.</p>
  <section class="summary">
    <div class="metric"><div class="label">Frames</div><div class="value">{len(sequence_frames)}</div></div>
    <div class="metric"><div class="label">WCS OK</div><div class="value">{wcs_ok}/{len(sequence_frames)}</div></div>
    <div class="metric"><div class="label">Time OK</div><div class="value">{time_ok}/{len(sequence_frames)}</div></div>
    <div class="metric"><div class="label">Max Registration Shift</div><div class="value">{max_shift} px</div></div>
  </section>
  <p class="warning">Quality flags: {escape(", ".join(flags) if flags else "OK")}</p>
  <h2>Tracklet Evidence</h2>
  <table>
    <thead>
      <tr>
        <th>ID</th><th>Frames</th><th>Motion</th><th>PA deg</th><th>Fit RMS px</th>
        <th>Median SNR</th><th>Confidence</th><th>SkyBoT</th><th>Review</th>
      </tr>
    </thead>
    <tbody>
      {"".join(rows) if rows else "<tr><td colspan='9'>No tracklets detected.</td></tr>"}
    </tbody>
  </table>
</body>
</html>
"""
