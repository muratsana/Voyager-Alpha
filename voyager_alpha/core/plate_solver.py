from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from astropy.io import fits
from astropy.wcs import WCS


DEFAULT_ASTAP_PATHS = (
    Path(r"C:\Program Files\astap\astap.exe"),
    Path(r"C:\Program Files\astap\astap_cli.exe"),
    Path(r"C:\Program Files (x86)\astap\astap.exe"),
)


@dataclass
class PlateSolveResult:
    success: bool
    solved_file: str | None = None
    message: str = ""
    return_code: int | None = None
    header: fits.Header | None = None
    wcs_path: str | None = None
    solution_status: str = ""
    warning: str = ""
    log_path: str | None = None
    command: list[str] = field(default_factory=list)


class AstapPlateSolver:
    """Run ASTAP against a temporary FITS copy and return the solved WCS header."""

    def __init__(
        self,
        executable: str | None = None,
        timeout_seconds: int = 180,
        search_radius_deg: float = 30.0,
        downsample: int = 0,
    ):
        self.executable = executable or _discover_astap()
        self.timeout_seconds = int(timeout_seconds)
        self.search_radius_deg = float(search_radius_deg)
        self.downsample = max(0, int(downsample))

    def is_available(self) -> bool:
        return bool(self.executable and Path(self.executable).is_file())

    def solve(
        self,
        fits_path: str,
        *,
        fov_deg: float | None = None,
        ra_hours: float | None = None,
        south_pole_distance_deg: float | None = None,
    ) -> PlateSolveResult:
        if not self.is_available():
            return PlateSolveResult(False, message="ASTAP executable not found.")

        source = Path(fits_path)
        if not source.is_file():
            return PlateSolveResult(False, message=f"FITS file not found: {fits_path}")

        with tempfile.TemporaryDirectory(prefix="voyager-alpha-astap-") as temp_dir:
            working = Path(temp_dir) / source.name
            shutil.copy2(source, working)
            command = [
                str(self.executable),
                "-f",
                str(working),
                "-r",
                f"{self.search_radius_deg:g}",
                "-z",
                str(self.downsample),
                "-wcs",
                "-update",
                "-log",
            ]
            if fov_deg is not None:
                command.extend(["-fov", f"{float(fov_deg):.8g}"])
            if ra_hours is not None:
                command.extend(["-ra", f"{float(ra_hours):.10g}"])
            if south_pole_distance_deg is not None:
                command.extend(["-spd", f"{float(south_pole_distance_deg):.10g}"])

            try:
                completed = subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    timeout=self.timeout_seconds,
                    check=False,
                    cwd=temp_dir,
                    creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
                )
            except subprocess.TimeoutExpired:
                return PlateSolveResult(False, message="ASTAP solve timed out.", command=command)
            except OSError as exc:
                return PlateSolveResult(False, message=str(exc), command=command)

            output = "\n".join(part for part in (completed.stdout, completed.stderr) if part).strip()
            header, wcs_file = _read_astap_solution(working)
            status, warning, log_file = _read_astap_status(working, Path(temp_dir))
            success = (
                completed.returncode == 0
                and header is not None
                and _header_has_wcs(header)
                and status.upper() != "F"
            )
            if not success and not output:
                output = warning or "ASTAP did not produce a usable .wcs solution."
            return PlateSolveResult(
                success=success,
                solved_file=str(source),
                message=output,
                return_code=completed.returncode,
                header=header,
                wcs_path=str(wcs_file) if wcs_file else None,
                solution_status=status,
                warning=warning,
                log_path=str(log_file) if log_file else None,
                command=command,
            )


def merge_wcs_header(original: fits.Header, solved: fits.Header) -> fits.Header:
    """Return an in-memory header with solved astrometric cards merged in."""

    merged = original.copy()
    for card in solved.cards:
        key = card.keyword
        if key in {"", "COMMENT", "HISTORY", "END"}:
            continue
        merged[key] = (card.value, card.comment)
    return merged


def _read_astap_solution(working: Path) -> tuple[fits.Header | None, Path | None]:
    candidates = [working.with_suffix(".wcs"), Path(str(working) + ".wcs")]
    for candidate in candidates:
        if not candidate.is_file():
            continue
        try:
            header = fits.Header.fromtextfile(candidate)
            if _header_has_wcs(header):
                return header, candidate
        except Exception:
            try:
                raw = candidate.read_text(encoding="ascii", errors="ignore")
                header = fits.Header.fromstring(raw, sep="\n")
                if _header_has_wcs(header):
                    return header, candidate
            except Exception:
                pass
    try:
        header = fits.getheader(working, 0)
        if _header_has_wcs(header):
            return header, None
    except Exception:
        pass
    return None, None


def _read_astap_status(working: Path, temp_dir: Path) -> tuple[str, str, Path | None]:
    values: dict[str, str] = {}
    ini_candidates = (working.with_suffix(".ini"), Path(str(working) + ".ini"), temp_dir / "astap.ini")
    for candidate in ini_candidates:
        if not candidate.is_file():
            continue
        for raw_line in candidate.read_text(encoding="utf-8", errors="ignore").splitlines():
            if "=" not in raw_line:
                continue
            key, value = raw_line.split("=", 1)
            values[key.strip().upper()] = value.strip().strip("'\"")
        if values:
            break

    log_candidates = (working.with_suffix(".log"), Path(str(working) + ".log"), temp_dir / "astap.log")
    log_file = next((candidate for candidate in log_candidates if candidate.is_file()), None)
    warning = values.get("WARNING", "")
    if not warning and log_file is not None:
        log_text = log_file.read_text(encoding="utf-8", errors="ignore")
        warning_line = next((line.strip() for line in log_text.splitlines() if "warning" in line.lower()), "")
        warning = warning_line
    return values.get("PLTSOLVD", ""), warning, log_file


def _header_has_wcs(header: fits.Header) -> bool:
    try:
        return bool(WCS(header).has_celestial)
    except Exception:
        return False


def _looks_solved(output: str) -> bool:
    text = output.lower()
    if not text:
        return True
    failure_markers = ("not solved", "no solution", "failed", "error")
    if any(marker in text for marker in failure_markers):
        return False
    return any(marker in text for marker in ("solution found", "solved", "wcs", "solution"))


def _discover_astap() -> str | None:
    from_path = shutil.which("astap.exe") or shutil.which("astap_cli.exe") or shutil.which("astap")
    if from_path:
        return from_path
    for candidate in DEFAULT_ASTAP_PATHS:
        if candidate.is_file():
            return str(candidate)
    return None
