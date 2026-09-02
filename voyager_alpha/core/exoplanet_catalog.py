from __future__ import annotations

import csv
import io
import os
import sqlite3
import urllib.parse
import urllib.request
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Callable, Iterable

import numpy as np


TAP_URL = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"


@dataclass(frozen=True)
class CatalogSource:
    key: str
    display_name: str
    query: str
    parser: Callable[[dict[str, str]], "TransitCatalogObject | None"]


@dataclass(frozen=True)
class TransitCatalogObject:
    source_key: str
    source_id: str
    name: str
    host_name: str
    disposition: str
    ra_deg: float
    dec_deg: float
    period_days: float | None = None
    epoch_bjd: float | None = None
    duration_hours: float | None = None
    depth_ppm: float | None = None
    magnitude: float | None = None
    source_updated: str = ""
    separation_arcsec: float | None = None


@dataclass(frozen=True)
class TransitPrediction:
    center_bjd: float
    offset_minutes: float
    overlaps_observation: bool


CATALOG_SOURCES = (
    CatalogSource(
        "nea_confirmed",
        "NASA Confirmed Transits",
        "select pl_name,hostname,ra,dec,pl_orbper,pl_tranmid,pl_trandur,pl_trandep,disc_facility "
        "from pscomppars where tran_flag=1",
        lambda row: _parse_confirmed(row),
    ),
    CatalogSource(
        "tess_toi",
        "TESS TOI / ExoFOP",
        "select toidisplay,tid,toi,ra,dec,tfopwg_disp,pl_orbper,pl_tranmid,pl_trandurh,"
        "pl_trandep,st_tmag,rowupdate from toi",
        lambda row: _parse_toi(row),
    ),
    CatalogSource(
        "kepler_koi",
        "Kepler KOI",
        "select kepoi_name,kepler_name,kepid,ra,dec,koi_disposition,koi_period,koi_time0bk,"
        "koi_duration,koi_depth,koi_kepmag,koi_vet_date from cumulative",
        lambda row: _parse_koi(row),
    ),
    CatalogSource(
        "k2_candidates",
        "K2 Planets & Candidates",
        "select pl_name,epic_candname,k2_name,hostname,ra,dec,disposition,pl_orbper,pl_tranmid,"
        "pl_trandur,pl_trandep,sy_kepmag,rowupdate from k2pandc where default_flag=1",
        lambda row: _parse_k2(row),
    ),
)


class ExoplanetCatalog:
    def __init__(self, database_path: str | Path | None = None):
        self.database_path = Path(database_path) if database_path else default_catalog_path()
        self.database_path.parent.mkdir(parents=True, exist_ok=True)
        self._initialize()

    def update_all(
        self,
        *,
        progress_callback: Callable[[int, str, str], None] | None = None,
        stop_callback: Callable[[], bool] | None = None,
        downloader: Callable[[CatalogSource], str] | None = None,
    ) -> dict[str, int]:
        fetch = downloader or _download_source
        counts: dict[str, int] = {}
        total = len(CATALOG_SOURCES)
        for index, source in enumerate(CATALOG_SOURCES):
            if stop_callback and stop_callback():
                raise RuntimeError("Katalog güncellemesi durduruldu.")
            start = int(index / total * 100)
            if progress_callback:
                progress_callback(start, source.key, f"{source.display_name} indiriliyor")
            try:
                payload = fetch(source)
                rows = list(_parse_payload(source, payload))
                if not rows:
                    raise ValueError("Kaynak sıfır geçerli kayıt döndürdü.")
                stored_count = self.replace_source(source, rows)
                counts[source.key] = stored_count
                if progress_callback:
                    progress_callback(
                        int((index + 1) / total * 100),
                        source.key,
                        f"{source.display_name}: {stored_count:,} kayıt",
                    )
            except Exception as exc:
                self.mark_source_error(source, str(exc))
                raise RuntimeError(f"{source.display_name} güncellenemedi: {exc}") from exc
        return counts

    def replace_source(self, source: CatalogSource, rows: Iterable[TransitCatalogObject]):
        unique = {(item.source_key, item.source_id): item for item in rows}
        objects = list(unique.values())
        updated = datetime.now(timezone.utc).isoformat(timespec="seconds")
        with self._connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            connection.execute("DELETE FROM transit_objects WHERE source_key = ?", (source.key,))
            connection.executemany(
                """
                INSERT INTO transit_objects (
                    source_key, source_id, name, host_name, disposition, ra_deg, dec_deg,
                    period_days, epoch_bjd, duration_hours, depth_ppm, magnitude, source_updated
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    (
                        item.source_key,
                        item.source_id,
                        item.name,
                        item.host_name,
                        item.disposition,
                        item.ra_deg,
                        item.dec_deg,
                        item.period_days,
                        item.epoch_bjd,
                        item.duration_hours,
                        item.depth_ppm,
                        item.magnitude,
                        item.source_updated,
                    )
                    for item in objects
                ],
            )
            connection.execute(
                """
                INSERT INTO catalog_sources (source_key, display_name, record_count, updated_utc, status, error)
                VALUES (?, ?, ?, ?, 'ok', '')
                ON CONFLICT(source_key) DO UPDATE SET
                    display_name=excluded.display_name,
                    record_count=excluded.record_count,
                    updated_utc=excluded.updated_utc,
                    status='ok',
                    error=''
                """,
                (source.key, source.display_name, len(objects), updated),
            )
            connection.commit()
        return len(objects)

    def mark_source_error(self, source: CatalogSource, message: str):
        with self._connect() as connection:
            connection.execute(
                """
                INSERT INTO catalog_sources (source_key, display_name, record_count, updated_utc, status, error)
                VALUES (?, ?, 0, '', 'error', ?)
                ON CONFLICT(source_key) DO UPDATE SET status='error', error=excluded.error
                """,
                (source.key, source.display_name, message[:500]),
            )

    def source_status(self) -> list[dict[str, str | int]]:
        with self._connect() as connection:
            stored = {
                row["source_key"]: dict(row)
                for row in connection.execute("SELECT * FROM catalog_sources").fetchall()
            }
        result = []
        for source in CATALOG_SOURCES:
            row = stored.get(source.key, {})
            result.append(
                {
                    "source_key": source.key,
                    "display_name": source.display_name,
                    "record_count": int(row.get("record_count", 0)),
                    "updated_utc": str(row.get("updated_utc", "")),
                    "status": str(row.get("status", "missing")),
                    "error": str(row.get("error", "")),
                }
            )
        return result

    def is_stale(self, max_age_days: int = 7) -> bool:
        limit = datetime.now(timezone.utc) - timedelta(days=max_age_days)
        for row in self.source_status():
            if row["status"] != "ok" or not row["updated_utc"]:
                return True
            try:
                updated = datetime.fromisoformat(str(row["updated_utc"]))
            except ValueError:
                return True
            if updated < limit:
                return True
        return False

    def total_records(self) -> int:
        with self._connect() as connection:
            row = connection.execute("SELECT COUNT(*) AS n FROM transit_objects").fetchone()
        return int(row["n"])

    def cone_search(
        self,
        ra_deg: float,
        dec_deg: float,
        *,
        radius_arcmin: float = 5.0,
        limit: int = 30,
    ) -> list[TransitCatalogObject]:
        radius_deg = max(float(radius_arcmin), 0.01) / 60.0
        with self._connect() as connection:
            rows = connection.execute(
                "SELECT * FROM transit_objects WHERE dec_deg BETWEEN ? AND ?",
                (float(dec_deg) - radius_deg, float(dec_deg) + radius_deg),
            ).fetchall()
        matches = []
        for row in rows:
            separation = angular_separation_arcsec(float(ra_deg), float(dec_deg), row["ra_deg"], row["dec_deg"])
            if separation > radius_deg * 3600.0:
                continue
            matches.append(
                TransitCatalogObject(
                    source_key=row["source_key"],
                    source_id=row["source_id"],
                    name=row["name"],
                    host_name=row["host_name"],
                    disposition=row["disposition"],
                    ra_deg=row["ra_deg"],
                    dec_deg=row["dec_deg"],
                    period_days=row["period_days"],
                    epoch_bjd=row["epoch_bjd"],
                    duration_hours=row["duration_hours"],
                    depth_ppm=row["depth_ppm"],
                    magnitude=row["magnitude"],
                    source_updated=row["source_updated"],
                    separation_arcsec=separation,
                )
            )
        priority = {"confirmed": 0, "candidate": 1, "unverified": 2, "false_positive": 3}
        matches.sort(key=lambda item: (item.separation_arcsec or 0.0, priority.get(item.disposition, 4), item.name))
        return matches[: max(1, int(limit))]

    def _initialize(self):
        with self._connect() as connection:
            connection.executescript(
                """
                CREATE TABLE IF NOT EXISTS catalog_sources (
                    source_key TEXT PRIMARY KEY,
                    display_name TEXT NOT NULL,
                    record_count INTEGER NOT NULL DEFAULT 0,
                    updated_utc TEXT NOT NULL DEFAULT '',
                    status TEXT NOT NULL DEFAULT 'missing',
                    error TEXT NOT NULL DEFAULT ''
                );
                CREATE TABLE IF NOT EXISTS transit_objects (
                    id INTEGER PRIMARY KEY,
                    source_key TEXT NOT NULL,
                    source_id TEXT NOT NULL,
                    name TEXT NOT NULL,
                    host_name TEXT NOT NULL DEFAULT '',
                    disposition TEXT NOT NULL,
                    ra_deg REAL NOT NULL,
                    dec_deg REAL NOT NULL,
                    period_days REAL,
                    epoch_bjd REAL,
                    duration_hours REAL,
                    depth_ppm REAL,
                    magnitude REAL,
                    source_updated TEXT NOT NULL DEFAULT '',
                    UNIQUE(source_key, source_id)
                );
                CREATE INDEX IF NOT EXISTS idx_transit_objects_dec ON transit_objects(dec_deg);
                CREATE INDEX IF NOT EXISTS idx_transit_objects_position ON transit_objects(ra_deg, dec_deg);
                """
            )

    @contextmanager
    def _connect(self):
        connection = sqlite3.connect(self.database_path, timeout=30.0)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA journal_mode=WAL")
        connection.execute("PRAGMA foreign_keys=ON")
        try:
            yield connection
        finally:
            connection.close()


def default_catalog_path() -> Path:
    local = os.environ.get("LOCALAPPDATA")
    root = Path(local) if local else Path.home() / ".voyager-alpha"
    return root / "Astrohub" / "SkySearch" / "catalogs" / "exoplanets.sqlite3"


def predict_nearest_transit(
    item: TransitCatalogObject,
    observation_start_jd: float,
    observation_end_jd: float,
) -> TransitPrediction | None:
    if not item.period_days or not item.epoch_bjd or item.period_days <= 0:
        return None
    start = float(min(observation_start_jd, observation_end_jd))
    end = float(max(observation_start_jd, observation_end_jd))
    midpoint = (start + end) / 2.0
    cycle = int(round((midpoint - item.epoch_bjd) / item.period_days))
    center = float(item.epoch_bjd + cycle * item.period_days)
    half_duration_days = max(float(item.duration_hours or 0.0), 0.0) / 48.0
    overlaps = bool(start - half_duration_days <= center <= end + half_duration_days)
    return TransitPrediction(center, (center - midpoint) * 1440.0, overlaps)


def angular_separation_arcsec(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    ra1_rad, dec1_rad, ra2_rad, dec2_rad = np.deg2rad([ra1, dec1, ra2, dec2])
    cosine = np.sin(dec1_rad) * np.sin(dec2_rad) + np.cos(dec1_rad) * np.cos(dec2_rad) * np.cos(ra1_rad - ra2_rad)
    return float(np.rad2deg(np.arccos(np.clip(cosine, -1.0, 1.0))) * 3600.0)


def _download_source(source: CatalogSource) -> str:
    url = TAP_URL + "?" + urllib.parse.urlencode({"query": source.query, "format": "csv"})
    request = urllib.request.Request(url, headers={"User-Agent": "VoyagerAlpha/1.0 (transit catalog updater)"})
    with urllib.request.urlopen(request, timeout=120) as response:
        payload = response.read()
    text = payload.decode("utf-8-sig", errors="strict")
    if text.lstrip().lower().startswith(("<!doctype", "<html")):
        raise ValueError("Katalog servisi CSV yerine beklenmeyen içerik döndürdü.")
    return text


def _parse_payload(source: CatalogSource, payload: str):
    reader = csv.DictReader(io.StringIO(payload))
    if not reader.fieldnames:
        raise ValueError("CSV başlığı bulunamadı.")
    for row in reader:
        item = source.parser(row)
        if item is not None:
            yield item


def _parse_confirmed(row: dict[str, str]) -> TransitCatalogObject | None:
    return _object(
        "nea_confirmed",
        row.get("pl_name"),
        row.get("pl_name"),
        row.get("hostname"),
        "confirmed",
        row,
        period="pl_orbper",
        epoch="pl_tranmid",
        duration="pl_trandur",
        depth="pl_trandep",
        depth_scale=10000.0,
    )


def _parse_toi(row: dict[str, str]) -> TransitCatalogObject | None:
    disposition = {
        "CP": "confirmed",
        "KP": "confirmed",
        "PC": "candidate",
        "APC": "candidate",
        "FP": "false_positive",
        "FA": "false_positive",
    }.get((row.get("tfopwg_disp") or "").strip().upper(), "unverified")
    return _object(
        "tess_toi",
        row.get("toi") or row.get("tid"),
        row.get("toidisplay") or f"TOI-{row.get('toi', '')}",
        f"TIC {row.get('tid', '')}".strip(),
        disposition,
        row,
        period="pl_orbper",
        epoch="pl_tranmid",
        duration="pl_trandurh",
        depth="pl_trandep",
        magnitude="st_tmag",
        updated="rowupdate",
    )


def _parse_koi(row: dict[str, str]) -> TransitCatalogObject | None:
    disposition = {
        "CONFIRMED": "confirmed",
        "CANDIDATE": "candidate",
        "FALSE POSITIVE": "false_positive",
        "NOT DISPOSITIONED": "unverified",
    }.get((row.get("koi_disposition") or "").strip().upper(), "unverified")
    epoch = _number(row.get("koi_time0bk"))
    if epoch is not None:
        epoch += 2454833.0
    item = _object(
        "kepler_koi",
        row.get("kepoi_name"),
        row.get("kepler_name") or row.get("kepoi_name"),
        f"KIC {row.get('kepid', '')}".strip(),
        disposition,
        row,
        period="koi_period",
        duration="koi_duration",
        depth="koi_depth",
        magnitude="koi_kepmag",
        updated="koi_vet_date",
    )
    if item is None:
        return None
    return TransitCatalogObject(**{**item.__dict__, "epoch_bjd": epoch})


def _parse_k2(row: dict[str, str]) -> TransitCatalogObject | None:
    disposition = {
        "CONFIRMED": "confirmed",
        "CANDIDATE": "candidate",
        "FALSE POSITIVE": "false_positive",
        "REFUTED": "false_positive",
    }.get((row.get("disposition") or "").strip().upper(), "unverified")
    name = row.get("k2_name") or row.get("pl_name") or row.get("epic_candname")
    return _object(
        "k2_candidates",
        row.get("epic_candname") or name,
        name,
        row.get("hostname"),
        disposition,
        row,
        period="pl_orbper",
        epoch="pl_tranmid",
        duration="pl_trandur",
        depth="pl_trandep",
        depth_scale=10000.0,
        magnitude="sy_kepmag",
        updated="rowupdate",
    )


def _object(
    source_key: str,
    source_id,
    name,
    host_name,
    disposition: str,
    row: dict[str, str],
    *,
    period: str,
    epoch: str | None = None,
    duration: str,
    depth: str,
    depth_scale: float = 1.0,
    magnitude: str | None = None,
    updated: str | None = None,
) -> TransitCatalogObject | None:
    ra = _number(row.get("ra"))
    dec = _number(row.get("dec"))
    identifier = str(source_id or "").strip()
    object_name = str(name or identifier).strip()
    if not identifier or not object_name or ra is None or dec is None:
        return None
    depth_value = _number(row.get(depth))
    return TransitCatalogObject(
        source_key=source_key,
        source_id=identifier,
        name=object_name,
        host_name=str(host_name or "").strip(),
        disposition=disposition,
        ra_deg=ra,
        dec_deg=dec,
        period_days=_number(row.get(period)),
        epoch_bjd=_number(row.get(epoch)) if epoch else None,
        duration_hours=_number(row.get(duration)),
        depth_ppm=depth_value * depth_scale if depth_value is not None else None,
        magnitude=_number(row.get(magnitude)) if magnitude else None,
        source_updated=str(row.get(updated) or "").strip() if updated else "",
    )


def _number(value) -> float | None:
    if value is None or str(value).strip() in {"", "null", "None", "nan"}:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if np.isfinite(result) else None
