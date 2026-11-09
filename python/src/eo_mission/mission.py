"""Timed accesses for the selected thesis constellation using Orekit dynamics.

This analysis is separate from the MATLAB design sweep and the older Python
tiling approximation. Its revisit statistics come from propagated point
accesses; they are finite-horizon observations, not a guarantee of future access.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize_scalar

from .orekit_runtime import initialize_runtime

EARTH_RADIUS_M = 6_378_137.0
EARTH_FLATTENING = 1 / 298.257223563
MU = 3.986004418e14
J2 = 1.082629989052e-3
J3 = -2.5327e-6
EVENT_COLUMNS = ["satellite", "target", "start_s", "end_s", "duration_s"]


@dataclass(frozen=True)
class SelectedMission:
    epoch: str = "2024-04-28T00:00:00Z"
    altitude_km: float = 690.0
    inclination_deg: float = 98.1472649550129
    eccentricity: float = 0.0010448573934173355
    argument_of_perigee_deg: float = 90.0
    satellites: int = 2
    ltan_h: float = 6.0
    acquired_swath_km: float = 240.0
    overlap_fraction: float = 0.05
    gsd_m: float = 80.0
    duration_days: float = 7.0
    step_s: float = 10.0
    latitude_min_deg: float = 25.0
    latitude_max_deg: float = 49.0
    longitude_min_deg: float = -125.0
    longitude_max_deg: float = -66.0
    target_latitude_step_deg: float = 2.0
    target_longitude_step_deg: float = 3.0
    station_name: str = "NOAA/Fairbanks geographic baseline"
    station_latitude_deg: float = 64.84
    station_longitude_deg: float = -147.712
    station_altitude_m: float = 0.0
    minimum_elevation_deg: float = 10.0
    require_daylight: bool = True
    reference_rate_mbps: float = 10.0

    def validate(self) -> None:
        numeric = [v for v in asdict(self).values() if isinstance(v, (int, float))]
        if not all(math.isfinite(v) for v in numeric):
            raise ValueError("Mission numbers must be finite.")
        if self.altitude_km <= 0 or not 0 <= self.eccentricity < 0.1:
            raise ValueError("The zonal model requires positive altitude and eccentricity below 0.1.")
        if self.satellites < 1 or int(self.satellites) != self.satellites:
            raise ValueError("satellites must be a positive integer.")
        if not 0 <= self.inclination_deg <= 180 or not 0 <= self.ltan_h < 24:
            raise ValueError("Inclination must be in [0, 180] degrees and LTAN in [0, 24) hours.")
        if self.duration_days <= 0 or not 0 < self.step_s <= 30:
            raise ValueError("Use a positive horizon and a sampling step in (0, 30] seconds.")
        if self.step_s > self.duration_days * 86400:
            raise ValueError("The sampling step exceeds the simulation horizon.")
        if self.acquired_swath_km <= 0 or not 0 <= self.overlap_fraction < 1 or self.gsd_m <= 0:
            raise ValueError("Swath and GSD must be positive; overlap must be in [0, 1).")
        if not -90 <= self.latitude_min_deg < self.latitude_max_deg <= 90:
            raise ValueError("Target latitude limits are invalid.")
        if not -180 <= self.longitude_min_deg < self.longitude_max_deg <= 180:
            raise ValueError("Target longitude limits are invalid.")
        if self.target_latitude_step_deg <= 0 or self.target_longitude_step_deg <= 0:
            raise ValueError("Target-grid spacing must be positive.")
        if not -90 <= self.station_latitude_deg <= 90 or not -180 <= self.station_longitude_deg <= 180:
            raise ValueError("Ground-station coordinates are invalid.")
        if not 0 <= self.minimum_elevation_deg < 90 or self.reference_rate_mbps <= 0:
            raise ValueError("Minimum elevation must be in [0, 90); reference rate must be positive.")
        parsed = datetime.fromisoformat(self.epoch.replace("Z", "+00:00"))
        if parsed.utcoffset() is None or parsed.utcoffset().total_seconds() != 0:
            raise ValueError("Specify the epoch explicitly in UTC (ending in Z or +00:00).")


def load_selected_mission(
    config_path: str | Path, selected_design_path: str | Path | None = None,
) -> SelectedMission:
    """Load timing settings and take the adopted design from the MATLAB CSV."""
    path = Path(config_path)
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    mission = SelectedMission(**raw.get("mission", {}))
    selected = selected_design_path or raw.get("selected_design")
    if selected:
        selected_path = Path(selected)
        if not selected_path.is_absolute():
            selected_path = path.parent / selected_path
        table = pd.read_csv(selected_path)
        if len(table) != 1:
            raise ValueError("The selected-design CSV must contain exactly one design.")
        row = table.iloc[0]
        altitude = float(row["Altura_km"])
        inclination = float(row["Inclinacion_SSO_deg"])
        eccentricity = -J3 / (2 * J2) * EARTH_RADIUS_M / (EARTH_RADIUS_M + altitude * 1000) * math.sin(math.radians(inclination))
        mission = replace(
            mission, altitude_km=altitude, inclination_deg=inclination,
            eccentricity=eccentricity, satellites=int(row["Num_Satelites"]),
            acquired_swath_km=float(row["Swath_Adquirido_km"]),
            gsd_m=float(row["GSD"]), ltan_h=float(row["LTAN_h"]),
        )
    mission.validate()
    return mission


def _inclusive_grid(start: float, stop: float, step: float) -> np.ndarray:
    values = np.arange(start, stop, step, dtype=float)
    return np.append(values, stop) if len(values) == 0 or not np.isclose(values[-1], stop) else values


def target_grid(mission: SelectedMission) -> pd.DataFrame:
    """Sample the same rectangular study-region bounds used by the thesis."""
    mission.validate()
    rows = []
    for lat in _inclusive_grid(mission.latitude_min_deg, mission.latitude_max_deg, mission.target_latitude_step_deg):
        for lon in _inclusive_grid(mission.longitude_min_deg, mission.longitude_max_deg, mission.target_longitude_step_deg):
            rows.append({"target": f"P{len(rows) + 1:04d}", "latitude_deg": float(lat), "longitude_deg": float(lon)})
    return pd.DataFrame(rows)


def geodetic_position(latitude_deg: float, longitude_deg: float, altitude_m: float = 0) -> np.ndarray:
    """WGS84 geodetic coordinates to Earth-fixed metres."""
    lat, lon = np.radians([latitude_deg, longitude_deg])
    e2 = EARTH_FLATTENING * (2 - EARTH_FLATTENING)
    n = EARTH_RADIUS_M / math.sqrt(1 - e2 * math.sin(lat)**2)
    return np.array([(n + altitude_m) * math.cos(lat) * math.cos(lon),
                     (n + altitude_m) * math.cos(lat) * math.sin(lon),
                     (n * (1 - e2) + altitude_m) * math.sin(lat)])


def _up(latitude_deg: float | np.ndarray, longitude_deg: float | np.ndarray) -> np.ndarray:
    lat, lon = np.radians(latitude_deg), np.radians(longitude_deg)
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1)


def geodetic_coordinates(xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Earth-fixed metres to WGS84 geodetic degrees and ellipsoidal metres."""
    x, y, z = np.asarray(xyz, dtype=float).T
    a, f = EARTH_RADIUS_M, EARTH_FLATTENING
    b, e2 = a*(1 - f), f*(2 - f)
    ep2 = (a*a - b*b)/(b*b)
    p = np.hypot(x, y)
    theta = np.arctan2(z*a, p*b)
    latitude = np.arctan2(z + ep2*b*np.sin(theta)**3, p - e2*a*np.cos(theta)**3)
    for _ in range(3):
        n = a/np.sqrt(1 - e2*np.sin(latitude)**2)
        latitude = np.arctan2(z + e2*n*np.sin(latitude), p)
    height = p*np.cos(latitude) + z*np.sin(latitude) - a*np.sqrt(1 - e2*np.sin(latitude)**2)
    return np.degrees(latitude), np.degrees(np.arctan2(y, x)), height


def interval_bounds(times: np.ndarray, margins: np.ndarray) -> list[tuple[float, float]]:
    """Find positive-margin intervals, linearly interpolating sampled crossings.

    A crossing outside the horizon is clipped to its boundary. Sub-step grazing
    passes can be missed, so the configured sampling step is exported with every
    analysis and convergence must be checked before operational use.
    """
    times = np.asarray(times, dtype=float)
    margins = np.asarray(margins, dtype=float)
    if times.ndim != 1 or times.shape != margins.shape or len(times) < 2:
        raise ValueError("Times and margins must be matching 1-D arrays with at least two samples.")
    if not np.all(np.isfinite(times)) or not np.all(np.isfinite(margins)) or np.any(np.diff(times) <= 0):
        raise ValueError("Samples must be finite and times strictly increasing.")
    inside = margins >= 0
    transitions = np.flatnonzero(inside[1:] != inside[:-1])
    crossings = {int(i): float(times[i] - margins[i] * (times[i + 1] - times[i]) / (margins[i + 1] - margins[i])) for i in transitions}
    starts = ([float(times[0])] if inside[0] else []) + [crossings[int(i)] for i in transitions if not inside[i]]
    ends = [crossings[int(i)] for i in transitions if inside[i]] + ([float(times[-1])] if inside[-1] else [])
    return [(start, end) for start, end in zip(starts, ends) if end > start]


def merge_intervals(intervals: list[tuple[float, float]]) -> list[tuple[float, float]]:
    merged: list[tuple[float, float]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def refine_intervals(times: np.ndarray, margins: np.ndarray, evaluate, angular_gate: float) -> list[tuple[float, float]]:
    """Search closest approaches between samples and refine event roots.

    A near-threshold local maximum is searched even when both surrounding
    samples are outside the footprint. This recovers grazing passes that a
    simple sampled boolean mask misses. The angular gate covers the possible
    LEO line-of-sight motion between neighboring samples.
    """
    positive = np.flatnonzero(margins >= 0)
    peaks = np.flatnonzero((margins[1:-1] >= margins[:-2]) & (margins[1:-1] >= margins[2:])) + 1
    peaks = peaks[margins[peaks] >= -angular_gate]
    candidates = np.union1d(positive, peaks)
    windows = merge_intervals([(float(times[max(0, i - 1)]), float(times[min(len(times) - 1, i + 1)])) for i in candidates])
    intervals = []
    for left, right in windows:
        peak = minimize_scalar(lambda t: -evaluate(t), bounds=(left, right), method="bounded", options={"xatol": 0.001})
        maximum = max([(float(peak.x), float(-peak.fun)), (left, evaluate(left)), (right, evaluate(right))], key=lambda p: p[1])
        if maximum[1] <= 0:
            continue
        middle = maximum[0]
        start = left if evaluate(left) >= 0 else brentq(evaluate, left, middle, xtol=0.001)
        end = right if evaluate(right) >= 0 else brentq(evaluate, middle, right, xtol=0.001)
        if end > start:
            intervals.append((float(start), float(end)))
    return merge_intervals(intervals)


def target_metrics(targets: pd.DataFrame, accesses: pd.DataFrame, horizon_s: float) -> pd.DataFrame:
    """Combine both satellites and retain censored/no-access targets explicitly."""
    rows = []
    for target in targets.to_dict("records"):
        events = accesses[accesses["target"] == target["target"]]
        intervals = merge_intervals(list(zip(events.start_s, events.end_s)))
        gaps = [b[0] - a[1] for a, b in zip(intervals, intervals[1:])]
        left = intervals[0][0] if intervals else horizon_s
        right = horizon_s - intervals[-1][1] if intervals else horizon_s
        rows.append({**target, "access_count": len(intervals),
                     "observed_max_gap_days": max(gaps) / 86400 if gaps else None,
                     "max_gap_days": max([left, right, *gaps]) / 86400,
                     "left_censored_gap_days": left / 86400,
                     "right_censored_gap_days": right / 86400,
                     "never_observed": not bool(intervals)})
    return pd.DataFrame(rows)


def analyze_accesses(ephemeris: pd.DataFrame, sun_positions: np.ndarray, mission: SelectedMission, *, targets: pd.DataFrame | None = None, include_contacts: bool = True) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute clear-sky nadir-cone accesses and station-elevation contacts."""
    targets = target_grid(mission) if targets is None else targets
    access_rows, contact_rows = [], []
    effective_swath = mission.acquired_swath_km * (1 - mission.overlap_fraction)
    half_angle = math.atan(effective_swath / (2 * mission.altitude_km))
    station_xyz = geodetic_position(mission.station_latitude_deg, mission.station_longitude_deg, mission.station_altitude_m)
    station_up = _up(mission.station_latitude_deg, mission.station_longitude_deg)
    for satellite, frame in ephemeris.groupby("satellite", sort=True):
        frame = frame.sort_values("time_s")
        times = frame.time_s.to_numpy()
        xyz = frame[["x_m", "y_m", "z_m"]].to_numpy()
        if len(sun_positions) != len(frame):
            raise ValueError("Sun and satellite sample counts differ.")
        nadir = -_up(frame.latitude_deg.to_numpy(), frame.longitude_deg.to_numpy())
        position_at = CubicSpline(times, xyz, axis=0)
        sun_at = CubicSpline(times, sun_positions, axis=0)
        angular_gate = max(0.05, np.diff(times).max()*8500/(mission.altitude_km*1000)*0.75)
        station_los = xyz - station_xyz
        elevation = np.arcsin(np.clip(station_los @ station_up / np.linalg.norm(station_los, axis=1), -1, 1))
        def station_margin(t):
            los = position_at(t) - station_xyz
            return float(np.arcsin(np.clip(los @ station_up / np.linalg.norm(los), -1, 1)) - math.radians(mission.minimum_elevation_deg))
        if include_contacts:
            for start, end in refine_intervals(times, elevation - math.radians(mission.minimum_elevation_deg), station_margin, angular_gate):
                contact_rows.append({"satellite": int(satellite), "target": mission.station_name, "start_s": start, "end_s": end, "duration_s": end - start})
        for target in targets.to_dict("records"):
            point = geodetic_position(target["latitude_deg"], target["longitude_deg"])
            target_up = _up(target["latitude_deg"], target["longitude_deg"])
            los = point - xyz
            distance = np.linalg.norm(los, axis=1)
            off_axis = np.arccos(np.clip(np.einsum("ij,ij->i", los, nadir) / distance, -1, 1))
            margin = half_angle - off_axis
            horizon = np.arcsin(np.clip((-los) @ target_up / distance, -1, 1))
            margin = np.minimum(margin, horizon)
            if mission.require_daylight:
                to_sun = sun_positions - point
                daylight = np.arcsin(np.clip(to_sun @ target_up / np.linalg.norm(to_sun, axis=1), -1, 1))
                margin = np.minimum(margin, daylight)
            def target_margin(t):
                position = position_at(t)
                lat, lon, _ = geodetic_coordinates(position[None, :])
                direction = -_up(lat[0], lon[0])
                line = point - position
                distance = np.linalg.norm(line)
                off_axis = math.acos(float(np.clip(line @ direction/distance, -1, 1)))
                horizon = math.asin(float(np.clip((-line) @ target_up/distance, -1, 1)))
                result = min(half_angle - off_axis, horizon)
                if mission.require_daylight:
                    solar_line = sun_at(t) - point
                    result = min(result, math.asin(float(np.clip(solar_line @ target_up/np.linalg.norm(solar_line), -1, 1))))
                return result
            for start, end in refine_intervals(times, margin, target_margin, angular_gate):
                access_rows.append({"satellite": int(satellite), **target, "start_s": start, "end_s": end, "duration_s": end - start})
    accesses = pd.DataFrame(access_rows, columns=[*EVENT_COLUMNS, "latitude_deg", "longitude_deg"])
    contacts = pd.DataFrame(contact_rows, columns=EVENT_COLUMNS)
    metrics = target_metrics(targets, accesses, mission.duration_days * 86400)
    return accesses, contacts, metrics


def _vector(v) -> tuple[float, float, float]:
    return float(v.getX()), float(v.getY()), float(v.getZ())


def propagate_selected_mission(mission: SelectedMission, data_path: str | Path | None = None) -> tuple[pd.DataFrame, np.ndarray, dict]:
    """Propagate the frozen mean elements with Orekit's J2/J3 zonal model."""
    mission.validate()
    data = initialize_runtime(data_path)
    import orekit
    from org.orekit.bodies import CelestialBodyFactory, OneAxisEllipsoid
    from org.orekit.attitudes import FrameAlignedProvider
    from org.orekit.frames import FramesFactory
    from org.orekit.orbits import KeplerianOrbit, PositionAngleType
    from org.orekit.propagation import PropagationType
    from org.orekit.propagation.analytical import EcksteinHechlerPropagator
    from org.orekit.time import AbsoluteDate, TimeScalesFactory
    from org.orekit.utils import IERSConventions

    epoch = mission.epoch.replace("Z", "").replace("+00:00", "")
    start = AbsoluteDate(epoch, TimeScalesFactory.getUTC())
    inertial = FramesFactory.getEME2000()
    earth_frame = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
    earth = OneAxisEllipsoid(EARTH_RADIUS_M, EARTH_FLATTENING, earth_frame)
    sun = CelestialBodyFactory.getSun()
    solar_position = _vector(sun.getPosition(start, inertial))
    sun_ra = math.atan2(solar_position[1], solar_position[0])
    raan = sun_ra + math.radians((mission.ltan_h - 12) * 15)
    times = _inclusive_grid(0, mission.duration_days * 86400, mission.step_s)
    dates = [start.shiftedBy(float(t)) for t in times]
    sun_positions = np.array([_vector(sun.getPosition(date, earth_frame)) for date in dates])
    rows = []
    for satellite in range(mission.satellites):
        orbit = KeplerianOrbit(
            float(EARTH_RADIUS_M + mission.altitude_km * 1000), float(mission.eccentricity),
            math.radians(mission.inclination_deg), math.radians(mission.argument_of_perigee_deg),
            float(raan), 2 * math.pi * satellite / mission.satellites,
            PositionAngleType.TRUE, inertial, start, MU,
        )
        propagator = EcksteinHechlerPropagator(
            orbit, FrameAlignedProvider.of(inertial), 1.0, EARTH_RADIUS_M,
            MU, -J2, -J3, 0.0, 0.0, 0.0, PropagationType.MEAN,
        )
        for index, date in enumerate(dates):
            state = propagator.propagate(date)
            position = state.getPVCoordinates(earth_frame).getPosition()
            xyz = np.array(_vector(position))
            geodetic = earth.transform(position, earth_frame, date)
            solar_direction = sun_positions[index] / np.linalg.norm(sun_positions[index])
            behind = float(xyz @ solar_direction)
            perpendicular = np.linalg.norm(xyz - behind * solar_direction)
            sunlit = not (behind < 0 and perpendicular < EARTH_RADIUS_M)
            rows.append({"time_s": float(times[index]), "satellite": satellite + 1,
                         "x_m": xyz[0], "y_m": xyz[1], "z_m": xyz[2],
                         "latitude_deg": math.degrees(geodetic.getLatitude()),
                         "longitude_deg": math.degrees(geodetic.getLongitude()),
                         "altitude_km": geodetic.getAltitude() / 1000, "sunlit": sunlit})
    provenance = {"orekit_version": str(orekit.VERSION), "orekit_data_path": str(data),
                  "orekit_data_sha256": hashlib.sha256(data.read_bytes()).hexdigest() if data.is_file() else None,
                  "raan_epoch_deg": math.degrees(raan) % 360,
                  "elements": "mean Eckstein-Hechler; true-anomaly plane spacing",
                  "earth_frame": "ITRF, IERS 2010, simple EOP interpolation",
                  "gravity": "central attraction + unnormalized C20=-J2, C30=-J3; C40/C50/C60=0",
                  "mass": "unit mass; zonal-only dynamics are mass independent",
                  "eclipse_display": "cylindrical Earth shadow; acquisition uses target solar elevation",
                  "atmospheric_drag": False}
    return pd.DataFrame(rows), sun_positions, provenance


def _json_records(frame: pd.DataFrame) -> list[dict]:
    return json.loads(frame.to_json(orient="records"))


def add_cloud_analyses(summary: dict, acquisitions: pd.DataFrame, targets: pd.DataFrame,
                       settings: dict, output_dir: str | Path) -> None:
    """Apply shared daily Markov chains to the actual propagated intervals."""
    if not settings.get("enabled", True):
        return
    from .mission_clouds import apply_markov_clouds, cloud_hardware_floor
    mission = summary["mission"]
    parameters = dict(cloud_fraction=float(settings.get("cloud_fraction", 1/6)),
                      persistence=float(settings.get("persistence", 0.5)),
                      mc_samples=int(settings.get("mc_samples", 200)),
                      seed=int(settings.get("seed", 12345)),
                      epoch=mission["epoch"], requirement_days=float(settings.get("requirement_days", 7)))
    horizon = mission["duration_days"]*86400
    out = Path(output_dir)
    clouds, realizations = apply_markov_clouds(acquisitions, targets, horizon, **parameters)
    clouds["hardware_floor"] = cloud_hardware_floor(horizon, **parameters)
    summary["clouds"] = clouds
    realizations.to_csv(out / "cloud_realizations.csv", index=False)
    sensitivity = []
    for fraction in settings.get("sensitivity_fractions", []):
        alternate = dict(parameters, cloud_fraction=float(fraction))
        stats, records = apply_markov_clouds(acquisitions, targets, horizon, **alternate)
        stats["hardware_floor"] = cloud_hardware_floor(horizon, **alternate)
        sensitivity.append(stats)
        records.to_csv(out / f"cloud_realizations_p{float(fraction):.6f}.csv", index=False)
    summary["cloud_sensitivity"] = sensitivity
    (out / "cloud_summary.json").write_text(json.dumps({"baseline": clouds, "sensitivity": sensitivity}, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def run_selected_mission(config_path: str | Path, output_dir: str | Path, *, data_path: str | Path | None = None, html: bool = True, backend: str = "orekit") -> Path:
    mission = load_selected_mission(config_path)
    configuration = yaml.safe_load(Path(config_path).read_text(encoding="utf-8")) or {}
    if backend == "orekit":
        ephemeris, sun_positions, provenance = propagate_selected_mission(mission, data_path)
        model = "Orekit Eckstein-Hechler J2/J3"
        provenance.update({"backend": "orekit", "orekit_executed": True})
    elif backend == "scipy-reference":
        from .mission_reference import propagate_reference
        ephemeris, sun_positions, provenance = propagate_reference(mission)
        model = "SciPy numerical J2/J3 reference — Orekit not executed"
    else:
        raise ValueError("backend must be 'orekit' or 'scipy-reference'; no automatic fallback is used.")
    accesses, contacts, targets = analyze_accesses(ephemeris, sun_positions, mission)
    provenance["configuration_sha256"] = hashlib.sha256(Path(config_path).read_bytes()).hexdigest()
    selected_path = configuration.get("selected_design")
    if selected_path:
        selected_path = Path(selected_path)
        if not selected_path.is_absolute():
            selected_path = Path(config_path).parent / selected_path
        provenance["selected_design_path"] = str(selected_path.resolve())
        provenance["selected_design_sha256"] = hashlib.sha256(selected_path.read_bytes()).hexdigest()
    summary = {
        "mission": asdict(mission),
        "propagation": {"model": model, "duration_days": mission.duration_days, "step_s": mission.step_s},
        "metrics": {"target_count": len(targets), "access_count": len(accesses), "contact_count": len(contacts),
                    "unobserved_targets": int(targets.never_observed.sum()),
                    "worst_revisit_days": float(targets.max_gap_days.max()),
                    "total_contact_s": float(contacts.duration_s.sum())},
        "targets": _json_records(targets), "provenance": provenance,
        "definitions": {
            "revisit": "Largest finite-horizon sampled gap without clear-sky geometric access, including censored horizon boundaries. Sampling can miss grazing passes; this is not a long-term bound or weekly-coverage certification.",
            "target_grid": "Rectangular thesis study-region bounds; not a CONUS land-boundary polygon.",
            "footprint": "Nadir-aligned circular cone, with effective cross-track width acquired_swath*(1-overlap); an approximation to the pushbroom footprint.",
            "illumination": "Acquisition requires positive geometric target solar elevation when require_daylight is true; no atmospheric refraction or cloud simulation.",
            "contacts": "Pure geometry above the configured station elevation; no antenna scheduling or contracted service availability.",
            "resolution": "Cubic interpolation of propagated positions, closest-approach search and event-root refinement to 1 ms. Near-threshold sampled extrema are searched for grazing passes. Compare finer ephemeris steps before using mission conclusions.",
            "thesis_baseline": "Selected MATLAB design is input only. Orekit uses a different propagator; contact/revisit results need not equal SGP4 or analytical MATLAB values.",
        },
    }
    summary["provenance"]["result_role"] = "Adopted MATLAB design validation case; not an independently optimized solution"
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    ephemeris.to_csv(out / "ephemeris.csv", index=False)
    accesses.to_csv(out / "acquisitions.csv", index=False)
    contacts.to_csv(out / "contacts.csv", index=False)
    targets.to_csv(out / "target_metrics.csv", index=False)
    if configuration.get("clouds"):
        add_cloud_analyses(summary, accesses, targets, configuration["clouds"], out)
    summary["provenance"]["ephemeris_sha256"] = hashlib.sha256((out / "ephemeris.csv").read_bytes()).hexdigest()
    (out / "mission_summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if html:
        from .mission_visualization import build_mission_explorer
        return build_mission_explorer(ephemeris, accesses, contacts, summary, out / "mission_explorer.html")
    return out / "mission_summary.json"
