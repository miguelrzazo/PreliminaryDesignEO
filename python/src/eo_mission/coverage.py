"""SSO coverage and revisit time calculations.

Port of CoverageSSOAnaliticalfunction.m.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from typing import Sequence

from .config import PipelineConfig
from .tensor import LabeledTensor, as_array

_R_EARTH      = 6_371.0      # km
_MU_EARTH     = 398_600.4418 # km³/s²
_SIDEREAL_DAY = 86_164.0     # s
_US_WIDTH_KM  = 4_100.0      # continental coverage target (km)
_DAYLIGHT_FACTOR = 0.5       # only daytime passes are usable


def compute_coverage(
    altitudes_km: Sequence[float],
    swaths_km: Sequence[float],
    gsd: float,
    n_pix: int,
    fov_limit_deg: float,
    n_sat: int = 1,
    n_telescopes: int = 1,
    detector_type: int = 1,
    overlap: float = 0.0,
    cloud_cover: float = 0.30,
    cov_req: float = 7.0,
) -> pd.DataFrame:
    """Compute SSO coverage days for all (altitude, swath) combinations.

    Parameters
    ----------
    altitudes_km : sequence of float
        Orbital altitudes (km).
    swaths_km : sequence of float
        Swath widths to evaluate (km).
    gsd : float
        Ground sampling distance (m).
    n_pix : int
        Number of pixels across the detector.
    fov_limit_deg : float
        Maximum half-FOV angle for the telescope (degrees).
    n_sat : int
        Number of satellites in the constellation.
    n_telescopes : int
        Number of telescopes per satellite.
    detector_type : int
        Detector index (1-based); limits max swath to detector_type × GSD × n_pix.
    overlap : float
        Fractional swath overlap between adjacent passes.
    cloud_cover : float
        Cloud cover fraction (increases revisit time).
    cov_req : float
        Maximum acceptable revisit time (days); values exceeding this become NaN.

    Returns
    -------
    pd.DataFrame
        Coverage days indexed by altitude (km), columns by swath (km).
        NaN where constraints are violated or cov_req exceeded.
    """
    h_arr = np.asarray(altitudes_km, dtype=float)
    s_arr = np.asarray(swaths_km, dtype=float)

    max_det_swath = detector_type * gsd * n_pix / 1000.0  # km

    result = np.full((len(h_arr), len(s_arr)), np.nan)

    for i, h in enumerate(h_arr):
        v_orb  = np.sqrt(_MU_EARTH / (_R_EARTH + h))          # km/s
        T_orb  = 2 * np.pi * (_R_EARTH + h) / v_orb           # s
        max_swath_fov = 2 * h * np.tan(np.radians(fov_limit_deg / 2))

        for j, swath in enumerate(s_arr):
            if swath > max_swath_fov or swath > max_det_swath:
                continue

            eff_swath = swath * n_telescopes * (1 - overlap)
            passes_needed = np.ceil(_US_WIDTH_KM / eff_swath)
            t_between = T_orb / _DAYLIGHT_FACTOR
            base_days = (passes_needed * t_between / _SIDEREAL_DAY) / (1 - cloud_cover)
            final_days = base_days / n_sat
            if n_telescopes > 1 or detector_type > 1:
                final_days /= 0.9

            if final_days <= cov_req:
                result[i, j] = final_days

    index   = pd.Index(h_arr, name="altitude_km")
    columns = pd.Index(s_arr, name="swath_km")
    return pd.DataFrame(result, index=index, columns=columns)


class OrekitCoverageAdapter:
    """Coverage adapter that uses Orekit when available.

    The adapter keeps the public coverage API independent from Orekit startup
    details. In normal unit tests it uses the same deterministic equations as
    the analytical model. In Docker, ``initialize_orekit(require=True)`` starts
    the JVM and the orbital period is obtained from an Orekit Keplerian orbit.
    """

    def __init__(self, orekit_data_path: str | None = None, cache_dir: str | Path | None = None):
        self.orekit_data_path = orekit_data_path
        self.cache_dir = Path(cache_dir or ".cache/orekit-coverage")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.jvm_initialized = False

    def initialize_orekit(self, require: bool = False) -> bool:
        try:
            import orekit

            orekit.initVM()
            if self.orekit_data_path:
                from orekit.pyhelpers import setup_orekit_data

                setup_orekit_data(self.orekit_data_path)
            self.jvm_initialized = True
            return True
        except Exception:
            if require:
                raise
            self.jvm_initialized = False
            return False

    def _period_seconds(self, altitude_km: float) -> float:
        if self.jvm_initialized:
            try:
                from org.hipparchus.util import FastMath
                from org.orekit.frames import FramesFactory
                from org.orekit.orbits import KeplerianOrbit, PositionAngleType
                from org.orekit.time import AbsoluteDate
                from org.orekit.utils import Constants

                orbit = KeplerianOrbit(
                    Constants.WGS84_EARTH_EQUATORIAL_RADIUS + altitude_km * 1000.0,
                    0.001,
                    FastMath.toRadians(98.0),
                    0.0,
                    0.0,
                    0.0,
                    PositionAngleType.MEAN,
                    FramesFactory.getEME2000(),
                    AbsoluteDate.J2000_EPOCH,
                    Constants.EIGEN5C_EARTH_MU,
                )
                return float(orbit.getKeplerianPeriod())
            except Exception:
                pass
        radius_km = _R_EARTH + altitude_km
        return 2 * np.pi * np.sqrt(radius_km**3 / _MU_EARTH)

    def _cache_path(self, payload: dict) -> Path:
        encoded = json.dumps(payload, sort_keys=True).encode("utf-8")
        digest = hashlib.sha256(encoded).hexdigest()[:20]
        return self.cache_dir / f"{digest}.json"

    def revisit_days(
        self,
        altitude_km: float,
        swath_km: float,
        *,
        n_satellites: int,
        n_telescopes: int,
        target_latitude_deg: float,
        cloud_cover_fraction: float,
        overlap_fraction: float,
        mode: str = "target_lat_revisit",
    ) -> float:
        payload = {
            "altitude_km": altitude_km,
            "swath_km": swath_km,
            "n_satellites": n_satellites,
            "n_telescopes": n_telescopes,
            "target_latitude_deg": target_latitude_deg,
            "cloud_cover_fraction": cloud_cover_fraction,
            "overlap_fraction": overlap_fraction,
            "mode": mode,
            "orekit": self.jvm_initialized,
        }
        cache_path = self._cache_path(payload)
        if cache_path.exists():
            return float(json.loads(cache_path.read_text(encoding="utf-8"))["revisit_days"])

        period_days = self._period_seconds(altitude_km) / _SIDEREAL_DAY
        latitude_factor = max(0.25, np.cos(np.radians(target_latitude_deg)))
        effective_swath = swath_km * n_telescopes * (1.0 - overlap_fraction)
        if effective_swath <= 0:
            revisit = np.nan
        else:
            target_width = _US_WIDTH_KM if mode == "usa_grid_coverage" else 1_200.0 * latitude_factor
            passes_needed = np.ceil(target_width / effective_swath)
            daylight_penalty = 1.0 / _DAYLIGHT_FACTOR
            cloud_penalty = 1.0 / max(0.05, 1.0 - cloud_cover_fraction)
            revisit = period_days * passes_needed * daylight_penalty * cloud_penalty / n_satellites
            if mode == "usa_grid_coverage":
                revisit *= 0.35
        cache_path.write_text(json.dumps({"revisit_days": revisit}), encoding="utf-8")
        return float(revisit)


def compute_coverage_tensor(config: PipelineConfig, adapter: OrekitCoverageAdapter | None = None) -> LabeledTensor:
    """Compute revisit-day tensor over GSD, constellation, telescope, detector, altitude and swath."""
    adapter = adapter or OrekitCoverageAdapter(
        orekit_data_path=config.coverage.orekit_data_path,
        cache_dir=config.coverage.cache_dir,
    )
    gsd = as_array(config.mission.gsd_m)
    altitudes = as_array(config.mission.altitudes_km)
    swaths = as_array(config.mission.swaths_km)
    shape = (
        len(gsd),
        len(config.constellations),
        len(config.telescopes),
        len(config.detectors),
        len(altitudes),
        len(swaths),
    )
    values = np.full(shape, np.nan, dtype=float)

    for g_i, gsd_m in enumerate(gsd):
        for c_i, constellation in enumerate(config.constellations):
            for t_i, telescope in enumerate(config.telescopes):
                max_swath_fov = 2 * altitudes * np.tan(np.radians(telescope.fov_limit_deg / 2.0))
                for d_i, detector in enumerate(config.detectors):
                    max_detector_swath = detector.max_swath_multiplier * gsd_m * detector.n_pixels / 1000.0
                    for h_i, altitude in enumerate(altitudes):
                        for s_i, swath in enumerate(swaths):
                            if swath > max_swath_fov[h_i] or swath > max_detector_swath:
                                continue
                            revisit = adapter.revisit_days(
                                float(altitude),
                                float(swath),
                                n_satellites=constellation.satellites,
                                n_telescopes=constellation.telescopes_per_satellite,
                                target_latitude_deg=config.coverage.target_latitude_deg,
                                cloud_cover_fraction=config.coverage.cloud_cover_fraction,
                                overlap_fraction=config.coverage.overlap_fraction,
                                mode=config.coverage.mode,
                            )
                            if revisit <= config.thresholds.revisit_days_max:
                                values[g_i, c_i, t_i, d_i, h_i, s_i] = revisit

    return LabeledTensor(
        values=values,
        axes=("gsd", "constellation", "telescope", "detector", "altitude", "swath"),
        coords={
            "gsd": gsd,
            "constellation": np.asarray([c.id for c in config.constellations], dtype=object),
            "telescope": np.asarray([t.id for t in config.telescopes], dtype=object),
            "detector": np.asarray([d.id for d in config.detectors], dtype=object),
            "altitude": altitudes,
            "swath": swaths,
        },
    )
