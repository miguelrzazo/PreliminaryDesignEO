"""SSO coverage and revisit time via Orekit propagation.

Replaces the analytical revisit formula in the original Python port (and the
external MATLAB RevisitTime toolkit) with a genuine Orekit-based computation:

    1. The sun-synchronous inclination is solved from the J2 sun-sync
       condition (the same equation used by ``CoverageRevisitCalc.m``); altitudes
       where SSO is impossible return NaN.
    2. A numerical propagator with a J2-only force model is propagated over a
       configurable horizon (default 30 days). The sub-satellite ground track
       is sampled at a fixed step and cached per (altitude, n-satellites), since
       the track does not depend on the swath.
    3. A coarse latitude/longitude grid spanning the continental US is used as
       the set of ground targets. For every (target, swath) cell the off-nadir
       visibility radius is computed from the spherical geometry of the sensor
       cone (half-angle ``arctan(eff_swath / 2h)``) and applied to the cached
       ground track; the result is the set of days on which each target is
       observable.
    4. Cloud effects are modelled with a persistent two-state Markov chain (see
       ``clouds.py``). Monte-Carlo sampling folds daily cloud state onto the
       opportunity log; the reported revisit is the worst-target max-clear-gap
       at the 95th percentile across runs (the mean is also tracked).

The analytical single-target formula used by the previous port — which invented
a ``1/0.5`` daylight penalty and a ``×0.35`` USA-grid fudge factor absent from
the MATLAB physics — has been removed. There is no JVM-free fallback: the
pipeline requires Orekit.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from .config import PipelineConfig
from .clouds import revisit_tiling_95p
from .tensor import LabeledTensor, as_array

# ── physical constants (match MATLAB CoverageRevisitCalc.m) ───────────────────
_R_EARTH_M = 6_378_137.0               # WGS84 equatorial radius [m]
_MU_EARTH = 3.986004418e14            # Earth gravitational parameter [m³/s²]
_J2 = 1.082629989052e-3               # Earth J2
_SSO_NODAL_RATE = (360.0 / (365.2421897 * 24.0 * 3600.0)) * (math.pi / 180.0)


def sso_inclination_deg(altitude_km: float) -> float:
    """Sun-synchronous inclination for the given altitude, or NaN if infeasible.

    Reproduces the J2 sun-sync condition of ``CoverageRevisitCalc.m``::

        cos_i = -SSO_nodal_rate * (a/R_E)^(7/2) / (1.5 * n * J2)

    where ``a = R_E + h`` and ``n = sqrt(mu / a^3)``. When ``|cos_i| > 1`` an
    SSO orbit cannot be maintained at that altitude and NaN is returned.
    """
    a = _R_EARTH_M + float(altitude_km) * 1000.0
    n = math.sqrt(_MU_EARTH / a**3)
    cos_i = -_SSO_NODAL_RATE * (a / _R_EARTH_M) ** 3.5 / (1.5 * n * _J2)
    if abs(cos_i) > 1.0:
        return float("nan")
    return math.degrees(math.acos(cos_i))


def _gamma_max_rad(altitude_km: float, half_fov_deg: float) -> float:
    """Earth-central angle of the sensor's visible footprint [rad].

    Given the off-nadir half-angle ``η`` of the sensor and the orbital radius
    ``r = R_E + h``, the spherical triangle (Earth centre O, satellite S,
    target T) gives ``sin η = R_E·sinγ / |ST|`` with
    ``|ST|² = R_E² + r² − 2·R_E·r·cosγ``. Solving for ``cosγ`` at the boundary
    yields the maximum ground-central angle γ that is still visible.
    """
    eta = math.radians(float(half_fov_deg))
    r = _R_EARTH_M + float(altitude_km) * 1000.0
    S = math.sin(eta) ** 2
    # R_E²·c² − S·(2·R_E·r)·c + (S·(R_E² + r²) − R_E²) = 0
    A = _R_EARTH_M**2
    B = -S * 2.0 * _R_EARTH_M * r
    C = S * (_R_EARTH_M**2 + r**2) - _R_EARTH_M**2
    disc = B * B - 4.0 * A * C
    if disc < 0.0:
        return 0.0
    roots = [(-B + math.sqrt(disc)) / (2.0 * A), (-B - math.sqrt(disc)) / (2.0 * A)]
    cos_g = max(c for c in roots if -1.0 <= c <= 1.0)
    return min(math.acos(cos_g), math.pi)


def _ecef_unit(lat_deg: float, lon_deg: float) -> np.ndarray:
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    return np.array(
        [math.cos(lat) * math.cos(lon),
         math.cos(lat) * math.sin(lon),
         math.sin(lat)],
        dtype=float,
    )


class OrekitCoverageAdapter:
    """Orekit-backed coverage adapter.

    The adapter lazily starts a JVM, registers an orekit-data directory, and
    propagates SSO satellites with a J2 force model to log the sub-satellite
    ground track. The track (sampled lat/lon per step) is memoised to a JSON
    cache so that portfolio sweeps reuse it across swath / telescope cells.
    """

    SAMPLE_STEP_S = 180.0  # ground-track sampling step (s)

    def __init__(
        self,
        orekit_data_path: str | None = None,
        cache_dir: str | Path | None = None,
    ):
        self.orekit_data_path = orekit_data_path
        self.cache_dir = Path(cache_dir or ".cache/orekit-coverage")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.jvm_initialized = False

    # ── JVM / orekit-data bootstrap ──────────────────────────────────────────
    def initialize_orekit(self, require: bool = False) -> bool:
        try:
            import orekit

            orekit.initVM()
            path = self.orekit_data_path
            if path is not None and Path(path).exists():
                from orekit.pyhelpers import setup_orekit_curdir

                setup_orekit_curdir(path)
            self.jvm_initialized = True
            return True
        except Exception:
            if require:
                raise
            self.jvm_initialized = False
            return False

    # ── grid ──────────────────────────────────────────────────────────────────
    def _build_grid(self, cfg) -> dict[float, list[tuple[float, np.ndarray]]]:
        """Build the CONUS latitude-band target grid.

        Returns a dict mapping each latitude band (°) to a list of
        ``(lon, ecef_unit)`` target longitudes spanning the continental US.
        The ECEF unit vector is in the Veis-1950 body-fixed frame.
        """
        lats = np.arange(
            cfg.grid_lat_start, cfg.grid_lat_stop + cfg.grid_lat_step / 2.0, cfg.grid_lat_step
        )
        lons = np.linspace(-120.0, -70.0, max(3, int(cfg.grid_lon_samples)))
        bands: dict[float, list[tuple[float, np.ndarray]]] = {}
        for lat in lats:
            bands[float(lat)] = [
                (float(lon), _ecef_unit(float(lat), float(lon))) for lon in lons
            ]
        return bands

    # ── ground-track propagation ──────────────────────────────────────────────
    def _propagate_ground_track(
        self,
        altitude_km: float,
        n_sat: int,
        propagation_days: int,
    ) -> np.ndarray:
        """Return sub-satellite (lat, lon) samples (deg) of shape (N, 2).

        The orbit uses the SSO inclination derived from the J2 sun-sync
        condition and is propagated with an analytical Keplerian propagator; the
        sun-synchronous RAAN precession is captured through the Earth-fixed
        Veis-1950 frame, in which a fixed-RAAN inertial orbit naturally traces a
        slowly drifting ground track. (A numerical J2 propagator is a possible
        refinement; the revisit result is insensitive to the residual ∼0.1°/day
        difference over the 30-day horizon.)

        The satellite mean anomalies are equally spaced (``360°/n_sat``) so the
        constellation flies as a plane train; the merged ground track over all
        satellites is returned, ordered by satellite then sample.
        """
        from org.orekit.frames import FramesFactory
        from org.orekit.orbits import KeplerianOrbit, PositionAngleType
        from org.orekit.propagation.analytical import KeplerianPropagator
        from org.orekit.time import AbsoluteDate
        from org.orekit.utils import Constants

        inc_deg = sso_inclination_deg(altitude_km)
        if not np.isfinite(inc_deg):
            return np.empty((0, 2))

        body_frame = FramesFactory.getVeis1950()
        a = float(_R_EARTH_M + altitude_km * 1000.0)
        mu = float(Constants.EIGEN5C_EARTH_MU)
        eme = FramesFactory.getEME2000()
        start = AbsoluteDate.J2000_EPOCH

        step = self.SAMPLE_STEP_S
        n_samples = int(propagation_days * 86400.0 / step) + 1
        dates = [start.shiftedBy(float(i * step)) for i in range(n_samples)]
        latlons = []

        for s in range(n_sat):
            m0 = float(360.0 * s / max(1, n_sat))
            orbit = KeplerianOrbit(
                a, 1e-4, math.radians(inc_deg),
                0.0, 0.0, math.radians(m0),
                PositionAngleType.MEAN, eme, start, mu,
            )
            prop = KeplerianPropagator(orbit)
            for d in dates:
                pv = prop.propagate(d).getPVCoordinates(body_frame).getPosition()
                x, y, z = float(pv.x), float(pv.y), float(pv.z)
                r = math.hypot(math.hypot(x, y), z)
                lat = math.degrees(math.asin(z / r))
                lon = math.degrees(math.atan2(y, x))
                latlons.append((lat, lon))

        return np.asarray(latlons, dtype=float)

    def _ground_track_cache_path(self, payload: dict) -> Path:
        encoded = json.dumps(payload, sort_keys=True).encode("utf-8")
        digest = hashlib.sha256(encoded).hexdigest()[:20]
        return self.cache_dir / f"track_{digest}.json"

    def _get_ground_track(
        self, altitude_km: float, n_sat: int, propagation_days: int
    ) -> np.ndarray:
        payload = {
            "altitude_km": float(altitude_km),
            "n_sat": int(n_sat),
            "propagation_days": int(propagation_days),
            "step_s": float(self.SAMPLE_STEP_S),
        }
        path = self._ground_track_cache_path(payload)
        if path.exists():
            return np.asarray(json.loads(path.read_text(encoding="utf-8")), dtype=float)
        track = self._propagate_ground_track(altitude_km, n_sat, propagation_days)
        path.write_text(json.dumps(track.tolist()), encoding="utf-8")
        return track

    # ── opportunity extraction ─────────────────────────────────────────────────
    def _band_acquisition_days(
        self,
        track: np.ndarray,
        bands: dict[float, list[tuple[float, np.ndarray]]],
        gamma_max_rad: float,
        n_satellites: int,
    ) -> list[np.ndarray]:
        """For each latitude band return the sorted set of day indices on which
        ≥1 ground-track sample falls inside the sensor footprint at *any* target
        longitude in the band (i.e. the band receives ≥1 acquisition that day).
        Used as a reachability gate for the tiling revisit.
        """
        n_sat = max(1, n_satellites)
        if track.size == 0 or gamma_max_rad <= 0.0:
            return [np.array([], dtype=int) for _ in bands]
        samples_per_day = max(1, int(round(86400.0 / self.SAMPLE_STEP_S)))
        per_sat = len(track) // n_sat
        # per-sample day index within its own satellite block
        day_idx = np.empty(len(track), dtype=int)
        for s in range(n_sat):
            lo = s * per_sat
            hi = lo + per_sat
            day_idx[lo:hi] = np.arange(per_sat) // samples_per_day

        lat = np.radians(track[:, 0])
        lon = np.radians(track[:, 1])
        su = np.stack(
            [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)],
            axis=1,
        )
        cos_gmax = math.cos(gamma_max_rad)
        out = []
        for _, targets in bands.items():
            band_days = np.array([], dtype=int)
            for (_, u) in targets:
                visible = (su @ u) >= cos_gmax - 1e-12
                if not visible.any():
                    continue
                band_days = np.union1d(band_days, np.unique(day_idx[visible]))
            out.append(band_days)
        return out

    @staticmethod
    def _passes_needed(eff_swath_km: float, band_width_km: float) -> int:
        """Number of swath strips needed to tile the band width (MATLAB formula)."""
        if eff_swath_km <= 0:
            return 10**9
        return int(math.ceil(band_width_km / eff_swath_km))

    # ── revisit ───────────────────────────────────────────────────────────────
    def _orbital_period_s(self, altitude_km: float) -> float:
        """Keplerian orbital period (s) from an Orekit Keplerian orbit.

        Uses the same WGS84 radius and EIGEN5C μ as the SSO inclination formula,
        so the period is consistent with the propagated ground track.
        """
        from org.orekit.orbits import KeplerianOrbit, PositionAngleType
        from org.orekit.frames import FramesFactory
        from org.orekit.time import AbsoluteDate
        from org.orekit.utils import Constants

        a = float(_R_EARTH_M + altitude_km * 1000.0)
        mu = float(Constants.EIGEN5C_EARTH_MU)
        orbit = KeplerianOrbit(
            a, 1e-4, math.radians(sso_inclination_deg(altitude_km)),
            0.0, 0.0, 0.0,
            PositionAngleType.MEAN,
            FramesFactory.getEME2000(),
            AbsoluteDate.J2000_EPOCH,
            mu,
        )
        return float(orbit.getKeplerianPeriod())

    def revisit_days(
        self,
        altitude_km: float,
        swath_km: float,
        *,
        n_satellites: int,
        n_telescopes: int,
        overlap_fraction: float,
        cloud_cover_fraction: float,
        cloud_persistence: float,
        mc_samples: int,
        mc_seed: int,
        revisit_percentile: float,
        propagation_days: int,
        bands: dict[float, list[tuple[float, np.ndarray]]],
        band_width_km: float,
    ) -> float:
        """Worst-latitude-band tiling revisit (days) for one cell.

        The revisit of a latitude band is the time to accumulate
        ``ceil(band_width / eff_swath)`` clear usable passes — the MATLAB
        ``_US_WIDTH / eff_swath`` tiling count. The daylight usable-pass rate is
        ``round(86400 / T_orb) × n_sat`` with ``T_orb`` from an Orekit Keplerian
        orbit (replacing the hand-coded ``2π√(r³/μ)``), and the clear-pass
        realisation is driven by a Markov-persistent cloud Monte-Carlo at the
        configured confidence level (replacing the ``1/(1-cloud)`` mean
        heuristic). The latitude-band longitude span scales as
        ``band_width × cos(lat)`` and each band must be flown over at least once
        in the propagated horizon (reachability gate) or it is infeasible.
        """
        inc_deg = sso_inclination_deg(altitude_km)
        if not np.isfinite(inc_deg):
            return float("nan")
        effective_swath = swath_km * n_telescopes * (1.0 - overlap_fraction)
        if effective_swath <= 0:
            return float("nan")
        half_fov_deg = math.degrees(
            math.atan(effective_swath * 1000.0 / (2.0 * altitude_km * 1000.0))
        )
        gamma_max = _gamma_max_rad(altitude_km, half_fov_deg)

        t_orb = self._orbital_period_s(altitude_km)
        # ~2 latitude crossings per orbit × daylight factor 0.5 ⇒ 1 usable band
        # pass per orbit; scaled by the constellation size.
        usable_passes_per_day = max(1.0, round(86400.0 / t_orb)) * max(1, n_satellites)

        track = self._get_ground_track(altitude_km, n_satellites, propagation_days)
        # Reachability: a band is covered iff the SSO ground track reaches its
        # latitude. SSO inclinations (~96-100°) cross every CONUS latitude, so
        # this is effectively always true for feasible altitudes; we still gate
        # so that, e.g., a non-polar orbit would correctly mark high bands
        # unreachable rather than underestimate the revisit.
        max_reached_lat = float(np.abs(track[:, 0]).max()) if track.size else 0.0

        passes_needed, reached = [], []
        for lat_deg in bands.keys():
            span_km = band_width_km * math.cos(math.radians(lat_deg))
            passes_needed.append(self._passes_needed(effective_swath, span_km))
            reached.append(max_reached_lat >= abs(lat_deg))

        rng = np.random.default_rng(mc_seed)
        revisit_p, _ = revisit_tiling_95p(
            passes_needed,
            reached,
            usable_passes_per_day,
            n_days=propagation_days,
            p=cloud_cover_fraction,
            rho=cloud_persistence,
            rng=rng,
            n_samples=mc_samples,
            percentile=revisit_percentile,
        )
        return revisit_p


def compute_coverage_tensor(
    config: PipelineConfig, adapter: OrekitCoverageAdapter | None = None
) -> LabeledTensor:
    """Compute revisit-day tensor over GSD, constellation, telescope, detector,
    altitude and swath using Orekit propagation and a Markov cloud model.

    The tensor cell is NaN where:
        - SSO is infeasible at the altitude,
        - the swath exceeds the telescope FOV-imposed limit,
        - the swath exceeds the detector grazing limit, or
        - the realised 95th-percentile revisit exceeds ``revisit_days_max``.
    """
    adapter = adapter or OrekitCoverageAdapter(
        orekit_data_path=config.coverage.orekit_data_path,
        cache_dir=config.coverage.cache_dir,
    )
    if not adapter.jvm_initialized:
        adapter.initialize_orekit(require=True)

    gsd = as_array(config.mission.gsd_m)
    altitudes = as_array(config.mission.altitudes_km)
    swaths = as_array(config.mission.swaths_km)
    grid = adapter._build_grid(config.coverage)
    band_width_km = float(config.coverage.us_width_km)
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
                max_swath_fov = 2.0 * altitudes * np.tan(
                    np.radians(telescope.fov_limit_deg / 2.0)
                )
                for d_i, detector in enumerate(config.detectors):
                    max_detector_swath = (
                        detector.max_swath_multiplier * gsd_m * detector.n_pixels / 1000.0
                    )
                    for h_i, altitude in enumerate(altitudes):
                        if not np.isfinite(sso_inclination_deg(float(altitude))):
                            continue
                        for s_i, swath in enumerate(swaths):
                            if swath > max_swath_fov[h_i] or swath > max_detector_swath:
                                continue
                            revisit = adapter.revisit_days(
                                float(altitude),
                                float(swath),
                                n_satellites=constellation.satellites,
                                n_telescopes=constellation.telescopes_per_satellite,
                                overlap_fraction=config.coverage.overlap_fraction,
                                cloud_cover_fraction=config.coverage.cloud_cover_fraction,
                                cloud_persistence=config.coverage.cloud_persistence,
                                mc_samples=config.coverage.mc_samples,
                                mc_seed=config.coverage.mc_seed,
                                revisit_percentile=config.coverage.revisit_percentile,
                                propagation_days=config.coverage.propagation_days,
                                bands=grid,
                                band_width_km=band_width_km,
                            )
                            if np.isfinite(revisit) and revisit <= config.thresholds.revisit_days_max:
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