"""SSO coverage and revisit time calculations.

Port of CoverageSSOAnaliticalfunction.m.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Sequence

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
