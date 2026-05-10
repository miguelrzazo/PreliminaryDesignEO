"""Cross-data analysis and optimum configuration selection.

Port of CrossDataFunction.m and OptimumConfigurationAnalysis.m.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Sequence


def cross_data_analysis(
    snr_df: pd.DataFrame,
    mtf_df: pd.DataFrame,
    coverage_df: pd.DataFrame,
    diameters_mm: Sequence[float],
    altitudes_km: Sequence[float],
    fov_limit_deg: float,
    swaths_km: Sequence[float],
    telescope_idx: int = 0,
    snr_req: float = 400.0,
    mtf_req: float = 0.25,
    cov_req: float = 7.0,
    max_refractivo_diam_mm: float = 90.0,
) -> pd.DataFrame:
    """Find minimum pupil diameter satisfying SNR, MTF and coverage at each altitude.

    Parameters
    ----------
    snr_df : pd.DataFrame
        SNR heatmap (altitude × diameter), NaN where below requirement.
    mtf_df : pd.DataFrame
        MTF heatmap (altitude × diameter), NaN where below requirement.
    coverage_df : pd.DataFrame
        Coverage days (altitude × swath), NaN where constraint violated.
    diameters_mm : sequence of float
        Pupil diameters corresponding to DataFrame columns (mm).
    altitudes_km : sequence of float
        Orbital altitudes corresponding to DataFrame rows (km).
    fov_limit_deg : float
        Maximum FOV for the telescope (degrees).
    swaths_km : sequence of float
        Swath values corresponding to coverage_df columns (km).
    telescope_idx : int
        0-based telescope index (0 = Refractivo; enforces diameter ≤ 90 mm).
    snr_req, mtf_req, cov_req : float
        Requirement thresholds (already encoded in NaN masks of input DFs,
        kept here for documentation).
    max_refractivo_diam_mm : float
        Maximum aperture for refractive telescope.

    Returns
    -------
    pd.DataFrame
        Index = altitudes_km; one column 'Dmin_mm' (minimum valid diameter).
        NaN where no valid combination exists.
    """
    h_arr = np.asarray(altitudes_km, dtype=float)
    d_arr = np.asarray(diameters_mm, dtype=float)
    s_arr = np.asarray(swaths_km,    dtype=float)

    snr_mat = snr_df.to_numpy()
    mtf_mat = mtf_df.to_numpy()
    cov_mat = coverage_df.to_numpy()

    dmin = np.full(len(h_arr), np.nan)

    for i, h in enumerate(h_arr):
        if i >= snr_mat.shape[0] or i >= mtf_mat.shape[0] or i >= cov_mat.shape[0]:
            continue

        # Find swaths satisfying coverage requirement at this altitude
        cov_row = cov_mat[i, :]
        valid_swath_idx = np.where(~np.isnan(cov_row) & (cov_row <= cov_req))[0]
        if len(valid_swath_idx) == 0:
            continue

        # Select widest valid swath and check FOV
        best_swath_idx = valid_swath_idx[-1]
        swath = s_arr[best_swath_idx]
        fov_actual = 2 * np.degrees(np.arctan(swath / (2 * h)))
        if fov_actual > fov_limit_deg:
            continue

        # Find minimum diameter satisfying both MTF and SNR
        snr_row = snr_mat[i, :]
        mtf_row = mtf_mat[i, :]
        valid_d = np.where(
            ~np.isnan(snr_row) & (snr_row >= snr_req) &
            ~np.isnan(mtf_row) & (mtf_row >= mtf_req)
        )[0]
        if len(valid_d) == 0:
            continue

        min_d_idx = valid_d[0]
        if min_d_idx >= len(d_arr):
            continue

        d_val = d_arr[min_d_idx]
        if telescope_idx == 0 and d_val > max_refractivo_diam_mm:
            continue

        dmin[i] = d_val

    return pd.DataFrame({"Dmin_mm": dmin}, index=pd.Index(h_arr, name="altitude_km"))


def optimum_configuration(
    mass_results: dict[str, pd.DataFrame],
    configs: Sequence[tuple[int, int]],
) -> dict:
    """Select the global-minimum total-mass configuration.

    Parameters
    ----------
    mass_results : dict
        Keyed by ``"config_name/det_idx/tel_idx"`` (e.g. ``"1sat_1tel/1/Refractivo"``).
        Each value is a DataFrame with column ``Masa_total`` (kg per satellite).
    configs : sequence of (n_sat, n_tel) tuples
        Constellation configurations to evaluate.

    Returns
    -------
    dict with keys: best_key, n_sat, total_mass_kg, per_sat_mass_kg, row.
    """
    best_mass = np.inf
    best_key  = None
    best_row  = None
    best_n_sat = 1

    for n_sat, n_tel in configs:
        config_tag = f"{n_sat}sat_{n_tel}tel"
        for key, df in mass_results.items():
            if not key.startswith(config_tag):
                continue
            valid = df.dropna(subset=["Masa_total"])
            if valid.empty:
                continue
            per_sat = valid["Masa_total"].min()
            total   = per_sat * n_sat
            if total < best_mass:
                best_mass  = total
                best_key   = key
                best_n_sat = n_sat
                best_row   = valid.loc[valid["Masa_total"].idxmin()]

    return {
        "best_key":        best_key,
        "n_sat":           best_n_sat,
        "total_mass_kg":   best_mass,
        "per_sat_mass_kg": best_mass / best_n_sat if best_n_sat else np.nan,
        "row":             best_row,
    }
