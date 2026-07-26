"""Cross-data analysis and optimum configuration selection.

Port of CrossDataFunction.m and OptimumConfigurationAnalysis.m.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Sequence

from .config import PipelineConfig
from .tensor import LabeledTensor


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


def rank_feasible_solutions(
    config: PipelineConfig,
    snr: LabeledTensor,
    mtf: LabeledTensor,
    coverage: LabeledTensor,
    mass: LabeledTensor,
) -> pd.DataFrame:
    """Build a ranked table of feasible portfolio solutions.

    A row is feasible when SNR, MTF, coverage, telescope aperture constraints and
    the configured thresholds are all satisfied. Ranking prefers lower total dry
    mass, then lower revisit days, then lower altitude.
    """
    rows: list[dict] = []
    detector_index = {det_id: i for i, det_id in enumerate(snr.coord("detector"))}
    telescope_index = {tel_id: i for i, tel_id in enumerate(snr.coord("telescope"))}
    constellation_index = {c_id: i for i, c_id in enumerate(coverage.coord("constellation"))}

    for g_i, gsd in enumerate(snr.coord("gsd")):
        for b_i, band in enumerate(config.bands):
            for det_id in band.detector_ids:
                d_i = detector_index[det_id]
                detector = config.detector_by_id[det_id]
                for tel_id, t_i in telescope_index.items():
                    telescope = config.telescope_by_id[str(tel_id)]
                    for h_i, altitude in enumerate(snr.coord("altitude")):
                        for dia_i, diameter in enumerate(snr.coord("diameter")):
                            snr_value = snr.values[g_i, b_i, d_i, t_i, h_i, dia_i]
                            mtf_value = mtf.values[g_i, b_i, d_i, t_i, h_i, dia_i]
                            if not (np.isfinite(snr_value) and np.isfinite(mtf_value)):
                                continue
                            if telescope.max_refractive_diameter_mm is not None and diameter > telescope.max_refractive_diameter_mm:
                                continue

                            for const_id, c_i in constellation_index.items():
                                cov_row = coverage.values[g_i, c_i, t_i, d_i, h_i]
                                if not np.isfinite(cov_row).any():
                                    continue
                                best_swath_i = int(np.nanargmin(cov_row))
                                revisit = float(cov_row[best_swath_i])
                                total_mass = float(mass.values[c_i, dia_i])
                                if not np.isfinite(total_mass):
                                    continue
                                rows.append(
                                    {
                                        "gsd_m": float(gsd),
                                        "band": band.id,
                                        "detector": detector.id,
                                        "detector_name": detector.name,
                                        "telescope": telescope.id,
                                        "telescope_name": telescope.name,
                                        "constellation": str(const_id),
                                        "altitude_km": float(altitude),
                                        "diameter_mm": float(diameter),
                                        "swath_km": float(coverage.coord("swath")[best_swath_i]),
                                        "snr": float(snr_value),
                                        "mtf": float(mtf_value),
                                        "revisit_days": revisit,
                                        "total_dry_mass_kg": total_mass,
                                    }
                                )

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df = df.sort_values(
        ["gsd_m", "total_dry_mass_kg", "revisit_days", "altitude_km", "diameter_mm"],
        kind="mergesort",
    ).reset_index(drop=True)
    df["rank_within_gsd"] = df.groupby("gsd_m").cumcount() + 1
    return df


def best_solution_per_gsd(ranked: pd.DataFrame) -> pd.DataFrame:
    """Return the first ranked feasible row for each GSD."""
    if ranked.empty:
        return ranked.copy()
    return ranked.sort_values(["gsd_m", "rank_within_gsd"]).groupby("gsd_m", as_index=False).first()
