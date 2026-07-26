"""Tests for mass, electrical and analysis modules."""
import numpy as np
import pandas as pd
import pytest
from eo_mission.mass import compute_dry_mass, compute_total_mass
from eo_mission.electrical import compute_power_budget
from eo_mission.analysis import cross_data_analysis, optimum_configuration
from eo_mission.optics import compute_snr, compute_mtf


# ── mass.py ──────────────────────────────────────────────────────────────────

def test_dry_mass_scales_with_diameter():
    """Dry mass should increase with aperture diameter."""
    df = compute_dry_mass([50, 100, 200])
    masses = df["Masa_seca"].to_numpy()
    assert np.all(np.diff(masses) > 0), "Dry mass must increase with diameter"


def test_dry_mass_double_telescope_factor():
    """Two telescopes should give 1.5× dry mass vs one."""
    df1 = compute_dry_mass([100], n_telescopes=1)
    df2 = compute_dry_mass([100], n_telescopes=2)
    ratio = df2["Masa_seca"].iloc[0] / df1["Masa_seca"].iloc[0]
    assert abs(ratio - 1.5) < 1e-9


def test_dry_mass_positive():
    """All dry mass values must be strictly positive."""
    df = compute_dry_mass(list(range(20, 201, 10)))
    assert (df["Masa_seca"] > 0).all()


def test_total_mass_greater_than_dry():
    """Total mass must be ≥ dry mass (fuel adds weight)."""
    alts = np.array([400.0, 600.0, 800.0])
    dry  = np.array([5.0, 5.0, 5.0])
    A    = np.array([0.01, 0.01, 0.01])
    df = compute_total_mass(alts, dry, A)
    valid = df.dropna(subset=["Masa_total"])
    assert (valid["Masa_total"] >= valid["Masa_seca"]).all()


# ── coverage.py (analytical tests removed; see test_orekit_integration.py) ──


# ── electrical.py ─────────────────────────────────────────────────────────────

def test_dawn_dusk_no_eclipse():
    """Dawn-dusk orbit (LTAN=6) must have zero eclipse fraction."""
    res = compute_power_budget(altitude_km=760, ltan=6)
    assert res["eclipse_fraction"] == 0.0


def test_noon_midnight_max_eclipse():
    """Noon-midnight orbit (LTAN=12) must have positive eclipse fraction."""
    res = compute_power_budget(altitude_km=760, ltan=12)
    assert res["eclipse_fraction"] > 0.0


def test_power_budget_keys():
    """Return dict must contain all expected keys."""
    res = compute_power_budget()
    expected = {"min_panel_area_m2", "min_battery_cap_wh", "panel_mass_kg",
                "battery_mass_kg", "eclipse_fraction", "orbital_period_s",
                "max_dod_observed"}
    assert expected.issubset(res.keys())


def test_min_panel_area_positive():
    """Minimum panel area must be positive for non-dawn-dusk orbit."""
    res = compute_power_budget(altitude_km=700, ltan=12)
    assert res["min_panel_area_m2"] > 0


# ── analysis.py ───────────────────────────────────────────────────────────────

def _make_optics_dfs():
    alts = np.arange(400, 701, 50, dtype=float)
    diams = np.arange(20, 201, 10, dtype=float)
    snr_df = compute_snr(1.61e-6, 15e-6, 0.6, 0.8, 50.0, 0.0, alts, diams, snr_req=0.0)
    mtf_df = compute_mtf(1.61e-6, 15e-6, 0.45, 0.90, 50.0, 0.0, alts, diams, mtf_req=0.0)
    return snr_df, mtf_df, alts, diams


def test_cross_data_returns_dataframe():
    snr_df, mtf_df, alts, diams = _make_optics_dfs()
    swaths = np.arange(10, 200, 10, dtype=float)
    # Synthetic coverage DataFrame (altitude × swath) with a feasible band of
    # swaths at every altitude; the analytical compute_coverage() was removed
    # so cross_data_analysis is exercised here without an Orekit JVM.
    feasible = (swaths[None, :] >= 60) & (swaths[None, :] <= 150)
    cov = np.where(np.broadcast_to(feasible, (len(alts), len(swaths))), 3.0, np.nan)
    cov_df = pd.DataFrame(cov, index=pd.Index(alts, name="altitude_km"),
                          columns=pd.Index(swaths, name="swath_km"))
    result = cross_data_analysis(snr_df, mtf_df, cov_df, diams, alts, 10.0, swaths)
    assert "Dmin_mm" in result.columns
    assert len(result) == len(alts)


def test_optimum_configuration_finds_minimum():
    """optimum_configuration must pick the key with the lowest total mass."""
    frames = {
        "1sat_1tel/det1/tel1": pd.DataFrame({"Masa_total": [10.0, 20.0]}),
        "1sat_1tel/det1/tel2": pd.DataFrame({"Masa_total": [5.0, 15.0]}),
    }
    result = optimum_configuration(frames, [(1, 1)])
    assert result["per_sat_mass_kg"] == pytest.approx(5.0)
    assert result["best_key"] == "1sat_1tel/det1/tel2"
