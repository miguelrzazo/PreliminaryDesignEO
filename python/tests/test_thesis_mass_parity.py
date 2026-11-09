"""Mass laws checked against MATLAB outputs and physical boundary cases.

The calculations sweep the original thesis altitude/pupil grids first. The
reference report is used only afterwards to compare matching candidate keys;
its dry masses, drag areas and fuel results never enter a model calculation.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import eo_mission.mass as mass_model
from eo_mission.mass import compute_dry_mass, compute_total_mass


ROOT = Path(__file__).parents[2]
REFERENCE = ROOT / "Compiled data/Pipeline/OptimumConfigs/reporte_completo_configuraciones.csv"


@pytest.fixture(scope="module")
def calculated_mass_grid():
    # Original design inputs, independently declared in MATLAB master.m.
    altitudes = np.arange(200, 1001, 10, dtype=float)
    diameters = np.arange(20, 201, 1, dtype=float)
    cases = []
    for telescopes in (1, 2):
        dry = compute_dry_mass(diameters, n_telescopes=telescopes)
        h_grid = np.repeat(altitudes, len(diameters))
        dry_grid = np.tile(dry["dryMass_kg"].to_numpy(), len(altitudes))
        area_grid = np.tile(dry["satelliteArea_m2"].to_numpy(), len(altitudes))
        totals = compute_total_mass(h_grid, dry_grid, area_grid)
        totals["telescopes"] = telescopes
        totals["diameter_mm"] = np.tile(diameters, len(altitudes))
        for column in dry.columns:
            if column not in totals:
                totals[column] = np.tile(dry[column].to_numpy(), len(altitudes))
        cases.append(totals)
    return pd.concat(cases).set_index(["telescopes", "Altura_km", "diameter_mm"])


def test_all_matlab_candidate_mass_budgets_match_fresh_grid(calculated_mass_grid):
    assert REFERENCE.exists(), f"MATLAB reference report missing: {REFERENCE}"
    reference = pd.read_csv(REFERENCE)
    reference = reference.dropna(subset=["Num_Telescopios", "Altura_km", "Diametro_Pupila_mm"])
    assert len(reference) > 3000
    keys = pd.MultiIndex.from_frame(reference[["Num_Telescopios", "Altura_km", "Diametro_Pupila_mm"]])
    calculated = calculated_mass_grid.loc[keys]
    comparisons = {
        "Masa_Seca_Satelite_kg": "dryMass_kg",
        "Masa_Combustible_kg": "fuelMass_kg",
        "Masa_Total_Satelite_kg": "totalMass_kg",
        "Num_Impulsos_Mantenimiento": "impulseCounts",
        "DeltaV_Total_ms": "totalDeltaV_m_s",
        "Masa_Seca_TM_kg": "thematicMapperDryMass_kg",
        "Masa_Seca_SEOSAT_kg": "seosatDryMass_kg",
        "Longitud_Satelite_m": "satelliteLength_m",
        "Area_Satelite_m2": "satelliteArea_m2",
        "Volumen_Satelite_m3": "satelliteVolume_m3",
        "Estandar_U": "cubesatUnits",
    }
    for expected_column, calculated_column in comparisons.items():
        np.testing.assert_allclose(
            calculated[calculated_column].to_numpy(), reference[expected_column].to_numpy(),
            rtol=1e-9, atol=1e-10, equal_nan=True, err_msg=expected_column,
        )
    np.testing.assert_allclose(
        calculated["totalMass_kg"].to_numpy() * reference["Num_Satelites"].to_numpy(),
        reference["Masa_Total_Constelacion_kg"].to_numpy(), rtol=1e-9, atol=1e-10,
        equal_nan=True,
    )


def test_matlab_viable_instrument_fields_match_fresh_grid(calculated_mass_grid):
    reference = pd.read_csv(REFERENCE)
    reference = reference.loc[reference["Es_Viable"] == 1]
    assert len(reference) == 20
    keys = pd.MultiIndex.from_frame(reference[["Num_Telescopios", "Altura_km", "Diametro_Pupila_mm"]])
    calculated = calculated_mass_grid.loc[keys]
    comparisons = {
        "Masa_Instrumento_kg": "instrumentMass_kg",
        "Longitud_Instrumento_m": "instrumentLength_m",
        "Area_Instrumento_m2": "instrumentArea_m2",
        "Volumen_Instrumento_m3": "instrumentVolume_m3",
        "Potencia_Electrica_Instrumento_W": "instrumentPower_W",
    }
    for expected, actual in comparisons.items():
        np.testing.assert_allclose(calculated[actual], reference[expected], rtol=1e-12, atol=1e-13)


@pytest.mark.parametrize("altitude_m,base_density", [
    (0, 1.225), (100_000, 5.297e-7), (200_000, 2.789e-10),
    (350_000, 9.518e-12), (600_000, 1.454e-13),
    (700_000, 3.614e-14), (1000_000, 3.019e-15),
])
def test_density_uses_new_layer_at_breakpoint(altitude_m, base_density):
    assert mass_model._atm_density(altitude_m) == pytest.approx(base_density, rel=1e-14)


def test_density_preserves_layer_jumps_and_high_altitude_extrapolation():
    before = 1.454e-13 * np.exp(-(699.999 - 600) / 71.835)
    assert mass_model._atm_density(699_999) == pytest.approx(before, rel=1e-14)
    assert mass_model._atm_density(1_100_000) == pytest.approx(3.019e-15 * np.exp(-100 / 268), rel=1e-14)
    assert mass_model._atm_density(-100) == 1.225
    np.testing.assert_allclose(mass_model._atm_density(np.array([0, 700_000])), [1.225, 3.614e-14])


@pytest.mark.parametrize("aperture_mm,reference", [(125, mass_model._SEOSAT), (200, mass_model._TM)])
def test_small_aperture_mass_and_power_factor_includes_exact_half(aperture_mm, reference):
    _, _, _, instrument_mass, instrument_power = mass_model._scale_instrument(aperture_mm, reference)
    assert instrument_mass == 2 * reference["weight"] * 0.5**3
    assert instrument_power == 2 * reference["power"] * 0.5**3
    _, _, _, above_mass, above_power = mass_model._scale_instrument(np.nextafter(float(aperture_mm), np.inf), reference)
    assert above_mass == pytest.approx(instrument_mass / 2)
    assert above_power == pytest.approx(instrument_power / 2)


def test_telescope_factor_only_changes_adopted_instrument_mass_and_satellite_size():
    one = compute_dry_mass([38], 1).iloc[0]
    two = compute_dry_mass([38], 2).iloc[0]
    assert two.instrumentMass_kg == pytest.approx(1.5 * one.instrumentMass_kg)
    assert two.instrumentPower_W == one.instrumentPower_W
    assert two.instrumentArea_m2 == one.instrumentArea_m2
    assert two.thematicMapperDryMass_kg == one.thematicMapperDryMass_kg
    assert two.seosatDryMass_kg == one.seosatDryMass_kg
    assert two.satelliteArea_m2 == pytest.approx((two.dryMass_kg / 79) ** (2 / 3))
    assert two.Peso_medio == two.instrumentMass_kg


def test_zero_drag_area_has_no_maintenance_even_at_low_altitude():
    result = compute_total_mass([200, 690], 3.0, 0.0)
    assert np.isinf(result.cycleTimes_s).all()
    assert (result.impulseCounts == 0).all()
    assert (result.totalDeltaV_m_s == 0).all()
    assert (result.fuelMass_kg == 0).all()
    assert (result.totalMass_kg == 3).all()


def test_manoeuvre_exactly_at_mission_end_is_excluded():
    duration = mass_model._MISSION_DURATION
    assert mass_model._completed_cycle_count(duration) == 0
    assert mass_model._completed_cycle_count(duration / 2) == 1
    assert mass_model._completed_cycle_count(duration / (2 + 1e-12)) == 2
    assert mass_model._completed_cycle_count(duration / (2 - 1e-12)) == 1


def test_cycle_quadrature_splits_at_every_interior_atmosphere_boundary(monkeypatch):
    intervals = []
    original_quad = mass_model.quad

    def record_quad(function, lower, upper, **kwargs):
        intervals.append((lower, upper))
        return original_quad(function, lower, upper, **kwargs)

    monkeypatch.setattr(mass_model, "quad", record_quad)
    mass_model._decay_cycle_time.cache_clear()
    result = mass_model._decay_cycle_time(490_000, 710_000, 5, 0.1)
    assert np.isfinite(result)
    assert intervals == [(490_000, 500_000), (500_000, 600_000), (600_000, 700_000), (700_000, 710_000)]


def test_drag_margin_cycle_bounds_and_fuel_margin_are_applied(monkeypatch):
    seen = []

    def cycle_time(lower, upper, dry, area):
        seen.append((lower, upper, dry, area))
        return mass_model._MISSION_DURATION / 2.5

    monkeypatch.setattr(mass_model, "_decay_cycle_time", cycle_time)
    result = compute_total_mass([690], 3.0, 0.12).iloc[0]
    assert seen == [(676_200.0, 703_800.0, 3.0, 0.144)]
    assert result.impulseCounts == 2
    assert result.totalDeltaV_m_s == 2 * result.cycleDeltaV_m_s
    expected_unmargined = 3 * (np.exp(result.totalDeltaV_m_s / (220 * 9.80665)) - 1)
    assert result.unmarginedFuelMass_kg == expected_unmargined
    assert result.fuelMass_kg == expected_unmargined * 1.1
    assert result.totalMass_kg == 3 + result.fuelMass_kg


def test_mass_cutoff_masks_only_mass_outputs(monkeypatch):
    monkeypatch.setattr(mass_model, "_decay_cycle_time", lambda *args: 1.0)
    result = compute_total_mass([200], 3.0, 0.1).iloc[0]
    assert np.isnan(result.totalMass_kg)
    assert np.isnan(result.fuelMass_kg)
    assert np.isnan(result.unmarginedFuelMass_kg)
    assert np.isfinite(result.totalDeltaV_m_s)
    assert result.impulseCounts > 0


def test_mass_cutoff_accepts_exactly_100000_kg():
    result = compute_total_mass([690, 690], [100000, 100001], 0)
    assert result.totalMass_kg.iloc[0] == 100000
    assert np.isnan(result.totalMass_kg.iloc[1])


@pytest.mark.parametrize("dry,area", [([3, 4], [0.1]), ([3], [0.1, 0.2])])
def test_inputs_cannot_be_silently_truncated(dry, area):
    with pytest.raises(ValueError, match="match the altitude count"):
        compute_total_mass([690], dry, area)
