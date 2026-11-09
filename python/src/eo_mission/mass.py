"""Instrument scaling and eight-year maintenance budgets used by the thesis.

The numerical laws and constants follow MATLAB ``computeDryMass`` and
``computeTotalMass``. Inputs are design quantities; no reference result files
are read by these calculations.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.integrate import quad

from .config import PipelineConfig
from .tensor import LabeledTensor, as_array

_TM = {"aperture": 400.0, "length": 1.0, "weight": 240.0, "power": 280.0}
_SEOSAT = {"aperture": 250.0, "length": 1.0, "weight": 100.0, "power": 100.0}
_MU_EARTH = 3.986004418e14  # m³/s²
_R_EARTH = 6_378e3  # m
_CD = 2.5
_ISP = 220.0  # s
_G0 = 9.80665  # m/s²
_MISSION_YEARS = 8
_MISSION_DURATION = _MISSION_YEARS * 365.25 * 24 * 3600  # s

# Vallado/Curtis piecewise exponential atmosphere, as used by MATLAB.
_ATM_BASE_ALTITUDES_KM = np.array([
    0, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140,
    150, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000,
], dtype=float)
_ATM_BASE_DENSITIES_KG_M3 = np.array([
    1.225, 3.899e-2, 1.774e-2, 3.972e-3, 1.057e-3, 3.206e-4,
    8.770e-5, 1.905e-5, 3.396e-6, 5.297e-7, 9.661e-8, 2.438e-8,
    8.484e-9, 3.845e-9, 2.070e-9, 5.464e-10, 2.789e-10, 7.248e-11,
    2.418e-11, 9.518e-12, 3.725e-12, 1.585e-12, 6.967e-13,
    1.454e-13, 3.614e-14, 1.170e-14, 5.245e-15, 3.019e-15,
])
_ATM_SCALE_HEIGHTS_KM = np.array([
    7.249, 6.349, 6.682, 7.554, 8.382, 7.714, 6.549, 5.799,
    5.382, 5.877, 7.263, 9.473, 12.636, 16.149, 22.523, 29.740,
    37.105, 45.546, 53.628, 53.298, 58.515, 60.828, 63.822,
    71.835, 88.667, 124.64, 181.05, 268.00,
])


def _input_vector(values: Sequence[float] | float, name: str, *, nonnegative: bool = False) -> np.ndarray:
    array = np.atleast_1d(np.asarray(values, dtype=float))
    if array.ndim != 1 or not array.size or not np.isfinite(array).all():
        raise ValueError(f"{name} must be a nonempty finite vector or scalar.")
    invalid = array < 0 if nonnegative else array <= 0
    if invalid.any():
        qualifier = "nonnegative" if nonnegative else "positive"
        raise ValueError(f"{name} must contain {qualifier} values.")
    return array


def _scale_instrument(aperture_i: float, ref: dict) -> tuple[float, float, float, float, float]:
    """Scale reference length, area, volume, mass and power at the pupil size."""
    ratio = aperture_i / ref["aperture"]
    small_aperture_factor = 2.0 if ratio <= 0.5 else 1.0
    length = ratio * ref["length"]
    return (
        length, length**2, length**3,
        ratio**3 * ref["weight"] * small_aperture_factor,
        ratio**3 * ref["power"] * small_aperture_factor,
    )


def _atm_density(h_m: float | np.ndarray) -> float | np.ndarray:
    """Piecewise exponential density in kg/m³; input altitude is in metres.

    Layer breakpoints use the new layer at equality, including its tabulated
    base density. Negative heights clamp to sea level, and the final layer
    extrapolates exponentially above 1000 km, matching the MATLAB model.
    """
    altitude_km = np.maximum(np.asarray(h_m, dtype=float) / 1000, 0)
    layer = np.maximum(np.searchsorted(_ATM_BASE_ALTITUDES_KM, altitude_km, side="right") - 1, 0)
    density = _ATM_BASE_DENSITIES_KG_M3[layer] * np.exp(
        -(altitude_km - _ATM_BASE_ALTITUDES_KM[layer]) / _ATM_SCALE_HEIGHTS_KM[layer]
    )
    return float(density) if np.ndim(density) == 0 else density


def _decay_ode(t: float, a: float, Cd: float, A: float, mass: float) -> float:
    return -Cd * A * _atm_density(a - _R_EARTH) * np.sqrt(_MU_EARTH * a) / mass


def _hohmann_dv(r1: float, r2: float) -> float:
    transfer_semimajor_axis = (r1 + r2) / 2
    lower_speed = np.sqrt(_MU_EARTH / r1)
    transfer_perigee_speed = np.sqrt(_MU_EARTH * (2 / r1 - 1 / transfer_semimajor_axis))
    upper_speed = np.sqrt(_MU_EARTH / r2)
    transfer_apogee_speed = np.sqrt(_MU_EARTH * (2 / r2 - 1 / transfer_semimajor_axis))
    return float((transfer_perigee_speed - lower_speed) + (upper_speed - transfer_apogee_speed))


@lru_cache(maxsize=8192)
def _decay_cycle_time(lower_altitude_m: float, upper_altitude_m: float,
                      dry_mass_kg: float, drag_area_m2: float) -> float:
    if upper_altitude_m <= lower_altitude_m or drag_area_m2 <= 0:
        return float("inf")
    # Mass/area is constant over each model cycle, so its atmospheric integral
    # can be reused across pupil diameters without changing the decay law.
    return dry_mass_kg / (_CD * drag_area_m2) * _atmospheric_decay_integral(
        lower_altitude_m, upper_altitude_m
    )


@lru_cache(maxsize=512)
def _atmospheric_decay_integral(lower_altitude_m: float, upper_altitude_m: float) -> float:
    breakpoints_m = _ATM_BASE_ALTITUDES_KM * 1000
    bounds = np.concatenate((
        [lower_altitude_m],
        breakpoints_m[(breakpoints_m > lower_altitude_m) & (breakpoints_m < upper_altitude_m)],
        [upper_altitude_m],
    ))

    def integrand(altitude_m: float) -> float:
        return 1 / (_atm_density(altitude_m) * np.sqrt(_MU_EARTH * (_R_EARTH + altitude_m)))

    # Separate quadratures preserve the jumps in the tabulated atmosphere.
    return float(sum(
        quad(integrand, lower, upper, epsrel=1e-8, epsabs=1e-3)[0]
        for lower, upper in zip(bounds[:-1], bounds[1:])
    ))


def _completed_cycle_count(cycle_time_s: float) -> int:
    if not np.isfinite(cycle_time_s) or cycle_time_s <= 0:
        return 0
    # A manoeuvre exactly at the eight-year mission endpoint is not counted.
    return int(np.floor((_MISSION_DURATION - np.spacing(_MISSION_DURATION)) / cycle_time_s))


def compute_dry_mass(diameters_mm: Sequence[float], n_telescopes: int = 1) -> pd.DataFrame:
    """Compute per-satellite dry mass and instrument/satellite dimensions.

    Below or at half each reference aperture, reference mass and power include
    a factor of two. Two telescopes multiply the adopted instrument mass by
    1.5; they do not change the reference estimates or adopted dimensions.
    Satellite dimensions use the thesis density of 79 kg/m³.

    Legacy Spanish columns and explicit MATLAB-style SI-unit columns are both
    returned. ``Peso_medio`` is the adopted instrument mass, in kilograms.
    """
    diameters = _input_vector(diameters_mm, "diameters_mm")
    if n_telescopes not in (1, 2):
        raise ValueError("n_telescopes must be 1 or 2.")
    rows = []
    for diameter in diameters:
        tm = _scale_instrument(float(diameter), _TM)
        seosat = _scale_instrument(float(diameter), _SEOSAT)
        instrument_length, instrument_area, instrument_volume, instrument_mass, instrument_power = (
            (left + right) / 2 for left, right in zip(tm, seosat)
        )
        if n_telescopes == 2:
            instrument_mass *= 1.5
        dry_mass = 4 * instrument_mass
        satellite_volume = dry_mass / 79
        satellite_length = satellite_volume ** (1 / 3)
        satellite_area = satellite_length**2
        rows.append({
            "Diametro_pupila": diameter, "Longitud_media": instrument_length,
            "Sup_media": instrument_area, "Peso_medio": instrument_mass,
            "Volumen_medio": instrument_volume, "Potencia_media": instrument_power,
            "Masa_seca": dry_mass, "dryMass_kg": dry_mass,
            "thematicMapperDryMass_kg": 4 * tm[3], "seosatDryMass_kg": 4 * seosat[3],
            "instrumentLength_m": instrument_length, "instrumentArea_m2": instrument_area,
            "instrumentVolume_m3": instrument_volume, "instrumentMass_kg": instrument_mass,
            "instrumentPower_W": instrument_power, "satelliteVolume_m3": satellite_volume,
            "satelliteLength_m": satellite_length, "satelliteArea_m2": satellite_area,
            "cubesatUnits": satellite_volume / 0.1**3,
        })
    return pd.DataFrame(rows)


def dry_mass_array(diameters_mm: Sequence[float], n_telescopes: int = 1) -> np.ndarray:
    return compute_dry_mass(diameters_mm, n_telescopes=n_telescopes)["Masa_seca"].to_numpy()


def compute_mass_tensor(config: PipelineConfig) -> LabeledTensor:
    """Per-constellation dry mass over the configured pupil diameter grid.

    This tensor supplies the auxiliary portfolio's dry-mass proxy. The thesis
    design optimizer must additionally calculate maintenance fuel for every
    candidate altitude through ``compute_total_mass``.
    """
    diameters = as_array(config.mission.diameters_mm)
    values = np.full((len(config.constellations), len(diameters)), np.nan, dtype=float)
    for index, constellation in enumerate(config.constellations):
        values[index] = dry_mass_array(diameters, constellation.telescopes_per_satellite) * constellation.satellites
    return LabeledTensor(
        values=values, axes=("constellation", "diameter"),
        coords={"constellation": np.asarray([case.id for case in config.constellations], dtype=object),
                "diameter": diameters},
    )


def compute_total_mass(
    altitudes_km: Sequence[float],
    dry_mass_arr: Sequence[float] | float,
    cross_section_arr: Sequence[float] | float,
    maintenance_lower_factor: float = 0.98,
    maintenance_upper_factor: float = 1.02,
) -> pd.DataFrame:
    """Eight-year station-keeping budget from the original thesis model.

    Each cycle decays from ``maintenance_upper_factor * h`` to
    ``maintenance_lower_factor * h`` and returns with a complete two-burn
    Hohmann manoeuvre. The atmosphere is integrated separately at every layer
    boundary, drag area receives a 1.2 margin and fuel receives a 1.1 margin.
    Dry mass and reference drag area may be scalar or match the altitude grid.

    ``Num_impulsos``/``impulseCounts`` count complete manoeuvres, not individual
    burns. Results above 100,000 kg mask total and fuel mass with NaN while
    preserving cycle counts and delta-V metadata. No deorbit reserve is added.
    """
    altitudes = _input_vector(altitudes_km, "altitudes_km")
    dry_masses = _input_vector(dry_mass_arr, "dry_mass_arr")
    areas = _input_vector(cross_section_arr, "cross_section_arr", nonnegative=True)
    if dry_masses.size not in (1, altitudes.size) or areas.size not in (1, altitudes.size):
        raise ValueError("Dry mass and drag area must be scalar or match the altitude count.")
    for name, value in (("maintenance_lower_factor", maintenance_lower_factor),
                        ("maintenance_upper_factor", maintenance_upper_factor)):
        if not np.isscalar(value) or not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be a finite positive scalar.")
    dry_masses = np.broadcast_to(dry_masses, altitudes.shape)
    areas = np.broadcast_to(areas, altitudes.shape)
    rows = []
    for altitude_km, dry_mass, area in zip(altitudes, dry_masses, areas):
        nominal_altitude_m = float(altitude_km * 1000)
        lower_m = float(maintenance_lower_factor * nominal_altitude_m)
        upper_m = float(maintenance_upper_factor * nominal_altitude_m)
        cycle_time = _decay_cycle_time(lower_m, upper_m, float(dry_mass), float(area * 1.2))
        cycles = _completed_cycle_count(cycle_time)
        cycle_dv = _hohmann_dv(_R_EARTH + lower_m, _R_EARTH + upper_m)
        total_dv = cycles * cycle_dv
        with np.errstate(over="ignore", invalid="ignore"):
            unmargined_fuel = float(dry_mass * (np.exp(total_dv / (_ISP * _G0)) - 1))
        fuel = unmargined_fuel * 1.1
        total_mass = dry_mass + fuel
        if total_mass > 100000:
            total_mass = fuel = unmargined_fuel = np.nan
        rows.append({
            "Altura_km": altitude_km, "Diametro_pupila": np.nan, "Masa_seca": dry_mass,
            "Masa_combustible": fuel, "Num_impulsos": cycles, "DeltaV_total": total_dv,
            "Masa_total": total_mass, "totalMass_kg": total_mass, "fuelMass_kg": fuel,
            "impulseCounts": cycles, "totalDeltaV_m_s": total_dv,
            "cycleTimes_s": cycle_time, "cycleDeltaV_m_s": cycle_dv,
            "unmarginedFuelMass_kg": unmargined_fuel,
        })
    return pd.DataFrame(rows)
