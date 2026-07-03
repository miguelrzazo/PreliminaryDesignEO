"""Electrical power budget simulation.

Port of ElectricalSimulation.m.
"""
from __future__ import annotations

import numpy as np

_R_EARTH   = 6_371.0   # km
_MU_EARTH  = 398_600.0 # km³/s²
_IRRADIANCE = 1_366.0  # W/m²


def compute_power_budget(
    altitude_km: float = 760.0,
    ltan: int = 6,
    panel_area: float = 0.05,
    panel_efficiency: float = 0.30,
    panel_degradation_rate: float = 0.045,
    panel_density: float = 15.0,
    power_consumption: float = 20.0,
    battery_capacity_wh: float = 10.0,
    max_dod: float = 0.90,
    battery_energy_density: float = 200.0,
    charge_efficiency: float = 0.90,
    discharge_efficiency: float = 0.95,
    battery_degradation_rate: float = 0.02,
    mission_years: int = 8,
    pts_per_day: int = 100,
) -> dict:
    """Run an 8-year energy budget simulation for a LEO satellite.

    Parameters
    ----------
    altitude_km : float
        Orbital altitude (km).
    ltan : int
        Local Time of Ascending Node (hours). 6/18 → dawn-dusk (no eclipse);
        0/12 → noon-midnight (maximum eclipse).
    panel_area : float
        Solar panel area (m²).
    panel_efficiency : float
        BOL panel efficiency (fraction).
    panel_degradation_rate : float
        Annual efficiency loss (fraction/year).
    panel_density : float
        Panel areal mass density (kg/m²).
    power_consumption : float
        Mean power draw (W).
    battery_capacity_wh : float
        BOL battery capacity (Wh).
    max_dod : float
        Maximum depth of discharge (fraction).
    battery_energy_density : float
        Battery energy density (Wh/kg).
    charge_efficiency, discharge_efficiency : float
        Round-trip efficiency fractions.
    battery_degradation_rate : float
        Annual capacity loss (fraction/year).
    mission_years : int
        Simulation duration (years).
    pts_per_day : int
        Time resolution (points per day).

    Returns
    -------
    dict with keys:
        min_panel_area_m2     — minimum panel area for EOL power balance (m²)
        min_battery_cap_wh    — minimum battery capacity for eclipse (Wh)
        panel_mass_kg         — mass for min_panel_area
        battery_mass_kg       — mass for min_battery_cap
        eclipse_fraction      — orbit fraction in eclipse
        orbital_period_s      — orbital period (s)
        max_dod_observed      — maximum DoD observed in simulation
    """
    r_orb = _R_EARTH + altitude_km               # km
    T_orb = 2 * np.pi * np.sqrt(r_orb**3 / _MU_EARTH)  # s

    # Eclipse fraction by LTAN
    if ltan in (6, 18):
        eclipse_frac = 0.0
    else:
        beta = np.arcsin(_R_EARTH / r_orb)
        eclipse_frac = 2 * beta / (2 * np.pi)

    t_eclipse = T_orb * eclipse_frac  # s per orbit

    # ── Vectorised 8-year simulation ──────────────────────────────────────────
    total_days  = mission_years * 365
    n_pts       = total_days * pts_per_day
    t_days      = np.linspace(0, total_days, n_pts)
    dt_h        = (t_days[1] - t_days[0]) * 24

    T_orb_days  = T_orb / 86400
    frac_orb    = np.mod(t_days / T_orb_days, 1)
    angle_rad   = frac_orb * 2 * np.pi
    in_eclipse  = np.abs(np.sin(angle_rad)) < eclipse_frac / 2

    degradation = 1 - panel_degradation_rate * t_days / 365
    eta_current = panel_efficiency * degradation
    day_of_year = np.mod(t_days, 365)
    irr_var     = 1 + 0.034 * np.cos(2 * np.pi * (day_of_year - 4) / 365)

    pwr_gen           = panel_area * eta_current * _IRRADIANCE * irr_var
    pwr_gen[in_eclipse] = 0.0

    bat_energy  = np.zeros(n_pts)
    bat_energy[0] = battery_capacity_wh * 0.8
    bat_deg     = 1 - battery_degradation_rate * t_days / 365
    bat_cap     = battery_capacity_wh * bat_deg

    for k in range(1, n_pts):
        consumed = power_consumption * dt_h
        if not in_eclipse[k]:
            net = pwr_gen[k] * dt_h - consumed
            bat_energy[k] = bat_energy[k-1] + net * charge_efficiency
        else:
            bat_energy[k] = bat_energy[k-1] - consumed / discharge_efficiency
        bat_energy[k] = np.clip(bat_energy[k], 0, bat_cap[k])

    # ── Theoretical sizing ────────────────────────────────────────────────────
    sun_frac = 1 - eclipse_frac
    eta_eol  = panel_efficiency * (1 - panel_degradation_rate * mission_years)
    min_panel_area = (power_consumption / (sun_frac * eta_eol * _IRRADIANCE)
                      if sun_frac > 0 else np.inf)

    min_bat_cap = 0.0
    bat_mass    = 0.0
    if t_eclipse > 0:
        e_eclipse   = power_consumption * t_eclipse / 3600
        min_bat_cap = e_eclipse / (max_dod * discharge_efficiency)
        bat_mass    = min_bat_cap / battery_energy_density

    e_min = float(np.min(bat_energy))
    dod_max = (battery_capacity_wh - e_min) / battery_capacity_wh

    return {
        "min_panel_area_m2":   min_panel_area,
        "min_battery_cap_wh":  min_bat_cap,
        "panel_mass_kg":       min_panel_area * panel_density,
        "battery_mass_kg":     bat_mass,
        "eclipse_fraction":    eclipse_frac,
        "orbital_period_s":    T_orb,
        "max_dod_observed":    dod_max,
    }
