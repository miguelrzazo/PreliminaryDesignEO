"""Dry mass and total mass models.

Port of calcularMasaSeca.m and calcularMasaTotal.m.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from typing import Sequence

# ── Reference instrument data (from calcularMasaSeca.m) ───────────────────────
_TM       = {"aperture": 400.0, "length": 1.0, "weight": 240.0, "power": 280.0}
_SEOSAT   = {"aperture": 250.0, "length": 1.0, "weight": 100.0, "power": 100.0}

# ── Propulsion constants (from calcularMasaTotal.m) ───────────────────────────
_MU_EARTH = 3.986004418e14   # m³/s²
_R_EARTH  = 6_378e3          # m
_CD       = 2.5              # drag coefficient
_ISP      = 220.0            # s (specific impulse, ECAPS HPGP nominal)
_G0       = 9.80665          # m/s²
_MISSION_YEARS = 8
_MISSION_DURATION = _MISSION_YEARS * 365.25 * 24 * 3600  # s


def _scale_instrument(aperture_i: float, ref: dict) -> tuple[float, float, float, float, float]:
    """Scale a reference instrument to aperture_i (mm). Returns L, S, V, W, P."""
    R = aperture_i / ref["aperture"]
    L = R * ref["length"]
    S = L ** 2
    V = L ** 3
    W = R ** 3 * ref["weight"]
    P = R ** 3 * ref["power"]
    return L, S, V, W, P


def _atm_density(h_m: float) -> float:
    """Piecewise exponential atmospheric density (kg/m³) up to 2000 km."""
    h = h_m / 1e3  # km
    table = [
        (100,  6.7e3,  1.225),
        (150,  9.5e3,  4.79e-7),
        (200,  25.5e3, 1.81e-9),
        (250,  37.5e3, 2.53e-10),
        (300,  44.8e3, 6.24e-11),
        (350,  50.3e3, 7.40e-9),
        (400,  54.8e3, 6.98e-12),
        (450,  58.2e3, 2.72e-12),
        (500,  61.3e3, 1.13e-12),
        (600,  70e3,   5e-13),
        (700,  80e3,   1e-13),
        (800,  90e3,   2e-14),
        (900,  100e3,  5e-15),
        (1000, 110e3,  1e-15),
        (1500, 150e3,  1e-16),
        (2000, 200e3,  1e-17),
    ]
    for h_limit, H, rho0 in table:
        if h < h_limit:
            h_base = h_limit - (table[table.index((h_limit, H, rho0))] and
                                next((t[0] for t in reversed(table) if t[0] < h_limit), 0))
            # Simpler: use integer multiples of 50 km as in MATLAB original
            h_floor = (int(h) // 50) * 50
            return rho0 * np.exp(-(h - h_floor) / (H / 1e3))
    return 1e-18


def _decay_ode(t: float, a: float, Cd: float, A: float, mass: float) -> float:
    h = a - _R_EARTH
    rho = _atm_density(h)
    return -Cd * A * rho * np.sqrt(_MU_EARTH * a) / mass


def _hohmann_dv(r1: float, r2: float) -> float:
    a_t = (r1 + r2) / 2
    v1  = np.sqrt(_MU_EARTH / r1)
    vp  = np.sqrt(_MU_EARTH * (2 / r1 - 1 / a_t))
    v2  = np.sqrt(_MU_EARTH / r2)
    va  = np.sqrt(_MU_EARTH * (2 / r2 - 1 / a_t))
    return abs(vp - v1) + abs(v2 - va)


def compute_dry_mass(
    diameters_mm: Sequence[float],
    n_telescopes: int = 1,
) -> pd.DataFrame:
    """Estimate dry mass for each pupil diameter.

    Parameters
    ----------
    diameters_mm : sequence of float
        Pupil diameters (mm).
    n_telescopes : int
        Number of telescopes per satellite (1 or 2).

    Returns
    -------
    pd.DataFrame
        Columns: Diametro_pupila, Longitud_media, Sup_media, Peso_medio,
        Volumen_medio, Potencia_media, Masa_seca (all SI/mm units as in MATLAB).
    """
    d = np.asarray(diameters_mm, dtype=float)
    n = len(d)
    rows = []
    for di in d:
        L_tm, S_tm, V_tm, W_tm, P_tm       = _scale_instrument(di, _TM)
        L_se, S_se, V_se, W_se, P_se       = _scale_instrument(di, _SEOSAT)
        L_m = (L_tm + L_se) / 2
        S_m = (S_tm + S_se) / 2
        V_m = (V_tm + V_se) / 2
        W_m = (W_tm + W_se) / 2
        P_m = (P_tm + P_se) / 2
        masa_seca = 4 * W_m * (1.5 if n_telescopes == 2 else 1.0)
        rows.append((di, L_m, S_m, W_m, V_m, P_m, masa_seca))

    return pd.DataFrame(rows, columns=[
        "Diametro_pupila", "Longitud_media", "Sup_media",
        "Peso_medio", "Volumen_medio", "Potencia_media", "Masa_seca",
    ])


def compute_total_mass(
    altitudes_km: Sequence[float],
    dry_mass_arr: Sequence[float],
    cross_section_arr: Sequence[float],
) -> pd.DataFrame:
    """Compute propellant and total mass for orbital station-keeping over 8 years.

    Uses Hohmann re-boost manoeuvres when altitude decays to 98% of target.

    Parameters
    ----------
    altitudes_km : sequence of float
        Target orbital altitudes (km).
    dry_mass_arr : sequence of float
        Dry mass per altitude (kg).
    cross_section_arr : sequence of float
        Frontal cross-section per altitude (m²); multiplied by 1.1 for margin.

    Returns
    -------
    pd.DataFrame
        Columns: Altura_km, Diametro_pupila (NaN), Masa_seca, Masa_combustible,
        Num_impulsos, DeltaV_total, Masa_total.
    """
    altitudes_km = np.asarray(altitudes_km, dtype=float)
    dry_mass_arr  = np.asarray(dry_mass_arr,  dtype=float)
    cross_section_arr = np.asarray(cross_section_arr, dtype=float)

    results = []
    for h0_km, m_dry, A_ref in zip(altitudes_km, dry_mass_arr, cross_section_arr):
        h_target = h0_km * 1e3
        h_thresh = 0.98 * h_target
        A = A_ref * 1.1

        t_elapsed = 0.0
        h_current = h_target
        total_dv  = 0.0
        impulses  = 0

        while t_elapsed < _MISSION_DURATION:
            a0 = _R_EARTH + h_current
            a_thresh = _R_EARTH + h_thresh
            t_remain = _MISSION_DURATION - t_elapsed

            def event(t, a, a_thresh=a_thresh):
                return a - a_thresh
            event.terminal  = True
            event.direction = -1

            try:
                sol = solve_ivp(
                    lambda t, a: [_decay_ode(t, a[0], _CD, A, m_dry)],
                    [0, t_remain],
                    [a0],
                    events=event,
                    rtol=1e-6, atol=1e-8,
                    max_step=86400,
                    dense_output=False,
                )
                t_elapsed += sol.t[-1]
                h_current  = sol.y[0, -1] - _R_EARTH
            except Exception:
                # Euler fallback
                dt = 3600.0
                a_cur = a0
                while a_cur > a_thresh and t_elapsed < _MISSION_DURATION:
                    da = _decay_ode(0, a_cur, _CD, A, m_dry)
                    a_cur     += da * dt
                    t_elapsed += dt
                h_current = a_cur - _R_EARTH

            if t_elapsed >= _MISSION_DURATION:
                break

            dv = _hohmann_dv(_R_EARTH + h_current, _R_EARTH + h_target)
            total_dv  += dv
            impulses  += 1
            h_current  = h_target

        fuel  = m_dry * (np.exp(total_dv / (_ISP * _G0)) - 1)
        total = m_dry + fuel
        if total > 1e5:
            total = np.nan
            fuel  = np.nan

        results.append({
            "Altura_km":       h0_km,
            "Diametro_pupila": np.nan,
            "Masa_seca":       m_dry,
            "Masa_combustible": fuel,
            "Num_impulsos":    impulses,
            "DeltaV_total":    total_dv,
            "Masa_total":      total,
        })

    return pd.DataFrame(results)
