"""Independent SciPy J2/J3 reference propagation; this backend does not use Orekit.

Earth rotation uses UTC as UT1 and an approximate sidereal angle, without EOP,
nutation or polar motion. Solar positions use a low-order solar ephemeris.
The thesis mean elements are interpreted as initial osculating elements here.
These choices make this a useful independent geometry check, not an ITRF/SGP4
or mean-element conversion reference.
"""
from __future__ import annotations

from datetime import datetime
import math

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from .mission import EARTH_RADIUS_M, J2, J3, MU, SelectedMission, _inclusive_grid, geodetic_coordinates


def acceleration(position: np.ndarray, *, j2: float = J2, j3: float = J3) -> np.ndarray:
    """Central acceleration and unnormalized zonal J2/J3, in SI units."""
    x, y, z = np.asarray(position, dtype=float)
    r2 = float(x*x + y*y + z*z)
    r = math.sqrt(r2)
    z2_r2 = z*z/r2
    central = -MU * np.array([x, y, z]) / (r2*r)
    zonal2 = 1.5*j2*MU*EARTH_RADIUS_M**2/r**5 * np.array([
        x*(5*z2_r2 - 1), y*(5*z2_r2 - 1), z*(5*z2_r2 - 3)])
    zonal3 = 0.5*j3*MU*EARTH_RADIUS_M**3/r**7 * np.array([
        x*z*(35*z2_r2 - 15), y*z*(35*z2_r2 - 15),
        35*z**4/r2 - 30*z*z + 3*r2])
    return central + zonal2 + zonal3


def initial_state(mission: SelectedMission, raan_rad: float, true_anomaly_rad: float) -> np.ndarray:
    a = EARTH_RADIUS_M + mission.altitude_km*1000
    e = mission.eccentricity
    p = a*(1 - e*e)
    f = true_anomaly_rad
    radius = p/(1 + e*math.cos(f))
    position = radius*np.array([math.cos(f), math.sin(f), 0])
    velocity = math.sqrt(MU/p)*np.array([-math.sin(f), e + math.cos(f), 0])
    omega, inc = math.radians(mission.argument_of_perigee_deg), math.radians(mission.inclination_deg)
    def rz(angle):
        c, s = math.cos(angle), math.sin(angle)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    rx = np.array([[1, 0, 0], [0, math.cos(inc), -math.sin(inc)], [0, math.sin(inc), math.cos(inc)]])
    rotation = rz(raan_rad) @ rx @ rz(omega)
    return np.concatenate([rotation @ position, rotation @ velocity])


def solar_position(julian_dates: np.ndarray) -> np.ndarray:
    """Low-order geocentric Sun vector in approximate equatorial-of-date axes."""
    days = np.asarray(julian_dates) - 2451545.0
    mean_longitude = np.radians((280.460 + 0.9856474*days) % 360)
    mean_anomaly = np.radians((357.528 + 0.9856003*days) % 360)
    longitude = mean_longitude + np.radians(1.915)*np.sin(mean_anomaly) + np.radians(0.020)*np.sin(2*mean_anomaly)
    obliquity = np.radians(23.439 - 0.0000004*days)
    radius = (1.00014 - 0.01671*np.cos(mean_anomaly) - 0.00014*np.cos(2*mean_anomaly))*149_597_870_700
    return radius[:, None]*np.column_stack([np.cos(longitude), np.cos(obliquity)*np.sin(longitude), np.sin(obliquity)*np.sin(longitude)])


def earth_rotation_angle(julian_dates: np.ndarray) -> np.ndarray:
    centuries = (np.asarray(julian_dates) - 2451545.0)/36525
    degrees = 280.46061837 + 360.98564736629*(julian_dates - 2451545.0) + 0.000387933*centuries**2 - centuries**3/38710000
    return np.radians(degrees % 360)


def rotate_to_earth_fixed(xyz: np.ndarray, angle: np.ndarray) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    return np.column_stack([c*xyz[:, 0] + s*xyz[:, 1], -s*xyz[:, 0] + c*xyz[:, 1], xyz[:, 2]])


def propagate_reference(mission: SelectedMission) -> tuple[pd.DataFrame, np.ndarray, dict]:
    """Numerically propagate each spacecraft; never silently replace Orekit."""
    mission.validate()
    seconds = _inclusive_grid(0, mission.duration_days*86400, mission.step_s)
    epoch = datetime.fromisoformat(mission.epoch.replace("Z", "+00:00"))
    jd = epoch.timestamp()/86400 + 2440587.5 + seconds/86400
    solar_inertial = solar_position(jd)
    solar_ra = math.atan2(solar_inertial[0, 1], solar_inertial[0, 0])
    raan = solar_ra + math.radians((mission.ltan_h - 12)*15)
    rotation = earth_rotation_angle(jd)
    sun_fixed = rotate_to_earth_fixed(solar_inertial, rotation)
    frames = []
    def derivative(_time, state):
        return np.concatenate([state[3:], acceleration(state[:3])])
    for satellite in range(mission.satellites):
        state0 = initial_state(mission, raan, 2*math.pi*satellite/mission.satellites)
        solution = solve_ivp(derivative, (0, seconds[-1]), state0, method="DOP853", t_eval=seconds,
                             rtol=1e-10, atol=np.array([1e-3]*3 + [1e-6]*3))
        if not solution.success:
            raise RuntimeError(f"Reference propagation failed: {solution.message}")
        xyz = rotate_to_earth_fixed(solution.y[:3].T, rotation)
        lat, lon, altitude = geodetic_coordinates(xyz)
        sun_direction = sun_fixed/np.linalg.norm(sun_fixed, axis=1)[:, None]
        behind = np.einsum("ij,ij->i", xyz, sun_direction)
        perpendicular = np.linalg.norm(xyz - behind[:, None]*sun_direction, axis=1)
        frames.append(pd.DataFrame({"time_s": seconds, "satellite": satellite + 1,
                                    "x_m": xyz[:, 0], "y_m": xyz[:, 1], "z_m": xyz[:, 2],
                                    "latitude_deg": lat, "longitude_deg": lon,
                                    "altitude_km": altitude/1000,
                                    "sunlit": ~((behind < 0) & (perpendicular < EARTH_RADIUS_M))}))
    return pd.concat(frames, ignore_index=True), sun_fixed, {
        "backend": "scipy-reference", "orekit_executed": False,
        "raan_epoch_deg": math.degrees(raan) % 360,
        "integrator": "SciPy solve_ivp DOP853, rtol=1e-10; atol=1mm/1um/s",
        "gravity": "central attraction + J2 + J3; no drag or higher zonal harmonics",
        "elements": "thesis mean elements interpreted directly as initial osculating elements; no mean-to-osculating conversion",
        "earth_frame": "approximate rotating equatorial frame, UTC=UT1; no EOP/polar motion/nutation; not ITRF",
        "sun": "low-order geocentric solar longitude and obliquity, approximate equatorial-of-date axes",
        "eclipse_display": "cylindrical Earth shadow; target solar elevation controls acquisitions",
        "atmospheric_drag": False,
    }
