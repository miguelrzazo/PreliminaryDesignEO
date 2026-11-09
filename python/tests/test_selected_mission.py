"""Scientific and event-accounting checks independent of an installed JVM."""
from dataclasses import replace
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.integrate import solve_ivp

from eo_mission.config import _expand_range
from eo_mission.mission import (
    EARTH_RADIUS_M, J2, J3, MU, SelectedMission, geodetic_position,
    analyze_accesses, interval_bounds, load_selected_mission, refine_intervals, target_grid, target_metrics,
)
from eo_mission.mission_reference import acceleration, geodetic_coordinates, initial_state
from eo_mission.orekit_runtime import resolve_data_path


def test_selected_design_uses_the_adopted_matlab_csv():
    mission = load_selected_mission(Path(__file__).parents[1] / "config/selected_mission.yaml")
    assert (mission.altitude_km, mission.acquired_swath_km, mission.gsd_m, mission.satellites) == (690, 240, 80, 2)
    assert mission.inclination_deg == pytest.approx(98.1472649550129)
    assert mission.eccentricity == pytest.approx(0.0010448573934173355)
    assert (mission.station_latitude_deg, mission.station_longitude_deg) == (64.84, -147.712)


@pytest.mark.parametrize("latitude,longitude,height", [(0, 0, 0), (64.84, -147.712, 0), (49, -66, 690000), (90, 30, 1000), (-90, -20, 10)])
def test_wgs84_geodetic_roundtrip(latitude, longitude, height):
    point = geodetic_position(latitude, longitude, height)
    lat, lon, h = geodetic_coordinates(point[None, :])
    assert lat[0] == pytest.approx(latitude, abs=1e-8)
    assert lon[0] == pytest.approx(longitude, abs=1e-8)
    assert h[0] == pytest.approx(height, abs=1e-3)


def test_zonal_acceleration_matches_independent_potential_gradient():
    point = np.array([5e6, -3e6, 4e6])
    def potential(xyz):
        r = np.linalg.norm(xyz)
        s = xyz[2]/r
        p2, p3 = (3*s*s - 1)/2, (5*s**3 - 3*s)/2
        return MU/r*(1 - J2*(EARTH_RADIUS_M/r)**2*p2 - J3*(EARTH_RADIUS_M/r)**3*p3)
    gradient = np.array([(potential(point + np.eye(3)[i]) - potential(point - np.eye(3)[i]))/2 for i in range(3)])
    np.testing.assert_allclose(acceleration(point), gradient, rtol=1e-8, atol=1e-8)


def test_retrograde_j2_node_precession_matches_sso_secular_rate():
    mission = SelectedMission()
    state = initial_state(mission, 0, 0)
    samples = np.linspace(0, 2*86400, 1000)
    result = solve_ivp(lambda _, y: np.r_[y[3:], acceleration(y[:3], j3=0)],
                       (0, samples[-1]), state, t_eval=samples, method="DOP853", rtol=1e-10, atol=1e-5)
    assert result.success
    angular_momentum = np.cross(result.y[:3].T, result.y[3:].T)
    node = np.cross(np.array([0, 0, 1]), angular_momentum)
    raan = np.unwrap(np.arctan2(node[:, 1], node[:, 0]))
    measured = np.polyfit(samples, raan, 1)[0]
    a = EARTH_RADIUS_M + mission.altitude_km*1000
    expected = -1.5*math.sqrt(MU/a**3)*J2*(EARTH_RADIUS_M/(a*(1 - mission.eccentricity**2)))**2*math.cos(math.radians(mission.inclination_deg))
    assert measured > 0
    assert measured == pytest.approx(expected, rel=0.025)


def test_sampled_access_crossings_are_interpolated_and_horizon_clipped():
    assert interval_bounds(np.arange(5)*10, np.array([-1, 1, 1, -1, 1])) == [(5.0, 25.0), (35.0, 40.0)]
    assert interval_bounds(np.arange(3)*10, np.ones(3)) == [(0.0, 20.0)]
    assert interval_bounds(np.arange(3)*10, -np.ones(3)) == []


def test_closest_approach_search_recovers_a_pass_between_negative_samples():
    times = np.array([0., 4., 8., 12.])
    margin = lambda t: 0.01 - (t - 5.0)**2/100
    values = np.array([margin(t) for t in times])
    values[1] -= 1e-6
    assert interval_bounds(times, values) == []
    intervals = refine_intervals(times, values, margin, 0.1)
    assert len(intervals) == 1
    assert intervals[0][0] == pytest.approx(4.0, abs=0.001)
    assert intervals[0][1] == pytest.approx(6.0, abs=0.001)


def test_target_metrics_merge_satellites_and_keep_never_observed_targets():
    targets = pd.DataFrame([{"target": "A", "latitude_deg": 40, "longitude_deg": -100}, {"target": "B", "latitude_deg": 42, "longitude_deg": -100}])
    accesses = pd.DataFrame({"target": ["A", "A", "A"], "start_s": [10, 15, 50], "end_s": [20, 25, 60]})
    metrics = target_metrics(targets, accesses, 100).set_index("target")
    assert metrics.loc["A", "access_count"] == 2
    assert metrics.loc["A", "observed_max_gap_days"]*86400 == pytest.approx(25)
    assert metrics.loc["A", "max_gap_days"]*86400 == pytest.approx(40)
    assert metrics.loc["B", "never_observed"]
    assert metrics.loc["B", "max_gap_days"]*86400 == 100
    assert pd.isna(metrics.loc["B", "observed_max_gap_days"])


def test_target_grid_includes_all_study_boundaries():
    grid = target_grid(SelectedMission())
    assert (grid.latitude_deg.min(), grid.latitude_deg.max()) == (25, 49)
    assert (grid.longitude_deg.min(), grid.longitude_deg.max()) == (-125, -66)


def test_nadir_access_requires_illumination_but_radio_contacts_do_not():
    mission = replace(SelectedMission(), duration_days=10/86400, satellites=1,
                      latitude_min_deg=-0.001, latitude_max_deg=0.001,
                      longitude_min_deg=-0.001, longitude_max_deg=0.001,
                      station_latitude_deg=0, station_longitude_deg=0)
    frame = pd.DataFrame({"time_s": [0., 10.], "satellite": [1, 1],
                          "x_m": [EARTH_RADIUS_M + 690000]*2, "y_m": [0., 0.], "z_m": [0., 0.],
                          "latitude_deg": [0., 0.], "longitude_deg": [0., 0.]})
    sun = np.tile(np.array([1.5e11, 0, 0]), (2, 1))
    accesses, contacts, _ = analyze_accesses(frame, sun, mission)
    assert len(accesses) == 4
    assert accesses.duration_s.eq(10).all()
    assert len(contacts) == 1
    night_accesses, night_contacts, metrics = analyze_accesses(frame, -sun, mission)
    assert night_accesses.empty
    assert metrics.never_observed.all()
    pd.testing.assert_frame_equal(contacts, night_contacts)
    over_wide = replace(mission, longitude_min_deg=3, longitude_max_deg=4)
    assert analyze_accesses(frame, sun, over_wide)[0].empty


@pytest.mark.parametrize("changes", [{"step_s": 0}, {"duration_days": -1}, {"overlap_fraction": 1}, {"target_latitude_step_deg": 0}, {"satellites": 0}, {"epoch": "2024-04-28T00:00:00"}])
def test_invalid_mission_inputs_fail_before_propagation(changes):
    with pytest.raises(ValueError):
        replace(SelectedMission(), **changes).validate()


def test_data_path_honors_the_environment_and_rejects_missing_data(tmp_path, monkeypatch):
    data = tmp_path / "orekit-data.zip"
    data.write_bytes(b"placeholder")
    monkeypatch.setenv("OREKIT_DATA_PATH", str(data))
    assert resolve_data_path() == data
    with pytest.raises(FileNotFoundError):
        resolve_data_path(tmp_path / "missing.zip")


def test_bad_grid_range_fails_instead_of_looping_forever():
    with pytest.raises(ValueError):
        _expand_range({"start": 400, "stop": 900, "step": 0})


def test_runtime_rejects_switching_data_behind_cached_java_time_scales(tmp_path, monkeypatch):
    from eo_mission import orekit_runtime
    old, new = tmp_path / "first.zip", tmp_path / "second.zip"
    old.write_bytes(b"first")
    new.write_bytes(b"second")
    monkeypatch.setattr(orekit_runtime, "_registered_paths", {old})
    with pytest.raises(RuntimeError, match="Restart Python"):
        orekit_runtime.initialize_runtime(new)
