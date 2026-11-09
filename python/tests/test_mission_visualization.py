"""Data integrity checks for the offline mission explorer."""
import json
import re

import numpy as np
import pandas as pd
import pytest

from eo_mission.mission_visualization import (
    _earth_position,
    _sample_rows,
    _target_frame,
    _validate_ephemeris,
    build_mission_explorer,
)


@pytest.fixture
def states():
    return pd.DataFrame({
        "time_s": [60.0, 0.0], "satellite": [1, 1],
        "x_m": [7_000_000.0, 7_000_000.0], "y_m": [0.0, 0.0], "z_m": [0.0, 0.0],
        "latitude_deg": [0.0, 0.0], "longitude_deg": [0.0, 0.0],
        "altitude_km": [621.863, 621.863], "sunlit": [True, True],
    })


def test_validation_does_not_modify_analysis_ephemeris(states):
    original = states.copy(deep=True)
    ordered = _validate_ephemeris(states)
    assert ordered["time_s"].tolist() == [0.0, 60.0]
    pd.testing.assert_frame_equal(states, original)


def test_duplicate_states_cannot_produce_ambiguous_animation(states):
    ambiguous = pd.concat([states, states.iloc[:1]], ignore_index=True)
    with pytest.raises(ValueError, match="one state per sample time"):
        _validate_ephemeris(ambiguous)


@pytest.mark.parametrize("column,value", [("x_m", np.nan), ("latitude_deg", 92), ("longitude_deg", -200), ("altitude_km", -1)])
def test_invalid_coordinates_are_rejected(states, column, value):
    states.loc[0, column] = value
    with pytest.raises(ValueError):
        _validate_ephemeris(states)


def test_display_sampling_preserves_first_and_last_states():
    data = pd.DataFrame({"time_s": np.arange(1000)})
    sampled = _sample_rows(data, 90)
    assert len(sampled) == 90
    assert sampled.iloc[0]["time_s"] == 0
    assert sampled.iloc[-1]["time_s"] == 999


def test_wgs84_surface_uses_equatorial_and_polar_radii():
    equator = _earth_position(0, 0)
    pole = _earth_position(90, 0)
    assert equator[0] == pytest.approx(6378.137)
    assert equator[2] == pytest.approx(0)
    assert pole[0] == pytest.approx(0, abs=1e-9)
    assert pole[2] == pytest.approx(6356.752314245, abs=1e-8)


def test_unobserved_targets_are_kept_in_target_grid():
    summary = {"targets": [
        {"target": "observed", "latitude_deg": 25.0, "longitude_deg": -100.0, "access_count": 2},
        {"target": "unobserved", "latitude_deg": 30.0, "longitude_deg": -90.0, "access_count": 0},
    ]}
    targets = _target_frame(summary, pd.DataFrame())
    assert targets["target"].tolist() == ["observed", "unobserved"]
    assert targets["access_count"].tolist() == [2, 0]


def test_html_is_self_contained_and_embeds_original_analysis_counts(states, tmp_path):
    target = tmp_path / "mission.html"
    summary = {"mission": {"altitude_km": 690, "acquired_swath_km": 240},
               "propagation": {"model": "Test data only", "step_s": 60},
               "metrics": {"access_count": 123, "target_count": 10}}
    assert build_mission_explorer(states, pd.DataFrame(), pd.DataFrame(), summary, target) == target
    html = target.read_text()
    assert not re.search(r'<script[^>]+src\s*=', html)
    assert not re.search(r'<link[^>]+href\s*=', html)
    assert "<strong>123</strong>" in html
    assert "Test data only" in html
    embedded = json.loads(re.search(r'<script type="application/json" id="mission-data">(.*?)</script>', html, re.S).group(1))
    assert embedded["summary"] == summary
    assert embedded["satellites"][0]["states"][0][0] == 0.0


def test_html_escapes_metadata_without_corrupting_embedded_data(states, tmp_path):
    label = '</script><script>alert("test")</script>'
    path = build_mission_explorer(states, pd.DataFrame(), pd.DataFrame(), {"mission": {"station_name": label}}, tmp_path / "mission.html")
    html = path.read_text()
    assert label not in html
    embedded = json.loads(re.search(r'<script type="application/json" id="mission-data">(.*?)</script>', html, re.S).group(1))
    assert embedded["summary"]["mission"]["station_name"] == label


def test_reversed_event_boundaries_are_rejected(states, tmp_path):
    events = pd.DataFrame({"satellite": [1], "target": ["station"], "start_s": [60], "end_s": [30], "duration_s": [30]})
    with pytest.raises(ValueError, match="end before it starts"):
        build_mission_explorer(states, pd.DataFrame(), events, {}, tmp_path / "mission.html")
