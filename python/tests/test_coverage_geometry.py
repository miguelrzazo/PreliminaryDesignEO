"""Pure-Python checks for the SSO geometry and ground-track cache."""
import json
import math

import numpy as np
import pytest

from eo_mission.coverage import OrekitCoverageAdapter, sso_inclination_deg


@pytest.mark.parametrize("altitude_km", [400.0, 630.0, 690.0, 1000.0])
def test_inclination_produces_one_solar_year_nodal_period(altitude_km):
    inclination = math.radians(sso_inclination_deg(altitude_km))
    radius = 6_378_137.0
    axis = radius + altitude_km * 1000
    mean_motion = math.sqrt(3.986004418e14 / axis**3)
    nodal_rate = -1.5 * mean_motion * 1.082629989052e-3 * (radius / axis)**2 * math.cos(inclination)
    assert math.isclose(2 * math.pi / nodal_rate / 86400, 365.2421897, rel_tol=1e-12)
    if altitude_km == 690:
        assert math.isclose(math.degrees(inclination), 98.147264955012915, abs_tol=1e-10)


def test_ground_track_does_not_reuse_old_inclination_cache(tmp_path, monkeypatch):
    adapter = OrekitCoverageAdapter(cache_dir=tmp_path)
    old_payload = {"altitude_km": 690.0, "n_sat": 2, "propagation_days": 1,
                   "step_s": float(adapter.SAMPLE_STEP_S)}
    old_path = adapter._ground_track_cache_path(old_payload)
    old_path.write_text(json.dumps([[99.0, -99.0]]))
    expected = np.array([[25.0, -100.0]])
    calls = []

    def propagate(*args):
        calls.append(args)
        return expected

    monkeypatch.setattr(adapter, "_propagate_ground_track", propagate)
    np.testing.assert_array_equal(adapter._get_ground_track(690, 2, 1), expected)
    np.testing.assert_array_equal(adapter._get_ground_track(690, 2, 1), expected)
    assert len(calls) == 1
    assert old_path.exists()
