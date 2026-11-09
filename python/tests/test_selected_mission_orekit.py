"""Exercise the selected-mission Java APIs with real physical data."""
from dataclasses import replace
import os
from pathlib import Path

import numpy as np
import pytest

from eo_mission.mission import SelectedMission, analyze_accesses, propagate_selected_mission


def test_selected_mission_propagates_a_real_zonal_orbit_and_computes_events():
    pytest.importorskip("orekit")
    data = os.environ.get("OREKIT_DATA_PATH")
    if not data or not Path(data).exists():
        pytest.skip("Selected-mission integration needs OREKIT_DATA_PATH.")
    mission = replace(SelectedMission(), duration_days=0.08, satellites=1, step_s=10)
    ephemeris, sun, provenance = propagate_selected_mission(mission, data)
    assert len(ephemeris) > 600
    assert ephemeris.altitude_km.between(650, 750).all()
    assert ephemeris.latitude_deg.max() > 80
    assert ephemeris.latitude_deg.min() < -80
    assert np.linalg.norm(sun, axis=1).min() > 1.4e11
    assert "ITRF" in provenance["earth_frame"]
    assert "C30" in provenance["gravity"]
    accesses, contacts, targets = analyze_accesses(ephemeris, sun, mission)
    assert len(accesses) > 0
    assert len(targets) == 273
    for table in (accesses, contacts):
        assert table.start_s.ge(0).all()
        assert table.end_s.le(mission.duration_days*86400).all()
        assert table.duration_s.gt(0).all()
