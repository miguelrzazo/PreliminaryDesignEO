import os
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from eo_mission.config import load_config
from eo_mission.coverage import OrekitCoverageAdapter, compute_coverage_tensor


pytestmark = pytest.mark.skipif(
    os.environ.get("EO_MISSION_RUN_OREKIT_TESTS") != "1",
    reason="Orekit integration tests run only in the Docker Orekit runtime",
)


def test_orekit_tiny_sso_coverage_grid():
    config = load_config(Path(__file__).parents[1] / "config" / "portfolio.yaml")
    config = replace(
        config,
        mission=replace(
            config.mission,
            altitudes_km=(550.0,),
            swaths_km=(180.0,),
            diameters_mm=(120.0,),
            gsd_m=(80.0,),
        ),
        bands=(replace(config.bands[0], detector_ids=(config.detectors[1].id,)),),
        detectors=(config.detectors[1],),
        telescopes=(config.telescopes[1],),
        constellations=(config.constellations[-1],),
    )
    adapter = OrekitCoverageAdapter(
        orekit_data_path=os.environ.get("OREKIT_DATA_PATH"),
        cache_dir=Path("/tmp/eo-mission-orekit-cache"),
    )
    adapter.initialize_orekit(require=True)
    coverage = compute_coverage_tensor(config, adapter=adapter)
    assert adapter.jvm_initialized
    assert coverage.values.shape == (1, 1, 1, 1, 1, 1)
    assert np.isfinite(coverage.values).any()
