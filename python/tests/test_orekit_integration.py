"""Orekit integration tests for the SSO coverage computation.

These tests run the real Orekit JVM (no longer gated behind an env var): the
coverage module is now Orekit-only. Set ``OREKIT_DATA_PATH`` to the orekit-data
zip or directory if not auto-discovered.
"""
import os
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from eo_mission.config import load_config
from eo_mission.coverage import (
    OrekitCoverageAdapter,
    compute_coverage_tensor,
    sso_inclination_deg,
)


def _orekit_data_path() -> str | None:
    p = os.environ.get("OREKIT_DATA_PATH")
    if p:
        return p
    # fall back to the well-known location used during environment setup
    candidate = Path("/tmp/orekit-data/orekit-data.zip")
    return str(candidate) if candidate.exists() else None


@pytest.fixture(scope="module")
def adapter():
    a = OrekitCoverageAdapter(
        orekit_data_path=_orekit_data_path(),
        cache_dir=Path("/tmp/eo-mission-orekit-cache"),
    )
    a.initialize_orekit(require=True)
    return a


# ── SSO inclination feasibility (pure-Python, no JVM) ─────────────────────────

def test_sso_inclination_typical_altitude():
    """520 km must give a retrograde SSO inclination near 97-98°."""
    inc = sso_inclination_deg(520.0)
    assert np.isfinite(inc)
    assert 95.0 < inc < 100.0


def test_sso_inclination_too_high_is_infeasible():
    """Very high altitudes cannot sustain SSO; inclination must be NaN.

    The J2 sun-sync condition has a solution up to roughly 5000 km; above it the
    required nodal precession exceeds what J2 can supply and |cos_i| > 1.
    """
    inc = sso_inclination_deg(6000.0)
    assert np.isnan(inc)


def test_sso_inclination_increases_with_altitude():
    """SSO inclination grows toward polar as altitude increases."""
    incs = [sso_inclination_deg(h) for h in (400.0, 600.0, 800.0)]
    incs = [i for i in incs if np.isfinite(i)]
    assert np.all(np.diff(incs) > 0)


# ── Orekit coverage tensor ────────────────────────────────────────────────────

def test_orekit_tiny_sso_coverage_grid(adapter):
    """A single feasible cell must produce a finite revisit ≤ the threshold."""
    config = load_config(Path(__file__).parents[1] / "config" / "portfolio.yaml")
    config = replace(
        config,
        mission=replace(
            config.mission,
            altitudes_km=(700.0,),
            swaths_km=(80.0,),
            diameters_mm=(120.0,),
            gsd_m=(80.0,),
        ),
        bands=(replace(config.bands[0], detector_ids=(config.detectors[1].id,)),),
        detectors=(config.detectors[1],),
        telescopes=(config.telescopes[3],),  # TMA (8° FOV ⇒ wide swath feasible)
        constellations=(config.constellations[3],),  # 2 sat × 2 telescopes
        coverage=replace(
            config.coverage,
            propagation_days=10,
            mc_samples=60,
            grid_lat_step=8.0,
            grid_lon_samples=3,
        ),
    )
    coverage = compute_coverage_tensor(config, adapter=adapter)
    assert adapter.jvm_initialized
    assert coverage.values.shape == (1, 1, 1, 1, 1, 1)
    assert np.isfinite(coverage.values).any()
    assert float(coverage.values[0, 0, 0, 0, 0, 0]) <= 7.0


def test_more_satellites_shorter_revisit(adapter):
    """A 2-satellite constellation must revisit no slower than 1 satellite."""
    base = load_config(Path(__file__).parents[1] / "config" / "portfolio.yaml")
    common = dict(
        propagation_days=10, mc_samples=60, grid_lat_step=8.0, grid_lon_samples=3,
    )
    cfg1 = replace(base, mission=replace(base.mission, altitudes_km=(700.0,),
                                         swaths_km=(80.0,), diameters_mm=(120.0,),
                                         gsd_m=(80.0,)),
                   bands=(replace(base.bands[0], detector_ids=(base.detectors[1].id,)),),
                   detectors=(base.detectors[1],), telescopes=(base.telescopes[3],),
                   constellations=(base.constellations[1],),  # 1 sat × 2 tel
                   coverage=replace(base.coverage, **common))
    cfg2 = replace(cfg1, constellations=(base.constellations[3],))  # 2 sat × 2 tel
    cov1 = compute_coverage_tensor(cfg1, adapter=adapter)
    cov2 = compute_coverage_tensor(cfg2, adapter=adapter)
    v1, v2 = float(cov1.values[0, 0, 0, 0, 0, 0]), float(cov2.values[0, 0, 0, 0, 0, 0])
    assert np.isfinite(v1) and np.isfinite(v2)
    assert v2 <= v1 + 1e-9


def test_coverage_nan_beyond_fov(adapter):
    """A swath wider than the telescope FOV limit must yield NaN."""
    base = load_config(Path(__file__).parents[1] / "config" / "portfolio.yaml")
    # TMA has fov_limit_deg=8 ⇒ at 600 km max swath ≈ 84 km; request 160 km ⇒ NaN.
    cfg = replace(base, mission=replace(base.mission, altitudes_km=(600.0,),
                                        swaths_km=(160.0,), diameters_mm=(80.0,),
                                        gsd_m=(80.0,)),
                  bands=(replace(base.bands[0], detector_ids=(base.detectors[1].id,)),),
                  detectors=(base.detectors[1],), telescopes=(base.telescopes[3],),
                  constellations=(base.constellations[0],),
                  coverage=replace(base.coverage, propagation_days=8, mc_samples=40,
                                   grid_lat_step=10.0, grid_lon_samples=3))
    cov = compute_coverage_tensor(cfg, adapter=adapter)
    assert np.isnan(cov.values[0, 0, 0, 0, 0, 0])
