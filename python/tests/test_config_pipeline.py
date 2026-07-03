from pathlib import Path

import numpy as np

from eo_mission.analysis import best_solution_per_gsd, rank_feasible_solutions
from eo_mission.config import load_config
from eo_mission.coverage import OrekitCoverageAdapter, compute_coverage_tensor
from eo_mission.mass import compute_mass_tensor
from eo_mission.optics import compute_optics_tensors


CONFIG = Path(__file__).parents[1] / "config" / "portfolio.yaml"


def test_yaml_config_parses_grid_and_catalogs():
    config = load_config(CONFIG)
    assert config.coverage.mode == "target_lat_revisit"
    assert config.thresholds.snr_min == 400
    assert len(config.mission.gsd_m) == 3
    assert "chroma_d" in config.detector_by_id
    assert config.telescope_by_id["refractive"].max_refractive_diameter_mm == 90


def test_optics_tensors_have_labeled_shapes():
    config = load_config(CONFIG)
    small = _small_config(config)
    snr, mtf = compute_optics_tensors(small)
    assert snr.axes == ("gsd", "band", "detector", "telescope", "altitude", "diameter")
    assert snr.values.shape == mtf.values.shape
    assert snr.values.shape[-2:] == (2, 3)


def test_coverage_masks_detector_and_fov_constraints():
    config = _small_config(load_config(CONFIG))
    adapter = OrekitCoverageAdapter(cache_dir=Path("/tmp/eo-mission-test-cache"))
    coverage = compute_coverage_tensor(config, adapter=adapter)
    assert coverage.values.shape[-1] == 3
    assert np.isfinite(coverage.values).any()


def test_ranked_feasible_solution_selection_per_gsd():
    config = _small_config(load_config(CONFIG))
    snr, mtf = compute_optics_tensors(config)
    coverage = compute_coverage_tensor(
        config, adapter=OrekitCoverageAdapter(cache_dir=Path("/tmp/eo-mission-test-cache"))
    )
    mass = compute_mass_tensor(config)
    ranked = rank_feasible_solutions(config, snr, mtf, coverage, mass)
    best = best_solution_per_gsd(ranked)
    assert not ranked.empty
    assert set(best["gsd_m"]) == set(config.mission.gsd_m)
    assert (best["rank_within_gsd"] == 1).all()


def _small_config(config):
    from dataclasses import replace

    band = replace(config.bands[0], detector_ids=tuple(det.id for det in config.detectors[:2]))
    return replace(
        config,
        mission=replace(
            config.mission,
            altitudes_km=(500.0, 600.0),
            diameters_mm=(80.0, 120.0, 180.0),
            swaths_km=(30.0, 50.0, 80.0),
            gsd_m=(50.0, 80.0),
        ),
        bands=(band,),
        detectors=config.detectors[:2],
        telescopes=config.telescopes[:2],
        constellations=config.constellations[:2],
    )
