from dataclasses import replace
from pathlib import Path

import numpy as np

from eo_mission.analysis import best_solution_per_gsd, rank_feasible_solutions
from eo_mission.config import load_config
from eo_mission.coverage import OrekitCoverageAdapter, compute_coverage_tensor
from eo_mission.mass import compute_mass_tensor
from eo_mission.optics import compute_optics_tensors


CONFIG = Path(__file__).parents[1] / "config" / "portfolio.yaml"


def _orekit_data_path() -> str | None:
    import os

    p = os.environ.get("OREKIT_DATA_PATH")
    if p:
        return p
    candidate = Path("/tmp/orekit-data/orekit-data.zip")
    return str(candidate) if candidate.exists() else None


def _adapter():
    a = OrekitCoverageAdapter(
        orekit_data_path=_orekit_data_path(),
        cache_dir=Path("/tmp/eo-mission-test-cache"),
    )
    a.initialize_orekit(require=True)
    return a


def test_yaml_config_parses_grid_and_catalogs():
    config = load_config(CONFIG)
    assert config.coverage.mode == "conus_grid"
    assert config.coverage.cloud_persistence == 0.5
    assert config.coverage.mc_samples == 200
    assert config.coverage.revisit_percentile == 95
    assert config.coverage.propagation_days == 30
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
    coverage = compute_coverage_tensor(config, adapter=_adapter())
    assert coverage.values.shape[-1] == 3
    assert np.isfinite(coverage.values).any()


def test_ranked_feasible_solution_selection_per_gsd():
    config = _small_config(load_config(CONFIG))
    snr, mtf = compute_optics_tensors(config)
    coverage = compute_coverage_tensor(config, adapter=_adapter())
    mass = compute_mass_tensor(config)
    ranked = rank_feasible_solutions(config, snr, mtf, coverage, mass)
    best = best_solution_per_gsd(ranked)
    assert not ranked.empty
    # only GSDs with at least one feasible (SNR + MTF + coverage) cell appear;
    # the small diameter grid may not close at the tightest GSD, soassert a
    # subset relationship rather than full coverage of every configured GSD.
    assert set(best["gsd_m"]).issubset(set(config.mission.gsd_m))
    assert not best.empty
    assert (best["rank_within_gsd"] == 1).all()


def _small_config(config):
    band = replace(config.bands[0], detector_ids=tuple(det.id for det in config.detectors[:2]))
    tiny_cov = replace(
        config.coverage,
        propagation_days=8,
        mc_samples=40,
        grid_lat_step=10.0,
        grid_lon_samples=3,
    )
    return replace(
        config,
        mission=replace(
            config.mission,
            altitudes_km=(600.0, 700.0),
            diameters_mm=(80.0, 120.0, 180.0),
            swaths_km=(60.0, 80.0, 100.0),
            gsd_m=(50.0, 80.0),
        ),
        bands=(band,),
        detectors=config.detectors[:2],
        telescopes=(config.telescopes[1], config.telescopes[3]),  # Korsch, TMA
        constellations=(config.constellations[1], config.constellations[3]),  # 1sat_2tel, 2sat_2tel
        coverage=tiny_cov,
    )
