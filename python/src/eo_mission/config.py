"""YAML configuration for the EO mission portfolio pipeline."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class Band:
    id: str
    wavelength_m: float
    detector_ids: tuple[str, ...]


@dataclass(frozen=True)
class Detector:
    id: str
    name: str
    eta: float
    mtf_detector: float
    n_pixels: int
    pixel_size_m: float
    max_swath_multiplier: float = 1.0


@dataclass(frozen=True)
class Telescope:
    id: str
    name: str
    mtf_alignment: float
    transmission: float
    obstruction_ratio: float
    fov_limit_deg: float
    max_refractive_diameter_mm: float | None = None


@dataclass(frozen=True)
class Constellation:
    id: str
    satellites: int
    telescopes_per_satellite: int


@dataclass(frozen=True)
class MissionGrid:
    altitudes_km: tuple[float, ...]
    diameters_mm: tuple[float, ...]
    swaths_km: tuple[float, ...]
    gsd_m: tuple[float, ...]


@dataclass(frozen=True)
class Thresholds:
    snr_min: float
    mtf_min: float
    revisit_days_max: float


@dataclass(frozen=True)
class CoverageSettings:
    mode: str = "conus_grid"
    target_latitude_deg: float = 39.5
    cloud_cover_fraction: float = 0.30
    cloud_persistence: float = 0.5
    overlap_fraction: float = 0.0
    mc_samples: int = 200
    mc_seed: int = 12345
    revisit_percentile: float = 95.0
    propagation_days: int = 30
    grid_lat_start: float = 25.0
    grid_lat_stop: float = 49.0
    grid_lat_step: float = 4.0
    grid_lon_samples: int = 5
    us_width_km: float = 4100.0
    cache_dir: str = ".cache/orekit-coverage"
    orekit_data_path: str | None = None
    parallel_workers: int = 1


@dataclass(frozen=True)
class PipelineConfig:
    mission: MissionGrid
    bands: tuple[Band, ...]
    detectors: tuple[Detector, ...]
    telescopes: tuple[Telescope, ...]
    constellations: tuple[Constellation, ...]
    thresholds: Thresholds
    coverage: CoverageSettings

    @property
    def detector_by_id(self) -> dict[str, Detector]:
        return {det.id: det for det in self.detectors}

    @property
    def telescope_by_id(self) -> dict[str, Telescope]:
        return {tel.id: tel for tel in self.telescopes}


def _expand_range(raw: Any) -> tuple[float, ...]:
    if isinstance(raw, dict):
        start = float(raw["start"])
        stop = float(raw["stop"])
        step = float(raw["step"])
        values: list[float] = []
        x = start
        while x <= stop + step * 1e-9:
            values.append(round(x, 10))
            x += step
        return tuple(values)
    return tuple(float(v) for v in raw)


def load_config(path: str | Path) -> PipelineConfig:
    with Path(path).open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)

    mission_raw = raw["mission"]
    mission = MissionGrid(
        altitudes_km=_expand_range(mission_raw["altitudes_km"]),
        diameters_mm=_expand_range(mission_raw["diameters_mm"]),
        swaths_km=_expand_range(mission_raw["swaths_km"]),
        gsd_m=_expand_range(mission_raw["gsd_m"]),
    )
    bands = tuple(
        Band(
            id=str(item["id"]),
            wavelength_m=float(item["wavelength_m"]),
            detector_ids=tuple(str(v) for v in item["detector_ids"]),
        )
        for item in raw["bands"]
    )
    detectors = tuple(
        Detector(
            id=str(item["id"]),
            name=str(item["name"]),
            eta=float(item["eta"]),
            mtf_detector=float(item["mtf_detector"]),
            n_pixels=int(item["n_pixels"]),
            pixel_size_m=float(item["pixel_size_m"]),
            max_swath_multiplier=float(item.get("max_swath_multiplier", 1.0)),
        )
        for item in raw["detectors"]
    )
    telescopes = tuple(
        Telescope(
            id=str(item["id"]),
            name=str(item["name"]),
            mtf_alignment=float(item["mtf_alignment"]),
            transmission=float(item["transmission"]),
            obstruction_ratio=float(item["obstruction_ratio"]),
            fov_limit_deg=float(item["fov_limit_deg"]),
            max_refractive_diameter_mm=(
                float(item["max_refractive_diameter_mm"])
                if item.get("max_refractive_diameter_mm") is not None
                else None
            ),
        )
        for item in raw["telescopes"]
    )
    constellations = tuple(
        Constellation(
            id=str(item["id"]),
            satellites=int(item["satellites"]),
            telescopes_per_satellite=int(item["telescopes_per_satellite"]),
        )
        for item in raw["constellations"]
    )
    thresholds_raw = raw["thresholds"]
    thresholds = Thresholds(
        snr_min=float(thresholds_raw["snr_min"]),
        mtf_min=float(thresholds_raw["mtf_min"]),
        revisit_days_max=float(thresholds_raw["revisit_days_max"]),
    )
    coverage_raw = raw.get("coverage", {})
    coverage = CoverageSettings(
        mode=str(coverage_raw.get("mode", "conus_grid")),
        target_latitude_deg=float(coverage_raw.get("target_latitude_deg", 39.5)),
        cloud_cover_fraction=float(coverage_raw.get("cloud_cover_fraction", 0.30)),
        cloud_persistence=float(coverage_raw.get("cloud_persistence", 0.5)),
        overlap_fraction=float(coverage_raw.get("overlap_fraction", 0.0)),
        mc_samples=int(coverage_raw.get("mc_samples", 200)),
        mc_seed=int(coverage_raw.get("mc_seed", 12345)),
        revisit_percentile=float(coverage_raw.get("revisit_percentile", 95.0)),
        propagation_days=int(coverage_raw.get("propagation_days", 30)),
        grid_lat_start=float(coverage_raw.get("grid_lat_start", 25.0)),
        grid_lat_stop=float(coverage_raw.get("grid_lat_stop", 49.0)),
        grid_lat_step=float(coverage_raw.get("grid_lat_step", 4.0)),
        grid_lon_samples=int(coverage_raw.get("grid_lon_samples", 5)),
        us_width_km=float(coverage_raw.get("us_width_km", 4100.0)),
        cache_dir=str(coverage_raw.get("cache_dir", ".cache/orekit-coverage")),
        orekit_data_path=coverage_raw.get("orekit_data_path"),
        parallel_workers=int(coverage_raw.get("parallel_workers", 1)),
    )

    return PipelineConfig(
        mission=mission,
        bands=bands,
        detectors=detectors,
        telescopes=telescopes,
        constellations=constellations,
        thresholds=thresholds,
        coverage=coverage,
    )
