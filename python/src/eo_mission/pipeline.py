"""Portfolio pipeline orchestration."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .analysis import best_solution_per_gsd, rank_feasible_solutions
from .artifacts import save_tensor, write_ranked_tables
from .config import PipelineConfig, load_config
from .coverage import OrekitCoverageAdapter, compute_coverage_tensor
from .mass import compute_mass_tensor
from .optics import compute_optics_tensors
from .plots import plot_summary_heatmaps


@dataclass(frozen=True)
class PipelineArtifacts:
    output_dir: Path
    ranked: pd.DataFrame
    best_per_gsd: pd.DataFrame


def run_pipeline(
    config_path: str | Path,
    output_dir: str | Path,
    *,
    initialize_orekit: bool = False,
) -> PipelineArtifacts:
    config: PipelineConfig = load_config(config_path)
    out = Path(output_dir)
    tensors_dir = out / "tensors"
    tables_dir = out / "tables"
    out.mkdir(parents=True, exist_ok=True)

    adapter = OrekitCoverageAdapter(
        orekit_data_path=config.coverage.orekit_data_path,
        cache_dir=out / config.coverage.cache_dir,
    )
    if initialize_orekit:
        adapter.initialize_orekit(require=True)

    snr, mtf = compute_optics_tensors(config)
    coverage = compute_coverage_tensor(config, adapter=adapter)
    mass = compute_mass_tensor(config)
    ranked = rank_feasible_solutions(config, snr, mtf, coverage, mass)
    best = best_solution_per_gsd(ranked)

    save_tensor(tensors_dir / "snr.npz", snr)
    save_tensor(tensors_dir / "mtf.npz", mtf)
    save_tensor(tensors_dir / "coverage.npz", coverage)
    save_tensor(tensors_dir / "mass.npz", mass)
    write_ranked_tables(ranked, tables_dir)
    best.to_csv(tables_dir / "best_per_gsd.csv", index=False)
    plot_summary_heatmaps(snr, mtf, coverage, ranked, out)

    return PipelineArtifacts(output_dir=out, ranked=ranked, best_per_gsd=best)
