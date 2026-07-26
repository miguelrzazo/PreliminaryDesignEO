"""Plot helpers for saved pipeline tensors."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .tensor import LabeledTensor


def _heatmap(matrix: np.ndarray, x: np.ndarray, y: np.ndarray, title: str, label: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    img = ax.imshow(
        matrix,
        aspect="auto",
        origin="lower",
        extent=[float(x[0]), float(x[-1]), float(y[0]), float(y[-1])],
        cmap="viridis",
    )
    ax.set_xlabel("Diameter or swath")
    ax.set_ylabel("Altitude (km)")
    ax.set_title(title)
    fig.colorbar(img, ax=ax, label=label)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_summary_heatmaps(
    snr: LabeledTensor,
    mtf: LabeledTensor,
    coverage: LabeledTensor,
    ranked: pd.DataFrame,
    output_dir: str | Path,
) -> None:
    out = Path(output_dir)
    gsd_i = 0
    band_i = 0
    detector_i = 0
    telescope_i = 0
    constellation_i = 0
    altitudes = snr.coord("altitude")
    diameters = snr.coord("diameter")
    swaths = coverage.coord("swath")

    _heatmap(
        snr.values[gsd_i, band_i, detector_i, telescope_i],
        diameters,
        altitudes,
        "SNR feasible tensor slice",
        "SNR",
        out / "heatmaps" / "snr_slice.png",
    )
    _heatmap(
        mtf.values[gsd_i, band_i, detector_i, telescope_i],
        diameters,
        altitudes,
        "MTF feasible tensor slice",
        "MTF",
        out / "heatmaps" / "mtf_slice.png",
    )
    _heatmap(
        coverage.values[gsd_i, constellation_i, telescope_i, detector_i],
        swaths,
        altitudes,
        "Coverage feasible tensor slice",
        "Revisit days",
        out / "heatmaps" / "coverage_slice.png",
    )

    if ranked.empty:
        return
    pivot = ranked.assign(feasible=1).pivot_table(
        index="altitude_km",
        columns="diameter_mm",
        values="feasible",
        aggfunc="max",
        fill_value=np.nan,
    )
    _heatmap(
        pivot.to_numpy(dtype=float),
        pivot.columns.to_numpy(dtype=float),
        pivot.index.to_numpy(dtype=float),
        "Feasible solution region",
        "Feasible",
        out / "heatmaps" / "feasible_region.png",
    )
