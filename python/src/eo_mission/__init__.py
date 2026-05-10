"""Earth Observation satellite mission design and optimization."""
from .optics import compute_snr, compute_mtf, plot_heatmap

__all__ = ["compute_snr", "compute_mtf", "plot_heatmap"]
