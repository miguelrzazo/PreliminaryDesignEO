"""Earth Observation satellite mission design and optimization."""

__all__ = ["compute_snr", "compute_mtf", "plot_heatmap", "load_config"]


def __getattr__(name):
    if name == "load_config":
        from .config import load_config

        return load_config
    if name in {"compute_snr", "compute_mtf", "plot_heatmap"}:
        from . import optics

        return getattr(optics, name)
    raise AttributeError(name)
