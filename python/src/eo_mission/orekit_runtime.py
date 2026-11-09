"""Load the JCC Orekit wrapper and register its physical data explicitly."""
from __future__ import annotations

import os
from pathlib import Path


def resolve_data_path(path: str | Path | None = None) -> Path:
    """Resolve the data archive/directory without silently accepting a bad path."""
    candidate = path or os.environ.get("OREKIT_DATA_PATH")
    if candidate is None:
        local = Path.cwd() / "orekit-data.zip"
        if local.is_file():
            candidate = local
    if candidate is None:
        raise FileNotFoundError(
            "Orekit physical data are required. Set OREKIT_DATA_PATH to the "
            "orekit-data.zip archive or an unpacked data directory."
        )
    resolved = Path(candidate).expanduser().resolve()
    if not resolved.exists():
        raise FileNotFoundError(f"Orekit data do not exist: {resolved}")
    if resolved.is_dir() and not any(resolved.iterdir()):
        raise ValueError(f"Orekit data directory is empty: {resolved}")
    return resolved


_registered_paths: set[Path] = set()


def initialize_runtime(path: str | Path | None = None) -> Path:
    """Start the JCC JVM, register data and verify that UTC can be constructed."""
    resolved = resolve_data_path(path)
    if _registered_paths and resolved not in _registered_paths:
        raise RuntimeError(
            "Changing Orekit physical-data paths in one JVM is not supported. "
            "Restart Python before selecting a different archive or directory."
        )
    try:
        import orekit
    except ImportError as exc:
        raise RuntimeError(
            "The JCC Orekit wrapper is not installed in this Python environment. "
            "Use the conda-forge environment in python/environment-orekit.yml."
        ) from exc
    if orekit.getVMEnv() is None:
        orekit.initVM()

    from java.io import File
    from org.orekit.data import DataContext, DirectoryCrawler, ZipJarCrawler
    from org.orekit.time import TimeScalesFactory

    if resolved not in _registered_paths:
        provider = DirectoryCrawler(File(str(resolved))) if resolved.is_dir() else ZipJarCrawler(File(str(resolved)))
        DataContext.getDefault().getDataProvidersManager().addProvider(provider)
        _registered_paths.add(resolved)
    try:
        TimeScalesFactory.getUTC()
    except Exception as exc:
        raise RuntimeError(
            f"Cannot load UTC from Orekit data at {resolved}. Check the archive contents."
        ) from exc
    return resolved
