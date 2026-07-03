"""Compatibility wrapper for the package CLI."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from eo_mission.cli import main


if __name__ == "__main__":
    main()
