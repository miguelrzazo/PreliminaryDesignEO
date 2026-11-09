"""Compatibility wrapper for the package CLI."""
from __future__ import annotations

import sys
import os
from pathlib import Path

# Small three-component geometry operations are slower when BLAS oversubscribes.
# Explicit user settings take precedence; library imports keep their own policy.
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(name,'1')

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from eo_mission.cli import main


if __name__ == "__main__":
    main()
