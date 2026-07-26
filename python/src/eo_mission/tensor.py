"""Small labeled-array helper for mission design tensors."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class LabeledTensor:
    """NumPy array with named axes and coordinate labels."""

    values: np.ndarray
    axes: tuple[str, ...]
    coords: Mapping[str, np.ndarray]

    def axis(self, name: str) -> int:
        return self.axes.index(name)

    def coord(self, name: str) -> np.ndarray:
        return np.asarray(self.coords[name])


def as_array(values: Sequence[float]) -> np.ndarray:
    return np.asarray(values, dtype=float)
