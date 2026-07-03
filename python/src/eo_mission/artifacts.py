"""Read and write pipeline artifacts."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .tensor import LabeledTensor


def save_tensor(path: str | Path, tensor: LabeledTensor) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    coords = {name: tensor.coord(name).tolist() for name in tensor.axes}
    np.savez_compressed(
        target,
        values=tensor.values,
        axes=np.asarray(tensor.axes, dtype=object),
        coords=json.dumps(coords),
    )


def load_tensor(path: str | Path) -> LabeledTensor:
    data = np.load(Path(path), allow_pickle=True)
    axes = tuple(str(v) for v in data["axes"])
    raw_coords = json.loads(str(data["coords"]))
    coords = {name: np.asarray(values) for name, values in raw_coords.items()}
    return LabeledTensor(values=data["values"], axes=axes, coords=coords)


def write_ranked_tables(ranked: pd.DataFrame, output_dir: str | Path) -> None:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(out / "ranked_solutions.csv", index=False)
    if ranked.empty:
        return
    per_gsd_dir = out / "per_gsd"
    per_gsd_dir.mkdir(exist_ok=True)
    for gsd, frame in ranked.groupby("gsd_m"):
        label = str(gsd).replace(".", "p")
        frame.to_csv(per_gsd_dir / f"solutions_gsd_{label}.csv", index=False)
    try:
        ranked.to_parquet(out / "ranked_solutions.parquet", index=False)
    except Exception:
        pass
