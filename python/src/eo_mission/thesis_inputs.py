"""Read the literal design inputs in MATLAB's master.m without executing MATLAB."""
from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np


@dataclass(frozen=True)
class ThesisInputs:
    gsd_m: float
    wavelengths_m: tuple[float, ...]
    altitudes_km: tuple[float, ...]
    diameters_mm: tuple[float, ...]
    swaths_km: tuple[float, ...]
    snr_min: float
    mtf_min: float
    revisit_days: float
    max_detectors: int
    overlap: float
    cloud_fraction: float
    ltan_h: float
    maintenance_lower: float
    maintenance_upper: float
    constellations: tuple[tuple[int, int], ...]
    co2: dict
    o2: dict
    telescope_names: tuple[str, ...]
    alignment: tuple[float, ...]
    transmission: tuple[float, ...]
    fov_deg: tuple[float, ...]
    obstruction: tuple[float, ...]
    source_path: str


def _scalar(expression: str) -> float:
    """Accept only numeric constants and elementary arithmetic, never eval/code."""
    def visit(node):
        if isinstance(node, ast.Constant) and isinstance(node.value, (float, int)):
            return node.value
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            return visit(node.operand) * (-1 if isinstance(node.op, ast.USub) else 1)
        if isinstance(node, ast.BinOp):
            a, b = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add): return a + b
            if isinstance(node.op, ast.Sub): return a - b
            if isinstance(node.op, ast.Mult): return a * b
            if isinstance(node.op, ast.Div): return a / b
        raise ValueError(f"Unsupported MATLAB numeric literal: {expression}")
    result = float(visit(ast.parse(expression.strip(), mode="eval").body))
    if not np.isfinite(result): raise ValueError("Design inputs must be finite")
    return result


def _literal(expression: str):
    expression = expression.strip()
    if expression.startswith("{"):
        return tuple(re.findall(r"'([^']*)'", expression))
    if expression.startswith("["):
        match = re.fullmatch(r"\[([^]]*)\]\s*(?:\*\s*(.*))?", expression)
        if not match: raise ValueError(f"Unsupported vector: {expression}")
        factor = _scalar(match[2]) if match[2] else 1
        rows = tuple(tuple(_scalar(v) * factor for v in re.split(r"[,\s]+", row.strip()) if v)
                     for row in match[1].split(";"))
        return rows if len(rows) > 1 else rows[0]
    if ":" in expression:
        parts = [_scalar(v) for v in expression.split(":")]
        start, step, stop = (parts[0], 1, parts[1]) if len(parts) == 2 else parts
        if step <= 0 or stop < start: raise ValueError("Invalid design range")
        return tuple(float(v) for v in np.arange(start, stop + step * 1e-8, step))
    return _scalar(expression)


def load_thesis_inputs(master_path: str | Path | None = None) -> ThesisInputs:
    master = Path(master_path) if master_path else Path(__file__).resolve().parents[3] / "MATLAB_Code/master.m"
    text = "\n".join(line.split("%", 1)[0] for line in master.read_text(encoding="utf-8").splitlines())
    def read(name):
        # The first assignment is the editable literal; later environment overrides are MATLAB execution settings.
        match = re.search(r"config\." + re.escape(name) + r"\s*=\s*(\[[^]]*\](?:\s*\*\s*[^;]+)?|\{[^}]*\}|[^;]+);", text)
        if not match: raise ValueError(f"Missing literal config.{name} in {master}")
        return _literal(match[1])
    def detector(prefix):
        return {key: read(prefix + "." + field) for key, field in {
            "names": "names", "eta": "quantumEfficiency", "mtf": "detectorMtf",
            "pixels": "pixelCount", "pitch_m": "pixelSize_m"}.items()}
    inputs = ThesisInputs(
        read("groundSampleDistance_m"), read("wavelengths_m"), read("altitudes_km"),
        read("pupilDiameters_mm"), read("swaths_km"), read("snrRequirement"),
        read("mtfRequirement"), read("revisitRequirement_days"), int(read("maxDetectorCount")),
        read("overlapFraction"), read("cloudFraction"), read("localTimeAscendingNode_h"),
        read("maintenanceLowerFactor"), read("maintenanceUpperFactor"),
        tuple(tuple(int(v) for v in row) for row in read("satelliteConfigurations")),
        detector("co2Detectors"), detector("o2Detectors"), read("telescopeNames"),
        read("telescopeAlignmentMtf"), read("telescopeTransmittance"),
        read("telescopeFov_deg"), read("obstructionRatios"), str(master.resolve()),
    )
    if not 0 <= inputs.cloud_fraction < 1 or inputs.gsd_m <= 0:
        raise ValueError("Invalid cloud fraction or GSD")
    return inputs
