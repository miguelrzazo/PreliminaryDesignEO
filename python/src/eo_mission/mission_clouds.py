"""Persistent daily clouds applied to propagated point-access intervals.

One daily chain is shared by every target and both satellites in a realization.
This assumes perfect spatial cloud correlation over the study region; its
parameters are illustrative and are not calibrated regional weather data.
"""
from __future__ import annotations

from collections.abc import Sequence
from datetime import datetime, timezone
import math

import numpy as np
import pandas as pd

from .clouds import sample_cloud_days, transition_matrix

DAY_S = 86400.0


def _utc_day_offset(epoch: str) -> float:
    instant = datetime.fromisoformat(epoch.replace("Z", "+00:00"))
    if instant.utcoffset() is None or instant.utcoffset().total_seconds() != 0:
        raise ValueError("The cloud epoch must explicitly use UTC.")
    utc = instant.astimezone(timezone.utc)
    midnight = utc.replace(hour=0, minute=0, second=0, microsecond=0)
    return (utc - midnight).total_seconds()


def _merge(intervals: list[tuple[float, float]]) -> list[tuple[float, float]]:
    merged: list[tuple[float, float]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def split_intervals_at_utc_midnight(
    intervals: Sequence[tuple[float, float]], horizon_s: float,
    epoch: str = "2024-04-28T00:00:00Z",
) -> list[tuple[int, float, float]]:
    """Return (UTC-day index, start seconds, end seconds), clipped to the horizon."""
    if not math.isfinite(horizon_s) or horizon_s <= 0:
        raise ValueError("The horizon must be positive and finite.")
    offset = _utc_day_offset(epoch)
    clipped = []
    for start, end in intervals:
        if not math.isfinite(start) or not math.isfinite(end) or end < start:
            raise ValueError("Access intervals must have finite, ordered endpoints.")
        start, end = max(0.0, float(start)), min(horizon_s, float(end))
        if end > start:
            clipped.append((start, end))
    chunks = []
    for start, end in _merge(clipped):
        day = math.floor((offset + start) / DAY_S)
        while start < end:
            boundary = (day + 1) * DAY_S - offset
            chunk_end = min(end, boundary)
            if chunk_end > start:
                chunks.append((day, start, chunk_end))
            start = chunk_end
            day += 1
    return chunks


def _gap_metrics(intervals: list[tuple[float, float]], horizon_s: float) -> tuple[float, float | None, bool]:
    intervals = _merge(intervals)
    if not intervals:
        return horizon_s / DAY_S, None, True
    interior = [right[0] - left[1] for left, right in zip(intervals, intervals[1:])]
    maximum = max(intervals[0][0], horizon_s - intervals[-1][1], *interior)
    return maximum / DAY_S, max(interior) / DAY_S if interior else None, False


def _distribution(values: Sequence[float]) -> dict:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if not len(finite):
        return {"mean": None, "p50": None, "p95": None, "valid_realizations": 0}
    return {"mean": float(finite.mean()), "p50": float(np.percentile(finite, 50)),
            "p95": float(np.percentile(finite, 95)), "valid_realizations": int(len(finite))}


def cloud_hardware_floor(
    horizon_s: float, *, cloud_fraction: float = 0.30, persistence: float = 0.5,
    mc_samples: int = 200, seed: int = 12345,
    epoch: str = "2024-04-28T00:00:00Z", requirement_days: float = 7.0,
) -> dict:
    """Bound any shared-cloud design by its longest fully cloudy time interval.

    No constellation or optical design can acquire a target during a regionwide
    cloudy interval. Each target's boundary-inclusive maximum no-access gap is
    therefore at least as long as that interval. The bound applies to these
    exact Monte Carlo chains and cloud semantics, not calibrated climate risk.
    """
    transition_matrix(cloud_fraction, persistence)
    if not math.isfinite(horizon_s) or horizon_s <= 0:
        raise ValueError("The horizon must be positive and finite.")
    if not math.isfinite(requirement_days) or requirement_days <= 0:
        raise ValueError("The requirement must be positive and finite.")
    if isinstance(mc_samples, bool) or not isinstance(mc_samples, (int, np.integer)) or mc_samples < 1:
        raise ValueError("mc_samples must be a positive integer.")
    offset = _utc_day_offset(epoch)
    n_days = math.ceil((offset + horizon_s) / DAY_S)
    clear = sample_cloud_days(n_days, cloud_fraction, persistence,
                              np.random.default_rng(seed), n_samples=mc_samples)
    day_starts = np.maximum(0.0, np.arange(n_days) * DAY_S - offset)
    day_ends = np.minimum(horizon_s, (np.arange(n_days) + 1) * DAY_S - offset)
    floors = []
    for realization in clear:
        cloudy = [(float(start), float(end)) for start, end, is_clear in zip(day_starts, day_ends, realization)
                  if not is_clear and end > start]
        merged = _merge(cloudy)
        floors.append(max((end - start for start, end in merged), default=0.0) / DAY_S)
    distribution = _distribution(floors)
    return {
        "parameters": {"cloud_fraction": float(cloud_fraction), "persistence": float(persistence),
                       "mc_samples": int(mc_samples), "seed": int(seed), "epoch": epoch,
                       "horizon_s": float(horizon_s), "requirement_days": float(requirement_days)},
        "regional_max_gap_floor_days": distribution,
        "per_realization_floor_days": floors,
        "empirical_p95_floor_exceeds_requirement": distribution["p95"] > requirement_days,
        "probability_cloud_floor_le_requirement": float(np.mean(np.asarray(floors) <= requirement_days)),
        "definition": "Longest contiguous regionwide cloudy interval inside the horizon; a lower bound on every target's maximum no-clear-access gap including censored edges. All hardware designs use these same cloud chains.",
        "scope": "Empirical bound for the configured shared-spatial-cloud Monte Carlo objective. If its p95 exceeds the requirement, no hardware design can meet that empirical p95 objective. This does not establish an infinite-horizon or calibrated meteorological probability.",
    }


def apply_markov_clouds(
    acquisitions: pd.DataFrame, targets: pd.DataFrame | Sequence[str], horizon_s: float,
    *, cloud_fraction: float = 0.30, persistence: float = 0.5,
    mc_samples: int = 200, seed: int = 12345,
    epoch: str = "2024-04-28T00:00:00Z", requirement_days: float = 7.0,
) -> tuple[dict, pd.DataFrame]:
    """Retain clear parts of actual accesses and summarize worst-target gaps.

    Gaps are measured between interval endpoints in elapsed seconds, including
    both horizon boundaries. The observed interior statistic excludes those
    censored edge gaps and is undefined when no target has two clear windows.
    """
    transition_matrix(cloud_fraction, persistence)
    if not math.isfinite(horizon_s) or horizon_s <= 0:
        raise ValueError("The horizon must be positive and finite.")
    if not math.isfinite(requirement_days) or requirement_days <= 0:
        raise ValueError("The requirement must be positive and finite.")
    if isinstance(mc_samples, bool) or not isinstance(mc_samples, (int, np.integer)) or mc_samples < 1:
        raise ValueError("mc_samples must be a positive integer.")
    identifiers = list(targets["target"] if isinstance(targets, pd.DataFrame) else targets)
    if not identifiers or len(set(identifiers)) != len(identifiers) or any(pd.isna(x) for x in identifiers):
        raise ValueError("Provide nonempty, unique target identifiers.")
    required = {"target", "start_s", "end_s"}
    if not required.issubset(acquisitions.columns):
        raise ValueError("Acquisitions require target, start_s and end_s columns.")
    if not set(acquisitions.target).issubset(identifiers):
        raise ValueError("An acquisition refers to an unknown target.")
    offset = _utc_day_offset(epoch)
    n_days = math.ceil((offset + horizon_s) / DAY_S)
    clear_days = sample_cloud_days(n_days, cloud_fraction, persistence,
                                  np.random.default_rng(seed), n_samples=mc_samples)
    grouped = acquisitions.groupby("target", sort=False)
    chunks = []
    for identifier in identifiers:
        events = grouped.get_group(identifier) if identifier in grouped.groups else acquisitions.iloc[:0]
        chunks.append(split_intervals_at_utc_midnight(
            list(zip(events.start_s, events.end_s)), horizon_s, epoch))

    baseline = [_gap_metrics([(start, end) for _, start, end in target], horizon_s) for target in chunks]
    records = []
    for realization, clear in enumerate(clear_days):
        target_results = [_gap_metrics([(start, end) for day, start, end in target if clear[day]], horizon_s)
                          for target in chunks]
        maximum = max(result[0] for result in target_results)
        interior = [result[1] for result in target_results if result[1] is not None]
        unobserved = sum(result[2] for result in target_results)
        count = sum(len(_merge([(start, end) for day, start, end in target if clear[day]])) for target in chunks)
        records.append({"realization": realization, "regional_max_gap_days": maximum,
                        "regional_observed_interior_max_gap_days": max(interior) if interior else None,
                        "unobserved_target_fraction": unobserved / len(identifiers),
                        "unobserved_target_count": unobserved, "clear_access_count": count,
                        "all_targets_observed": unobserved == 0,
                        "sampled_regional_max_gap_le_requirement": maximum <= requirement_days})
    frame = pd.DataFrame(records)
    summary = {
        "parameters": {"cloud_fraction": float(cloud_fraction), "persistence": float(persistence),
                       "mc_samples": int(mc_samples), "seed": int(seed), "epoch": epoch,
                       "horizon_s": float(horizon_s), "requirement_days": float(requirement_days),
                       "utc_calendar_days": n_days, "target_count": len(identifiers)},
        "regional_max_gap_days": _distribution(frame.regional_max_gap_days),
        "regional_observed_interior_max_gap_days": _distribution(frame.regional_observed_interior_max_gap_days),
        "unobserved_target_fraction": _distribution(frame.unobserved_target_fraction),
        "probability_regional_max_gap_le_requirement": float(frame.sampled_regional_max_gap_le_requirement.mean()),
        "probability_all_targets_observed": float(frame.all_targets_observed.mean()),
        "probability_all_targets_observed_and_max_gap_le_requirement": float(
            (frame.all_targets_observed & frame.sampled_regional_max_gap_le_requirement).mean()),
        "geometric_baseline": {"regional_max_gap_days": max(result[0] for result in baseline),
                               "unobserved_target_fraction": sum(result[2] for result in baseline) / len(identifiers)},
        "definitions": {
            "cloud_model": "Stationary binary daily Markov chain, clear=False means fully cloudy; one chain shared by all targets and satellites in each realization.",
            "spatial_correlation": "Perfect regionwide daily spatial correlation; illustrative, uncalibrated climate parameters.",
            "intervals": "Actual propagated point accesses, split at UTC midnight and clipped to the horizon; surviving intervals from both satellites are merged per target.",
            "regional_max_gap_days": "Each realization's largest target gap between clear-access endpoints, including left and right censored horizon edges; Monte Carlo percentiles describe these regional maxima.",
            "observed_interior_gap": "Largest fully observed gap between consecutive clear intervals; absent if no target has two intervals. Absence is not zero revisit.",
            "probability": "Monte Carlo fraction with finite-horizon sampled regional maximum gap <= requirement, not a future weekly-coverage probability or mission certification.",
            "censoring": "An unobserved target has a gap equal to the full horizon. If horizon <= requirement, the gap threshold alone is always satisfied, even with no accesses; inspect all-target observation probability too.",
        },
    }
    return summary, frame
