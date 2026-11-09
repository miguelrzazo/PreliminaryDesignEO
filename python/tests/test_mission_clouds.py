"""Cloud filtering preserves actual access timing and midnight boundaries."""
import json

import numpy as np
import pandas as pd
import pytest

from eo_mission.mission_clouds import apply_markov_clouds, cloud_hardware_floor, split_intervals_at_utc_midnight


def events(rows):
    return pd.DataFrame(rows, columns=["satellite", "target", "start_s", "end_s"])


def test_zero_clouds_match_merged_geometric_accesses_exactly():
    acquisitions = events([(1, "A", 10.25, 20.5), (2, "A", 18, 25.75),
                           (1, "A", 40.5, 45), (2, "B", 30.25, 35.5)])
    summary, records = apply_markov_clouds(acquisitions, ["A", "B"], 60, cloud_fraction=0, mc_samples=4)
    assert summary["geometric_baseline"]["regional_max_gap_days"] == pytest.approx(30.25 / 86400)
    assert np.allclose(records.regional_max_gap_days, 30.25 / 86400)
    assert records.clear_access_count.tolist() == [3] * 4
    assert records.regional_observed_interior_max_gap_days.tolist() == pytest.approx([14.75 / 86400] * 4)
    assert summary["probability_all_targets_observed"] == 1
    json.dumps(summary, allow_nan=False)


def test_full_clouds_keep_unobserved_targets_and_do_not_imply_weekly_coverage():
    summary, records = apply_markov_clouds(events([(1, "A", 100, 200)]), ["A", "B"],
                                            7 * 86400, cloud_fraction=1, mc_samples=3)
    assert records.regional_max_gap_days.tolist() == [7] * 3
    assert records.unobserved_target_fraction.tolist() == [1] * 3
    assert records.clear_access_count.tolist() == [0] * 3
    assert summary["regional_observed_interior_max_gap_days"]["mean"] is None
    assert summary["probability_regional_max_gap_le_requirement"] == 1
    assert summary["probability_all_targets_observed_and_max_gap_le_requirement"] == 0
    json.dumps(summary, allow_nan=False)


def test_midnight_splits_and_merges_both_satellites_without_rounding_seconds():
    chunks = split_intervals_at_utc_midnight([(86390.5, 86410.25), (86405, 86420.5)], 2 * 86400)
    assert chunks == [(0, 86390.5, 86400), (1, 86400, 86420.5)]


def test_cloudy_after_midnight_removes_only_cloudy_piece(monkeypatch):
    def deterministic(n_days, p, rho, rng, n_samples):
        assert n_days == 2
        return np.tile([True, False], (n_samples, 1))
    monkeypatch.setattr("eo_mission.mission_clouds.sample_cloud_days", deterministic)
    acquisitions = events([(1, "A", 86390.5, 86410.25), (2, "A", 86405, 86420.5)])
    summary, records = apply_markov_clouds(acquisitions, ["A"], 86430.25, mc_samples=2)
    assert records.clear_access_count.tolist() == [1, 1]
    assert summary["regional_max_gap_days"]["mean"] == pytest.approx(86390.5 / 86400)


def test_calendar_midnight_uses_epoch_time_not_elapsed_day_boundary():
    chunks = split_intervals_at_utc_midnight([(0.25, 20.5)], 30, "2024-04-28T23:59:50Z")
    assert chunks == [(0, 0.25, 10), (1, 10, 20.5)]


def test_cloud_chain_is_shared_regionwide(monkeypatch):
    def deterministic(n_days, p, rho, rng, n_samples):
        return np.tile([False, True], (n_samples, 1))
    monkeypatch.setattr("eo_mission.mission_clouds.sample_cloud_days", deterministic)
    acquisitions = events([(1, "A", 10, 20), (2, "B", 10, 20)])
    _, records = apply_markov_clouds(acquisitions, ["A", "B"], 2 * 86400, mc_samples=2)
    assert records.unobserved_target_count.tolist() == [2, 2]


def test_reproducible_seed_and_fraction_sensitivity():
    acquisitions = events([(1, "A", day * 86400 + 100.125, day * 86400 + 200.375) for day in range(30)])
    first, records = apply_markov_clouds(acquisitions, ["A"], 30 * 86400)
    second, repeat = apply_markov_clouds(acquisitions, ["A"], 30 * 86400)
    assert first == second
    pd.testing.assert_frame_equal(records, repeat)
    lower, _ = apply_markov_clouds(acquisitions, ["A"], 30 * 86400, cloud_fraction=1 / 6)
    assert lower["parameters"]["cloud_fraction"] == 1 / 6
    assert lower["regional_max_gap_days"]["mean"] <= first["regional_max_gap_days"]["mean"]


def test_interval_clipping_and_missing_targets_are_explicit():
    summary, records = apply_markov_clouds(events([(1, "A", -10, 20), (2, "A", 40, 100)]),
                                         pd.DataFrame({"target": ["A", "B"]}), 60, cloud_fraction=0, mc_samples=1)
    assert records.unobserved_target_fraction.iloc[0] == 0.5
    assert summary["regional_max_gap_days"]["mean"] == 60 / 86400


@pytest.mark.parametrize("kwargs", [{"horizon_s": 0}, {"cloud_fraction": -1},
                                   {"persistence": 1}, {"mc_samples": 0},
                                   {"requirement_days": float("nan")},
                                   {"epoch": "2024-04-28T00:00:00"}])
def test_invalid_parameters_fail_explicitly(kwargs):
    arguments = {"horizon_s": 60, **kwargs}
    with pytest.raises(ValueError):
        apply_markov_clouds(events([(1, "A", 10, 20)]), ["A"], **arguments)


def test_unknown_target_and_invalid_intervals_fail_explicitly():
    with pytest.raises(ValueError, match="unknown target"):
        apply_markov_clouds(events([(1, "B", 10, 20)]), ["A"], 60)
    with pytest.raises(ValueError, match="ordered endpoints"):
        apply_markov_clouds(events([(1, "A", 20, 10)]), ["A"], 60)


def test_cloud_hardware_floor_zero_and_full_clouds():
    clear = cloud_hardware_floor(30 * 86400, cloud_fraction=0, mc_samples=3)
    cloudy = cloud_hardware_floor(30 * 86400, cloud_fraction=1, mc_samples=3)
    assert clear["per_realization_floor_days"] == [0, 0, 0]
    assert cloudy["per_realization_floor_days"] == [30, 30, 30]
    assert cloudy["empirical_p95_floor_exceeds_requirement"]
    json.dumps(cloudy, allow_nan=False)


def test_cloud_floor_clips_partial_utc_days(monkeypatch):
    def deterministic(n_days, p, rho, rng, n_samples):
        assert n_days == 3
        return np.tile([False, False, True], (n_samples, 1))
    monkeypatch.setattr("eo_mission.mission_clouds.sample_cloud_days", deterministic)
    result = cloud_hardware_floor(2 * 86400, epoch="2024-04-28T23:59:50Z", mc_samples=1)
    assert result["per_realization_floor_days"] == [(86400 + 10) / 86400]


def test_shared_cloud_floor_bounds_actual_access_maxima_for_same_seed():
    horizon = 30 * 86400
    acquisitions = events([(1, "A", day * 86400 + 100, day * 86400 + 200) for day in range(30)])
    _, records = apply_markov_clouds(acquisitions, ["A"], horizon)
    floor = cloud_hardware_floor(horizon)
    assert np.all(records.regional_max_gap_days.to_numpy() >= floor["per_realization_floor_days"])
    assert floor["regional_max_gap_floor_days"]["p95"] == pytest.approx(11)
    sensitivity = cloud_hardware_floor(horizon, cloud_fraction=1 / 6)
    assert sensitivity["regional_max_gap_floor_days"]["p95"] == pytest.approx(7.05)
