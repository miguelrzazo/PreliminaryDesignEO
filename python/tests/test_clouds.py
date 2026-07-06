"""Unit tests for the Markov-persistent cloud model (pure Python, no JVM)."""
import numpy as np
import pytest

from eo_mission.clouds import (
    transition_matrix,
    sample_cloud_days,
    revisit_statistic,
)


def test_transition_matrix_row_stochastic():
    P = transition_matrix(0.3, 0.5)
    assert P.shape == (2, 2)
    assert np.allclose(P.sum(axis=1), 1.0)
    assert (P >= 0).all() and (P <= 1).all()


def test_transition_matrix_stationary_distribution():
    """The stationary cloudy probability must equal p for any rho."""
    for p, rho in [(0.3, 0.5), (0.5, 0.0), (0.7, 0.9), (0.1, 0.2)]:
        P = transition_matrix(p, rho)
        # stationary π satisfies π P = π; π_cloudy = p
        w, v = np.linalg.eig(P.T)
        idx = int(np.argmax(w.real))
        pi = np.real(v[:, idx])
        pi = np.abs(pi) / np.abs(pi).sum()  # eigenvector sign is arbitrary
        assert np.isclose(pi[1], p, atol=1e-6)


def test_transition_matrix_rejects_invalid_params():
    with pytest.raises(ValueError):
        transition_matrix(-0.1, 0.5)
    with pytest.raises(ValueError):
        transition_matrix(0.3, 1.0)


def test_sample_cloud_days_reproducible():
    rng1 = np.random.default_rng(42)
    rng2 = np.random.default_rng(42)
    a = sample_cloud_days(50, 0.3, 0.5, rng1, n_samples=3)
    b = sample_cloud_days(50, 0.3, 0.5, rng2, n_samples=3)
    assert a.shape == (3, 50)
    np.testing.assert_array_equal(a, b)


def test_sample_cloud_days_stationary_frequency():
    """With rho=0 (iid) the long-run clear-day frequency ≈ 1 - p."""
    rng = np.random.default_rng(0)
    samples = sample_cloud_days(2000, 0.30, 0.0, rng, n_samples=1)
    clear_freq = samples.mean()
    assert abs(clear_freq - 0.70) < 0.05


def test_sample_cloud_days_persistence_clusters_clouds():
    """Higher rho ⇒ longer cloudy spells ⇒ more autocorrelation."""
    rng = np.random.default_rng(0)
    iid = sample_cloud_days(3000, 0.5, 0.0, rng, n_samples=1)[0]
    rng = np.random.default_rng(1)
    pers = sample_cloud_days(3000, 0.5, 0.8, rng, n_samples=1)[0]
    # mean run length of cloudy days
    def mean_run(is_cloudy):
        runs, cur = [], 1
        for i in range(1, len(is_cloudy)):
            if is_cloudy[i] == is_cloudy[i-1] and is_cloudy[i]:
                cur += 1
            else:
                if is_cloudy[i-1]:
                    runs.append(cur)
                cur = 1
        if is_cloudy[-1]:
            runs.append(cur)
        return np.mean(runs) if runs else 0.0
    cloud_iid = ~iid
    cloud_pers = ~pers
    assert mean_run(cloud_pers) > mean_run(cloud_iid)


def test_revisit_statistic_increases_with_cloud_cover():
    """More clouds ⇒ longer 95th-percentile revisit."""
    rng = np.random.default_rng(0)
    # daily opportunities for one target across 20 days
    opps = [np.array([0, 1, 5, 8, 12, 15, 18])]
    rp_low, _ = revisit_statistic(opps, 20, 0.1, 0.5, rng, n_samples=200, percentile=95.0)
    rp_high, _ = revisit_statistic(opps, 20, 0.7, 0.5, rng, n_samples=200, percentile=95.0)
    assert rp_high > rp_low


def test_revisit_statistic_no_opportunities_is_horizon():
    """A target with no passes is never acquired ⇒ revisit = horizon days."""
    rng = np.random.default_rng(0)
    rp, _ = revisit_statistic([np.array([])], 25, 0.3, 0.5, rng, n_samples=50)
    assert rp == pytest.approx(25.0)


def test_revisit_statistic_returns_nan_for_empty_target_list():
    rng = np.random.default_rng(0)
    rp, mean = revisit_statistic([], 25, 0.3, 0.5, rng, n_samples=10)
    assert np.isnan(rp) and np.isnan(mean)


def test_revisit_percentile_is_at_least_mean():
    """The 95th percentile must not undercut the mean."""
    rng = np.random.default_rng(0)
    opps = [np.array([0, 3, 7, 11, 15, 19])]
    rp, mean = revisit_statistic(opps, 20, 0.4, 0.5, rng, n_samples=500, percentile=95.0)
    assert rp >= mean - 1e-9