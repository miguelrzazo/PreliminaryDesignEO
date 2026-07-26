"""Probabilistic cloud model for revisit-time accounting.

The MATLAB reference (``CoverageRevisitCalc.m``) inflates the analytical revisit
by a factor ``1 / (1 - cloud_cover)`` — a mean-value correction that ignores the
way cloudy days cluster. Real cloud cover is autocorrelated: a cloudy day is
more likely to be followed by another cloudy day than by a clear one. This
module models that persistence with a two-state Markov chain on the daily cloud
state and estimates, by Monte-Carlo, the revisit you actually achieve at a given
confidence level.

Parameterisation
-----------------
The chain is parameterised by the steady-state cloud-cover probability ``p``
(``cloud_cover_fraction`` from the config) and a *persistence* coefficient
``rho`` in ``[0, 1)`` — the lag-1 autocorrelation of the daily cloud state.
rho = 0 ⇒ independent days (the iid case the MATLAB heuristic approximates);
rho → 1 ⇒ maximally persistent (cloudy spells last many days). From (p, rho)
the transition probabilities are

    P(cloudy → cloudy) = p + rho·(1 - p)
    P(clear  → clear ) = (1 - p) + rho·p

which keeps the stationary distribution exactly P(cloudy) = p for every rho and
guarantees valid probabilities for p ∈ [0, 1], rho ∈ [0, 1).
"""
from __future__ import annotations

import numpy as np


def transition_matrix(p: float, rho: float) -> np.ndarray:
    """Return the 2×2 daily cloud transition matrix.

    States are ordered ``(clear=0, cloudy=1)`` and the matrix is row-stochastic
    ``P[state_t → state_{t+1}]``.
    """
    if not 0.0 <= p <= 1.0:
        raise ValueError(f"cloud_cover_fraction p must be in [0,1], got {p}")
    if not 0.0 <= rho < 1.0:
        raise ValueError(f"cloud_persistence rho must be in [0,1), got {rho}")
    p_cc = p + rho * (1.0 - p)          # cloudy → cloudy
    p_nn = (1.0 - p) + rho * p          # clear  → clear
    return np.array(
        [
            [p_nn, 1.0 - p_nn],         # clear → [clear, cloudy]
            [1.0 - p_cc, p_cc],         # cloudy → [clear, cloudy]
        ],
        dtype=float,
    )


def sample_cloud_days(
    n_days: int,
    p: float,
    rho: float,
    rng: np.random.Generator,
    n_samples: int = 1,
) -> np.ndarray:
    """Simulate ``n_samples`` daily cloud chains of length ``n_days``.

    Parameters
    ----------
    n_days : int
        Number of days in the propagation horizon.
    p, rho : float
        Cloud-cover probability and persistence (see module docstring).
    rng : numpy.random.Generator
        Seeded generator for reproducibility.
    n_samples : int
        Number of independent chains to draw.

    Returns
    -------
    numpy.ndarray
        Boolean array of shape ``(n_samples, n_days)`` where ``True`` means the
        day is *clear* (usable for imaging).
    """
    P = transition_matrix(p, rho)
    states = np.empty((n_samples, n_days), dtype=bool)
    # initial state drawn from the stationary distribution
    states[:, 0] = rng.random(n_samples) >= p
    for t in range(1, n_days):
        prev_cloudy = ~states[:, t - 1]
        # probability of staying/becoming cloudy from the previous state
        p_to_cloudy = np.where(prev_cloudy, P[1, 1], P[0, 1])
        states[:, t] = rng.random(n_samples) >= p_to_cloudy
    return states  # True = clear


def _max_clear_gap_days(
    opportunity_day_index: np.ndarray,
    clear_days: np.ndarray,
    n_days: int,
) -> float:
    """Largest gap (in days) between consecutive *clear* opportunities.

    ``opportunity_day_index`` are the integer day indices (0-based) on which the
    satellite could see a target. ``clear_days`` is the boolean cloud chain for
    one Monte-Carlo run. If no opportunity falls on a clear day, the gap is the
    full horizon (the target was never acquired).
    """
    if opportunity_day_index.size == 0:
        return float(n_days)
    usable = opportunity_day_index[clear_days[opportunity_day_index]]
    if usable.size == 0:
        return float(n_days)
    if usable.size == 1:
        return float(max(usable[0], n_days - 1 - usable[0]))
    gaps = np.diff(usable)
    edge_gap = max(usable[0], n_days - 1 - usable[-1])
    return float(max(gaps.max(), edge_gap))


def revisit_statistic(
    opportunities_per_target: list[np.ndarray],
    n_days: int,
    p: float,
    rho: float,
    rng: np.random.Generator,
    n_samples: int = 200,
    percentile: float = 95.0,
) -> tuple[float, float]:
    """Monte-Carlo worst-target max-clear-gap revisit (single-target metric).

    Parameters
    ----------
    opportunities_per_target : list of 1-D int arrays
        For each ground target, the integer day indices (0-based over the
        horizon) on which at least one satellite has a visibility window.
    n_days : int
        Propagation horizon in days.
    p, rho : float
        Cloud-cover probability and persistence.
    rng : numpy.random.Generator
        Seeded generator.
    n_samples : int
        Number of Monte-Carlo cloud realisations.
    percentile : float
        Confidence level reported (e.g. 95 for the 95th-percentile revisit).

    Returns
    -------
    (revisit_p, revisit_mean) : tuple of float
        The percentile and the mean of the worst-target max-clear-gap across
        runs, both in days.
    """
    if not opportunities_per_target:
        return float("nan"), float("nan")
    cloud = sample_cloud_days(n_days, p, rho, rng, n_samples=n_samples)
    worst = np.empty(n_samples, dtype=float)
    for r in range(n_samples):
        clear = cloud[r]
        worst[r] = max(
            _max_clear_gap_days(opps, clear, n_days)
            for opps in opportunities_per_target
        )
    return float(np.percentile(worst, percentile)), float(worst.mean())


def revisit_statistic_bands(
    band_acquisition_days: list[np.ndarray],
    passes_needed_per_band: list[int],
    n_days: int,
    p: float,
    rho: float,
    rng: np.random.Generator,
    n_samples: int = 200,
    percentile: float = 95.0,
) -> tuple[float, float]:
    """Monte-Carlo latitude-band tiling revisit under persistent clouds.

    This is the coverage semantics of the MATLAB reference
    (``CoverageRevisitCalc.m``): the revisit of a latitude band is the time to
    tile its full longitudinal width, i.e. to accumulate ``passes_needed`` clear
    acquisitions of the band. ``passes_needed = ceil(band_width / swath)`` is the
    same equation the MATLAB code uses (``_US_WIDTH / eff_swath``); here the
    acquisition cadence comes from a real Orekit-propagated ground track instead
    of a hand-coded orbital period, and the ``1/(1-cloud)`` mean heuristic is
    replaced by a Markov-persistent cloud Monte-Carlo at a configurable
    confidence level.

    For each Monte-Carlo run, a daily cloud chain is drawn; the *clear*
    acquisition days of each band are the band's acquisition days that fall on
    clear days. The band revisit is ``passes_needed × n_days / max(1,
    n_clear)`` (the expected wait to gather ``passes_needed`` clear acquisitions
    given a realised clear-acquisition rate). The reported revisit is the
    worst-band value at the ``percentile``-th percentile across runs.

    Returns
    -------
    (revisit_p, revisit_mean) : tuple of float, days
    """
    if not band_acquisition_days:
        return float("nan"), float("nan")
    cloud = sample_cloud_days(n_days, p, rho, rng, n_samples=n_samples)
    worst = np.empty(n_samples, dtype=float)
    for r in range(n_samples):
        clear = cloud[r]
        worst[r] = _worst_band_revisit(
            band_acquisition_days, passes_needed_per_band, clear, n_days
        )
    return float(np.percentile(worst, percentile)), float(worst.mean())


def _worst_band_revisit(
    bands: list[np.ndarray],
    passes_needed: list[int],
    clear_days: np.ndarray,
    n_days: int,
) -> float:
    """Worst (over bands) revisit for one Monte-Carlo cloud realisation."""
    worst = 0.0
    for acqui_days, k in zip(bands, passes_needed):
        if acqui_days.size == 0 or k <= 0:
            return float(n_days)
        clear_acqui = acqui_days[clear_days[acqui_days]]
        n_clear = int(clear_acqui.size)
        if n_clear == 0:
            return float(n_days)
        revisit = k * n_days / max(1, n_clear)
        if revisit > worst:
            worst = revisit
    return float(worst)


def revisit_tiling_95p(
    passes_needed_per_band: list[int],
    reached_per_band: list[bool],
    usable_passes_per_day: float,
    n_days: int,
    p: float,
    rho: float,
    rng: np.random.Generator,
    n_samples: int = 200,
    percentile: float = 95.0,
) -> tuple[float, float]:
    """Orbit-level tiling revisit under persistent clouds (95th percentile).

    The revisit (days) is computed from the single worst-band aggregate:
    the maximum ``passes_needed`` across all latitude bands divided by the
    clear-day pass rate. The clear-pass realisation is driven by a
    Markov-persistent daily cloud chain: a day either yields its full usable
    passes (clear) or none (cloudy). The reported revisit is the worst-band
    aggregate at the requested percentile across Monte-Carlo runs; bands
    flagged unreachable (``reached_per_band`` False) force the horizon as
    the revisit.

    This keeps the MATLAB ``_US_WIDTH / eff_swath`` tiling count but derives
    the cadence from the real Orekit period and replaces ``1/(1-cloud)``
    with the Markov Monte-Carlo at a configurable confidence level.
    """
    if not passes_needed_per_band or usable_passes_per_day <= 0:
        return float("nan"), float("nan")
    cloud = sample_cloud_days(n_days, p, rho, rng, n_samples=n_samples)
    n_clear = cloud.sum(axis=1)  # (n_samples,)
    # clear-day fraction per run; guard against zero
    frac = n_clear / float(n_days)
    frac = np.where(frac <= 0, 1e-6, frac)
    clear_passes_per_day = usable_passes_per_day * frac  # (n_samples,)
    # revisit per run = max over bands of passes_needed_band / clear_passes_per_day
    band_rev = np.array(passes_needed_per_band, dtype=float)  # (n_bands,)
    reached = np.array(reached_per_band, dtype=bool)
    band_rev = np.where(reached, band_rev, float(n_days))
    # (n_samples,) = max_band(band_rev) / clear_passes_per_day
    worst_band_needed = float(band_rev.max())
    revisit = worst_band_needed / clear_passes_per_day  # (n_samples,)
    # any unreachable band ⇒ revisit is horizon
    if (~reached).any():
        revisit = np.full_like(revisit, float(n_days))
    lo = 100.0 - percentile
    return float(np.percentile(revisit, percentile)), float(revisit.mean())