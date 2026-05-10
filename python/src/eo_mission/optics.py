"""SNR and MTF calculations for EO satellite optical systems.

Port of SNRfunction.m and MTFfunction.m. Formulas are identical to the MATLAB
originals; only the I/O structure differs (returns DataFrames instead of
writing files directly).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Sequence


# ── Physical constants ─────────────────────────────────────────────────────────
_H_PLANCK    = 6.626e-34   # J·s
_C_LIGHT     = 3.00e8      # m/s
_R_EARTH     = 6_371e3     # m
_MU_EARTH    = 3.986e14    # m³/s²  (gravitational parameter)
_BANDWIDTH_M = 20e-3       # m  (20 nm spectral bandwidth)
_TDI         = 1           # TDI stages
_RAD_REF     = 100         # W/(m²·sr·μm)  reference radiance

# ── Noise floor (electrons RMS) ────────────────────────────────────────────────
_NOISE = {
    "dark":    50,   # thermal dark current
    "read":   100,   # readout electronics
    "preamp":   5,   # preamplifier
    "video":   10,   # video line
    "jitter":   5,   # platform jitter
    "emc":      5,   # electromagnetic interference
    "quant":    2,   # ADC quantisation
    "nonlin":   2,   # non-linearity
}
_N_FIXED_SQ = sum(v**2 for v in _NOISE.values())


def compute_snr(
    lambda_c: float,
    pixel_size: float,
    eta: float,
    tau: float,
    gsd: float,
    r_obs: float,
    altitudes: Sequence[float],
    diameters: Sequence[float],
    snr_req: float = 400.0,
) -> pd.DataFrame:
    """Compute SNR heatmap over orbital altitude × pupil diameter.

    Parameters
    ----------
    lambda_c : float
        Central wavelength (m).
    pixel_size : float
        Detector pixel pitch (m).
    eta : float
        Quantum efficiency.
    tau : float
        Telescope transmission.
    gsd : float
        Ground sampling distance (m).
    r_obs : float
        Secondary mirror linear obstruction ratio.
    altitudes : sequence of float
        Orbital altitudes (km).
    diameters : sequence of float
        Pupil diameters (mm).
    snr_req : float
        Minimum SNR threshold; values below become NaN.

    Returns
    -------
    pd.DataFrame
        SNR values indexed by altitude (km), columns by diameter (mm).
        Entries below ``snr_req`` are NaN.
    """
    h_arr = np.asarray(altitudes, dtype=float)
    d_arr = np.asarray(diameters, dtype=float)

    # Orbital velocity and integration time per altitude (1-D vectors)
    h_m       = h_arr * 1e3
    focal_len = (pixel_size * h_m) / gsd                    # m
    v_orb     = np.sqrt(_MU_EARTH / (_R_EARTH + h_m))       # m/s
    t_int     = gsd / v_orb                                  # s

    # Effective diameter accounting for secondary obscuration (scalar per column)
    D_m       = d_arr / 1e3                                  # mm → m
    D_eff     = np.sqrt(D_m**2 * (1 - r_obs**2))            # m

    # Irradiance at detector: broadcast (n_alt,) × (n_diam,)
    f_over_D  = focal_len[:, None] / D_eff[None, :]         # dimensionless
    irradiance = (np.pi * tau * _BANDWIDTH_M * _RAD_REF) / (1 + 4 * f_over_D**2)

    # Signal electrons
    pixel_area = pixel_size**2
    Ne = (irradiance * pixel_area * eta * lambda_c * _TDI * t_int[:, None]) / (_H_PLANCK * _C_LIGHT)

    # Total noise and SNR
    N_total = np.sqrt(_N_FIXED_SQ + Ne)
    snr     = Ne / N_total

    snr[snr < snr_req] = np.nan

    index   = pd.Index(h_arr, name="altitude_km")
    columns = pd.Index(d_arr, name="diameter_mm")
    return pd.DataFrame(snr, index=index, columns=columns)


def _otf_obscured(X: float, R: float) -> float:
    """Compute diffraction OTF for a circular aperture with central obscuration.

    Implements equations 6-10 to 6-15 of the mission analysis.

    Parameters
    ----------
    X : float
        Normalised spatial frequency (f / f_cutoff).
    R : float
        Secondary mirror linear obscuration ratio.
    """
    Y = X / R if R > 0 else 0.0

    # Term A
    A = (2 / np.pi) * (np.arccos(X) - X * np.sqrt(1 - X**2)) if 0 <= X <= 1 else 0.0

    # Term B
    B = ((2 * R**2) / np.pi) * (np.arccos(Y) - Y * np.sqrt(1 - Y**2)) if 0 <= Y <= 1 else 0.0

    # alpha for term C
    arg = (1 + R**2 - 4 * X**2) / (2 * R) if R > 0 else 0.0
    alpha = np.arccos(np.clip(arg, -1, 1)) if R > 0 else 0.0

    # Term C
    if R == 0 or X <= 0:
        C = 0.0
    elif X <= (1 - R) / 2:
        C = -2 * R**2
    elif X < (1 + R) / 2:
        C = (
            (2 * R / np.pi) * np.sin(alpha)
            + ((1 + R**2) / np.pi) * alpha
            - ((2 * (1 - R**2)) / np.pi) * np.arctan((1 + R) / (1 - R) * np.tan(alpha / 2))
            - 2 * R**2
        )
    else:
        C = 0.0

    otf = (A + B + C) / (1 - R**2)
    return max(otf, 0.0) if np.isfinite(otf) else 0.0


def compute_mtf(
    lambda_: float,
    pixel_size: float,
    mtf_detector: float,
    mtf_alignment: float,
    gsd: float,
    r_obs: float,
    altitudes: Sequence[float],
    diameters: Sequence[float],
    mtf_req: float = 0.25,
) -> pd.DataFrame:
    """Compute system MTF heatmap over orbital altitude × pupil diameter.

    Parameters
    ----------
    lambda_ : float
        Wavelength (m).
    pixel_size : float
        Detector pixel pitch (m).
    mtf_detector : float
        Detector MTF (sampling/pixel-response factor).
    mtf_alignment : float
        Telescope alignment MTF.
    gsd : float
        Ground sampling distance (m).
    r_obs : float
        Secondary mirror linear obstruction ratio.
    altitudes : sequence of float
        Orbital altitudes (km).
    diameters : sequence of float
        Pupil diameters (mm).
    mtf_req : float
        Minimum MTF threshold; values below become NaN.

    Returns
    -------
    pd.DataFrame
        System MTF indexed by altitude (km), columns by diameter (mm).
        Entries below ``mtf_req`` are NaN.
    """
    # Static MTF factors (combined once)
    mtf_aberrations   = 0.95
    mtf_manufacturing = 0.98
    mtf_vibrations    = 0.99
    mtf_thermoelastic = 0.95
    mtf_margin        = 0.90
    mtf_static = (
        mtf_margin * mtf_vibrations * mtf_manufacturing
        * mtf_thermoelastic * mtf_aberrations * mtf_detector
    )

    h_arr = np.asarray(altitudes, dtype=float)
    d_arr = np.asarray(diameters, dtype=float)
    result = np.full((len(h_arr), len(d_arr)), np.nan)

    f_nyquist = 1 / (2 * pixel_size)   # cy/m

    for i, h_km in enumerate(h_arr):
        h_m       = h_km * 1e3
        focal_len = (h_m * pixel_size) / gsd        # m
        f_x_raw   = f_nyquist * focal_len            # cy/m, at focal plane

        for j, d_mm in enumerate(d_arr):
            d_m  = d_mm / 1e3
            f_co = d_m / lambda_                     # cy/m
            f_x  = min(f_x_raw, f_co)
            X    = f_x / f_co if f_co > 0 else 0.0

            otf   = _otf_obscured(X, r_obs)
            mtf_v = otf * mtf_alignment * mtf_static

            result[i, j] = mtf_v if mtf_v >= mtf_req else np.nan

    index   = pd.Index(h_arr, name="altitude_km")
    columns = pd.Index(d_arr, name="diameter_mm")
    return pd.DataFrame(result, index=index, columns=columns)


def plot_heatmap(
    df: pd.DataFrame,
    title: str,
    colorbar_label: str,
    output_path: str | Path,
    vmin: float | None = None,
    vmax: float | None = None,
) -> None:
    """Save a heatmap of a SNR or MTF DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Data with altitude (km) as index and diameter (mm) as columns.
    title : str
        Plot title (supports LaTeX via matplotlib).
    colorbar_label : str
        Colorbar axis label.
    output_path : str or Path
        Destination PNG file path.
    vmin, vmax : float, optional
        Colorscale limits; determined from data if not given.
    """
    fig, ax = plt.subplots(figsize=(10, 7))
    alts = df.index.to_numpy()
    diams = df.columns.to_numpy()

    img = ax.imshow(
        df.to_numpy(),
        aspect="auto",
        origin="lower",
        extent=[diams[0], diams[-1], alts[0], alts[-1]],
        cmap="viridis",
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_xlabel("Pupil diameter (mm)", fontsize=12)
    ax.set_ylabel("Orbital altitude (km)", fontsize=12)
    ax.set_title(title, fontsize=14)
    cb = fig.colorbar(img, ax=ax)
    cb.set_label(colorbar_label, fontsize=12)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
