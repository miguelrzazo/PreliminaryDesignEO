"""Run the complete EO mission analysis pipeline.

Equivalent to master_ver1.m (ver 1.0 scenario) + electrical budget.

Usage:
    cd python/
    pip install -e .
    python scripts/run_analysis.py [--output-dir OUTPUT] [--skip-mass] [--skip-electrical]
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))
from eo_mission.optics import compute_snr, compute_mtf, plot_heatmap
from eo_mission.mass import compute_dry_mass, compute_total_mass
from eo_mission.coverage import compute_coverage
from eo_mission.electrical import compute_power_budget
from eo_mission.analysis import cross_data_analysis, optimum_configuration

# ── Scenario: Analytical model ver 1.0 ────────────────────────────────────────
LAMBDA  = {"Lambda1": 1.61e-6, "Lambda2": 2.01e-6, "Lambda3": 0.76e-6}

# CO2-band detectors (det 1-3)
DETECTORS_12 = [
    {"name": "CAPYORK",     "eta": 0.60, "mtf_det": 0.45, "n_pix": 1024, "pixel_size": 15e-6},
    {"name": "CHROMA_D",    "eta": 0.80, "mtf_det": 0.40, "n_pix": 3000, "pixel_size": 18e-6},
    {"name": "SATURN_VISIR","eta": 0.60, "mtf_det": 0.45, "n_pix": 1000, "pixel_size": 30e-6},
]

# O2 A-band detectors (det 4-6)
DETECTORS_3 = [
    {"name": "SATURN_VISIR","eta": 0.60, "mtf_det": 0.50, "n_pix": 1000, "pixel_size": 30e-6},
    {"name": "CHROMA_D",    "eta": 0.80, "mtf_det": 0.40, "n_pix": 3000, "pixel_size": 18e-6},
    {"name": "CMOS",        "eta": 0.35, "mtf_det": 0.36, "n_pix":  512, "pixel_size": 25e-6},
]

TELESCOPES = [
    {"name": "Refractivo", "mtf_align": 0.90, "tau": 0.80, "r_obs": 0.0, "fov_deg": 4.0},
    {"name": "Korsch",     "mtf_align": 0.85, "tau": 0.65, "r_obs": 0.2, "fov_deg": 6.0},
    {"name": "Cassegrain", "mtf_align": 0.80, "tau": 0.70, "r_obs": 0.2, "fov_deg": 5.0},
    {"name": "TMA",        "mtf_align": 0.70, "tau": 0.60, "r_obs": 0.0, "fov_deg": 8.0},
]

GSD       = 50.0
SNR_REQ   = 400.0
MTF_REQ   = 0.25
ALTITUDES = np.arange(400, 1401, 10, dtype=float)
DIAMETERS = np.arange(20,   601,  1, dtype=float)
SWATHS_KM = np.arange(10,   500,  5, dtype=float)

BANDS = [
    ("Lambda1", LAMBDA["Lambda1"], DETECTORS_12, 0),
    ("Lambda2", LAMBDA["Lambda2"], DETECTORS_12, 0),
    ("Lambda3", LAMBDA["Lambda3"], DETECTORS_3,  3),
]

CONFIGS = [(1, 1), (1, 2), (2, 1), (2, 2)]


def main():
    parser = argparse.ArgumentParser(description="EO mission analysis pipeline")
    parser.add_argument("--output-dir", default="output", type=Path)
    parser.add_argument("--skip-mass",       action="store_true")
    parser.add_argument("--skip-electrical", action="store_true")
    args = parser.parse_args()

    out = args.output_dir
    for sub in ("SNR", "MTF", "Coverage", "Mass", "Electrical", "HvsDmin", "OptimumConfigs"):
        (out / sub).mkdir(parents=True, exist_ok=True)

    # ── Step 1: Optics (SNR + MTF) ────────────────────────────────────────────
    print("\n=== Step 1: SNR / MTF ===")
    snr_store: dict = {}
    mtf_store: dict = {}

    for band_tag, lam, detectors, id_offset in BANDS:
        for d_i, det in enumerate(detectors):
            det_id = d_i + 1 + id_offset
            for t_i, tel in enumerate(TELESCOPES):
                tel_id = t_i + 1
                tag    = f"{band_tag}_Detector{det_id}_Telescopio{tel_id}"

                print(f"  SNR  {tag}")
                snr_df = compute_snr(
                    lambda_c=lam, pixel_size=det["pixel_size"], eta=det["eta"],
                    tau=tel["tau"], gsd=GSD, r_obs=tel["r_obs"],
                    altitudes=ALTITUDES, diameters=DIAMETERS, snr_req=SNR_REQ,
                )
                snr_df.to_csv(out / "SNR" / f"SNR_{tag}_resultados.csv")
                plot_heatmap(snr_df,
                    title=f"SNR — {tel['name']} — Det {det_id} (λ={lam*1e6:.2f} μm)",
                    colorbar_label="SNR",
                    output_path=out / "SNR" / f"SNR_{tag}_heatmap.png",
                    vmin=SNR_REQ, vmax=2000)
                snr_store[tag] = snr_df

                print(f"  MTF  {tag}")
                mtf_df = compute_mtf(
                    lambda_=lam, pixel_size=det["pixel_size"],
                    mtf_detector=det["mtf_det"], mtf_alignment=tel["mtf_align"],
                    gsd=GSD, r_obs=tel["r_obs"],
                    altitudes=ALTITUDES, diameters=DIAMETERS, mtf_req=MTF_REQ,
                )
                mtf_df.to_csv(out / "MTF" / f"MTF_{tag}_resultados.csv")
                plot_heatmap(mtf_df,
                    title=f"MTF — {tel['name']} — Det {det_id} (λ={lam*1e6:.2f} μm)",
                    colorbar_label="MTF",
                    output_path=out / "MTF" / f"MTF_{tag}_heatmap.png",
                    vmin=MTF_REQ, vmax=0.35)
                mtf_store[tag] = mtf_df

    # ── Step 2: Coverage ──────────────────────────────────────────────────────
    print("\n=== Step 2: Coverage ===")
    cov_store: dict = {}
    all_dets = [(d_i, det, "12") for d_i, det in enumerate(DETECTORS_12)] + \
               [(d_i, det, "3")  for d_i, det in enumerate(DETECTORS_3)]

    for n_sat, n_tel in CONFIGS:
        for t_i, tel in enumerate(TELESCOPES):
            for d_i, det, band_group in all_dets:
                det_id = d_i + 1 + (3 if band_group == "3" else 0)
                tag    = f"{n_sat}sat_{n_tel}tel_{tel['name']}_Det{det_id}"
                print(f"  Coverage {tag}")
                cov_df = compute_coverage(
                    altitudes_km=ALTITUDES, swaths_km=SWATHS_KM,
                    gsd=GSD, n_pix=det["n_pix"], fov_limit_deg=tel["fov_deg"],
                    n_sat=n_sat, n_telescopes=n_tel, detector_type=d_i + 1,
                )
                cov_df.to_csv(out / "Coverage" / f"coverage_{tag}.csv")
                cov_store[tag] = cov_df

    # ── Step 3: Mass ─────────────────────────────────────────────────────────
    if not args.skip_mass:
        print("\n=== Step 3: Mass ===")
        dry_df = compute_dry_mass(DIAMETERS)
        dry_df.to_csv(out / "Mass" / "dry_mass.csv", index=False)

        # Cross-section approximation: circular aperture area × 1.1 margin
        A_cross = np.pi * (DIAMETERS / 2 / 1000) ** 2 * 1.1

        for n_sat, n_tel in CONFIGS:
            tag = f"{n_sat}sat_{n_tel}tel"
            dry_vals = dry_df["Masa_seca"].to_numpy() * (1.5 if n_tel == 2 else 1.0)
            # Sample total mass at representative altitudes (every 50 km)
            alt_sample = np.arange(400, 1401, 50, dtype=float)
            d_sample   = DIAMETERS[::10]
            dry_sample = dry_vals[::10]
            A_sample   = A_cross[::10]
            print(f"  Total mass {tag} ({len(alt_sample)} altitudes × {len(d_sample)} diameters)...")
            for j, (d_val, m_dry, A_val) in enumerate(zip(d_sample, dry_sample, A_sample)):
                tot_df = compute_total_mass(alt_sample, np.full(len(alt_sample), m_dry),
                                            np.full(len(alt_sample), A_val))
                tot_df["Diametro_pupila"] = d_val
                tot_df.to_csv(
                    out / "Mass" / f"MasaTotal_{tag}_D{d_val:.0f}mm.csv", index=False)
    else:
        print("\n[skip] Mass module")

    # ── Step 4: Electrical ────────────────────────────────────────────────────
    if not args.skip_electrical:
        print("\n=== Step 4: Electrical power budget ===")
        for alt in [500, 600, 700, 760, 800]:
            result = compute_power_budget(altitude_km=float(alt))
            print(f"  h={alt} km | panel_min={result['min_panel_area_m2']:.4f} m² "
                  f"| bat_min={result['min_battery_cap_wh']:.2f} Wh "
                  f"| eclipse={result['eclipse_fraction']:.3f}")
        import json
        with open(out / "Electrical" / "power_budget_summary.json", "w") as f:
            json.dump({
                str(alt): compute_power_budget(altitude_km=float(alt))
                for alt in [500, 600, 700, 760, 800]
            }, f, indent=2)
    else:
        print("\n[skip] Electrical module")

    print(f"\nDone. Results in: {out.resolve()}")


if __name__ == "__main__":
    main()
