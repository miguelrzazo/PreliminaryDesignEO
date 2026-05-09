"""Run the EO mission optical analysis (equivalent to master_ver1.m, ver 1.0 scenario).

Usage:
    cd python/
    pip install -e .
    python scripts/run_analysis.py [--output-dir OUTPUT]
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))
from eo_mission.optics import compute_snr, compute_mtf, plot_heatmap

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
    {"name": "Refractivo", "mtf_align": 0.90, "tau": 0.80, "r_obs": 0.0},
    {"name": "Korsch",     "mtf_align": 0.85, "tau": 0.65, "r_obs": 0.2},
    {"name": "Cassegrain", "mtf_align": 0.80, "tau": 0.70, "r_obs": 0.2},
    {"name": "TMA",        "mtf_align": 0.70, "tau": 0.60, "r_obs": 0.0},
]

GSD          = 50.0
SNR_REQ      = 400.0
MTF_REQ      = 0.25
ALTITUDES    = np.arange(400, 1401, 10, dtype=float)
DIAMETERS    = np.arange(20,   601,  1, dtype=float)

BANDS = [
    ("Lambda1", LAMBDA["Lambda1"], DETECTORS_12, 0),
    ("Lambda2", LAMBDA["Lambda2"], DETECTORS_12, 0),
    ("Lambda3", LAMBDA["Lambda3"], DETECTORS_3,  3),
]


def main():
    parser = argparse.ArgumentParser(description="EO mission optical analysis")
    parser.add_argument("--output-dir", default="output", type=Path)
    args = parser.parse_args()

    out = args.output_dir
    (out / "SNR").mkdir(parents=True, exist_ok=True)
    (out / "MTF").mkdir(parents=True, exist_ok=True)

    for band_tag, lam, detectors, id_offset in BANDS:
        for d_i, det in enumerate(detectors):
            det_id = d_i + 1 + id_offset
            for t_i, tel in enumerate(TELESCOPES):
                tel_id = t_i + 1
                tag = f"{band_tag}_Detector{det_id}_Telescopio{tel_id}"
                print(f"  Computing SNR  {tag} ...")
                snr_df = compute_snr(
                    lambda_c   = lam,
                    pixel_size = det["pixel_size"],
                    eta        = det["eta"],
                    tau        = tel["tau"],
                    gsd        = GSD,
                    r_obs      = tel["r_obs"],
                    altitudes  = ALTITUDES,
                    diameters  = DIAMETERS,
                    snr_req    = SNR_REQ,
                )
                snr_df.to_csv(out / "SNR" / f"SNR_{tag}_resultados.csv")
                plot_heatmap(
                    snr_df,
                    title=f"SNR — {tel['name']} — Detector {det_id} (λ={lam*1e6:.2f} μm)",
                    colorbar_label="SNR",
                    output_path=out / "SNR" / f"SNR_{tag}_heatmap.png",
                    vmin=SNR_REQ,
                    vmax=2000,
                )

                print(f"  Computing MTF  {tag} ...")
                mtf_df = compute_mtf(
                    lambda_       = lam,
                    pixel_size    = det["pixel_size"],
                    mtf_detector  = det["mtf_det"],
                    mtf_alignment = tel["mtf_align"],
                    gsd           = GSD,
                    r_obs         = tel["r_obs"],
                    altitudes     = ALTITUDES,
                    diameters     = DIAMETERS,
                    mtf_req       = MTF_REQ,
                )
                mtf_df.to_csv(out / "MTF" / f"MTF_{tag}_resultados.csv")
                plot_heatmap(
                    mtf_df,
                    title=f"MTF — {tel['name']} — Detector {det_id} (λ={lam*1e6:.2f} μm)",
                    colorbar_label="MTF",
                    output_path=out / "MTF" / f"MTF_{tag}_heatmap.png",
                    vmin=MTF_REQ,
                    vmax=0.35,
                )

    print(f"\nDone. Results in: {out.resolve()}")


if __name__ == "__main__":
    main()
