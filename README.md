# Preliminary Design of an Earth Observation Satellite

Bachelor's thesis (TFG) in Aerospace Engineering: preliminary design of a small-satellite constellation for CO₂/O₂ monitoring over the continental United States, covering payload sizing, orbit and mission analysis, launcher selection, and ground segment.

**📄 Read the full thesis: [`Latex_Code/00main.pdf`](Latex_Code/00main.pdf)** (Spanish)

## Result in brief

A constellation of two satellites in a 520 km Sun-synchronous orbit, each carrying two refractive telescopes with Teledyne H2RG detectors (28 mm pupil diameter), achieves 80 m GSD imaging of the continental US with weekly revisit. The design closes with a total constellation mass of a few kilograms of payload optics per satellite plus the propellant to maintain the orbit for eight years, downlinking through the NOAA Fairbanks (Alaska) X-band station.

## Repository layout

| Path | Contents |
|------|----------|
| `Latex_Code/` | Thesis LaTeX sources and the compiled PDF |
| `MATLAB_Code/` | The design pipeline: `master.m` + `lib/` (authoritative) |
| `MATLAB_Code/legacy/` | Archived earlier model versions, kept for traceability only |
| `Compiled data/` | Final simulation outputs (MTF, SNR, coverage, h–D_min, mass) and the GSD sensitivity summary |
| `python/` | Python port of the analysis pipeline |
| `Defence/` | Final defence material (simulation videos) |

## Running the MATLAB pipeline

The full sizing pipeline (optical MTF/SNR heatmaps, coverage/revisit simulation, configuration cross-check, dry and total mass) runs from `master.m`:

```bash
TFG_GSD=80 TFG_OUTPUT_ROOT="Compiled data" \
  matlab -batch "cd('MATLAB_Code'); master"
```

`TFG_GSD` sets the required ground sample distance in metres, so the whole study can be repeated for a different resolution requirement. `TFG_OUTPUT_ROOT` sets where results are written.

The revisit computation uses the RevisitTime toolkit (Crisp & Livadiotti, University of Manchester), which must be on the MATLAB path.

## Python port

`python/` contains a NumPy/SciPy port of the MATLAB pipeline (optics, coverage, mass, and electrical modules) with tests:

```bash
cd python
pip install -e ".[dev]"
pytest
python scripts/run_analysis.py
```

It was vibecoded after the thesis, for fun, to re-validate the sizing study in a friendlier stack (frozen dataclasses, NumPy, PyYAML, matplotlib) and to play with packaging. The configuration is a single YAML file (`python/config/portfolio.yaml`) that defines the design space — bands, detectors, telescopes, constellations, altitudes, swaths, GSDs — and each module sweeps over it, returning a `LabeledTensor` (a self-describing `np.ndarray` with named axes) so the 6-D result grid stays easy to slice.

A couple of small improvements over the MATLAB original snuck in:

- **Vectorisation.** The MATLAB scripts loop cell-by-cell and write results to disk. The Python port keeps the same physics but operates on whole `np.ndarray`s, so SNR/MTF/coverage/mass are evaluated as broadcasts over the full grid in one go. The `LabeledTensor` (`python/src/eo_mission/tensor.py`) keeps the axis names and coordinates alongside the array, so downstream code can do `snr.values[g_i, b_i, d_i, t_i, h_i, dia_i]` without losing track of what each index means.
- **Portfolio sweep.** The original MATLAB pipeline picks a single configuration and sizes it. The Python port is portfolio-first: `rank_feasible_solutions` (`python/src/eo_mission/analysis.py:159`) enumerates *every* (band, detector, telescope, constellation, altitude, diameter, swath) cell that satisfies SNR, MTF, coverage and aperture thresholds, then sorts the feasible set by dry mass → revisit → altitude, and `best_solution_per_gsd` keeps the lightest one per GSD. The result is `tables/ranked_solutions.csv` and `tables/best_per_gsd.csv` — a ranked design catalogue rather than a single point.

The coverage module is the interesting one: the MATLAB version (`CoverageSSOAnaliticalfunction.m`) is a closed-form expression for revisit time. Here, the same computation is wrapped behind `OrekitCoverageAdapter` (`python/src/eo_mission/coverage.py:103`), which lazily starts a JVM via JPype, points it at an `orekit-data` directory, and pulls the orbital period from an `org.orekit.orbits.KeplerianOrbit` (SSO inclination, J2000 epoch, EME2000 frame) instead of a hand-coded `2π·sqrt(r³/μ)`. The pure-analytical path stays as a fallback for unit tests and for environments where the JVM/orekit-data aren't installed, so the rest of the pipeline is unaware of Orekit. Revisit values are also memoised to a small JSON cache so grid sweeps don't recompute identical cells. A small smoke test (`python/tests/test_orekit_integration.py`) only runs inside the Docker Orekit runtime (`EO_MISSION_RUN_OREKIT_TESTS=1`).

## Building the thesis

Requires a TeX distribution with `latexmk` and `biber`:

```bash
cd Latex_Code
latexmk -pdf 00main.tex
```

## License

Code is MIT-licensed; the thesis text and figures are all rights reserved. See [LICENSE.md](LICENSE.md).
