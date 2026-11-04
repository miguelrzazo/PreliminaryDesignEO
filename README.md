# Preliminary Design of an Earth Observation Satellite

Bachelor's thesis (TFG) in Aerospace Engineering: preliminary design of a small-satellite constellation for CO₂/O₂ monitoring over the continental United States, covering payload sizing, orbit and mission analysis, launcher selection, and ground segment.

Read the thesis: [TFG_Miguel_Rosa.pdf](TFG_Miguel_Rosa.pdf) (Spanish).

## Result in brief

The selected design uses two satellites at 690 km in a Sun-synchronous frozen orbit, each with two refractive telescopes (38 mm pupil) and two H2RG detectors. At 80 m GSD, the sampled regional model gives a maximum effective revisit of 6.642 days with the assumed mean cloud correction. The estimated mass is 3.538 kg per satellite, including maintenance and deorbit reserves. Downlink sizing uses a 10 Mbps S-band reference for ASF AS1 equipment at Fairbanks; service availability is assumed.

## Repository layout

| Path | Contents |
|------|----------|
| `Latex_Code/` | Thesis LaTeX sources and figure sources |
| `MATLAB_Code/` | The design pipeline: `master.m` + `lib/` (authoritative) |
| `MATLAB_Code/legacy/` | Archived earlier model versions, kept for traceability only |
| `Compiled data/` | Regenerated MATLAB pipeline, selected mission budgets, and auxiliary GSD sensitivity data |
| `python/` | Python port of the analysis pipeline |
| `Defence/` | Final defence material (simulation videos) |

## Running the MATLAB pipeline

Run the optical and coverage sweep, selected mission budgets, and figure exports in MATLAB with Aerospace Toolbox:

```bash
matlab -batch "addpath('MATLAB_Code'); runThesisAnalysis"
```

`master.m` contains the editable design grid and requirements. `TFG_GSD` overrides the ground sample distance in metres. Results are written to `Compiled data/Pipeline` and `Compiled data/Mission`; the independent study scripts reuse `Mission/selected_design.csv`. Pipeline mass rankings include dry mass and maintenance propellant; `Mission/EndOfLife` adds the deorbit reserve.

The revisit computation uses the RevisitTime toolkit (Crisp & Livadiotti, University of Manchester), which must be on the MATLAB path.

## Python port

`python/` contains a NumPy/SciPy port of the MATLAB pipeline (optics, coverage, mass, and electrical modules) with tests. It requires Orekit (conda-forge only — not on PyPI):

```bash
mamba create -y -n eo -c conda-forge python=3.11 orekit openjdk jpype1 numpy scipy pandas matplotlib pyyaml pytest pytest-cov
mamba activate eo
cd python
pip install -e ".[dev]"
# Download orekit-data once (the adapter registers it for the JVM):
python -c "from orekit.pyhelpers import download_orekit_data_curdir; download_orekit_data_curdir()"
OREKIT_DATA_PATH="$PWD/orekit-data.zip" pytest
OREKIT_DATA_PATH="$PWD/orekit-data.zip" python scripts/run_analysis.py --orekit
```

This auxiliary implementation uses a different coverage and cloud model and does not supply the thesis results. The configuration is a single YAML file (`python/config/portfolio.yaml`) that defines the design space — bands, detectors, telescopes, constellations, altitudes, swaths, GSDs — and each module sweeps over it, returning a `LabeledTensor` (a self-describing `np.ndarray` with named axes) so the 6-D result grid stays easy to slice.

### What this is

A Python re-implementation of the MATLAB optical/coverage/mass pipeline behind the thesis. The physics (SNR, MTF, orbital geometry, SSO inclination, dry-mass empirical law, electrical budget) is faithful to `MATLAB_Code/lib/*.m`; the structure is portfolio-first and fully vectorised.

### How the coverage module uses Orekit

The coverage module (`python/src/eo_mission/coverage.py`) is the interesting one. The MATLAB version (`CoverageRevisitCalc.m`) wraps the external RevisitTime toolkit (Crisp & Livadiotti, University of Manchester) and divides by `1/(1−cloud)` as a mean-value cloud correction. The Python port replaces both pieces:

- **SSO inclination from J2.** The sun-synchronous inclination is solved from the J2 sun-sync condition (`cos_i = −SSO_nodal_rate·(a/R_E)^(7/2)/(1.5·n·J2)`, the same equation as `CoverageRevisitCalc.m`); altitudes where `|cos_i|>1` return NaN, exactly as the MATLAB code marks SSO-infeasible altitudes.
- **Orekit-propagated ground track.** A `KeplerianPropagator` initialised with an Orekit `KeplerianOrbit` (SSO inclination, J2000 epoch, EME2000 frame, EIGEN5C μ) is propagated over a configurable horizon (30 days default) and the sub-satellite ground track is sampled in the Veis-1950 Earth-fixed frame. The orbit period comes from Orebit (`getKeplerianPeriod`) instead of a hand-coded `2π·sqrt(r³/μ)`. Multi-satellite constellations are modelled as a plane train with equally spaced mean anomalies. Ground tracks are memoised to a JSON cache per (altitude, n-satellites) so portfolio sweeps reuse them.
- **CONUS latitude-band grid.** A coarse latitude/longitude grid spanning the continental US (25–49° N, a few longitude samples per band) provides the target set. The sensor's visible footprint radius is computed from the spherical geometry of the off-nadir cone (`γ_max` solved from the spherical triangle O-S-T). Revisit of a latitude band is the time to tile its longitudinal span: `passes_needed = ceil(band_width·cos(lat) / eff_swath)`, the MATLAB `_US_WIDTH / eff_swath` tiling count.
- **Markov-persistent clouds, 95th-percentile revisit.** The MATLAB `1/(1−cloud)` mean heuristic is replaced by a persistent two-state Markov chain (`python/src/eo_mission/clouds.py`) parameterised by steady-state cloud cover `p` and a persistence coefficient `rho` (lag-1 autocorrelation of the daily cloud state), which keeps the stationary distribution exactly P(cloudy)=p for any `rho`. Monte-Carlo sampling (200 runs default, seeded) folds daily cloud state onto the acquisition cadence; the reported revisit is the worst-latitude-band value at the 95th percentile across runs (the mean is also tracked). This reports the revisit you actually achieve at 95% confidence rather than a mean correction.

The pure-analytical single-target formula used by an earlier version of the port — which had invented a `1/0.5` daylight penalty and a `×0.35` USA-grid fudge factor absent from the MATLAB physics — has been removed. Coverage is Orekit-only; there is no JVM-free fallback. The pipeline (`run_pipeline`) initializes the JVM by default.

### Other improvements over the MATLAB original

- **SNR calibration margin.** The 0.9 factor of `SNRfunction.m` (`SNR_value = 0.9·Ne/N_total`) is now applied in `optics.py`, matching the MATLAB reference exactly (the previous port omitted it).
- **Vectorisation.** The MATLAB scripts loop cell-by-cell and write results to disk. The Python port keeps the same physics but operates on whole `np.ndarray`s, so SNR/MTF/mass are evaluated as broadcasts over the full grid in one go. The `LabeledTensor` (`python/src/eo_mission/tensor.py`) keeps the axis names and coordinates alongside the array, so downstream code can do `snr.values[g_i, b_i, d_i, t_i, h_i, dia_i]` without losing track of what each index means.
- **Portfolio sweep.** The original MATLAB pipeline picks a single configuration and sizes it. The Python port is portfolio-first: `rank_feasible_solutions` (`python/src/eo_mission/analysis.py:159`) enumerates *every* (band, detector, telescope, constellation, altitude, diameter, swath) cell that satisfies SNR, MTF, coverage and aperture thresholds, then sorts the feasible set by dry mass → revisit → altitude, and `best_solution_per_gsd` keeps the lightest one per GSD. The result is `tables/ranked_solutions.csv` and `tables/best_per_gsd.csv` — a ranked design catalogue rather than a single point.

### Tests

`pytest` runs pure-Python unit tests (Markov cloud sampler, SSO inclination feasibility, SNR margin, MTF, mass, electrical, analysis) plus Orekit integration tests that propagate a real SSO orbit and check the revisit tensor. The Orekit tests need `OREKIT_DATA_PATH` set to the `orekit-data.zip` (or directory) registered with `setup_orekit_curdir`.

## Building the thesis

Requires a TeX distribution with `latexmk` and `biber`:

```bash
cd Latex_Code
latexmk -pdf 00main.tex
```

## License

Code is MIT-licensed; the thesis text and figures are all rights reserved. See [LICENSE.md](LICENSE.md).
