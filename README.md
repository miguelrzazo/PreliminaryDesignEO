# Preliminary Design of an Earth Observation Satellite

An aerospace engineering bachelor's thesis on the preliminary design of a small-satellite constellation for CO₂/O₂ observation over the continental United States. The project connects optical performance and regional coverage requirements to orbit selection, spacecraft mass, propulsion, communications and onboard storage.

**[Read the latest built thesis (PDF, in Spanish)](TFG_Miguel_Rosa.pdf)** — 210 pages, built on 9 November 2026.

MATLAB supplies the thesis calculations. A later Python/Orekit extension recomputes the design trade space and explores the resulting mission through propagated orbits and interactive visualizations. Its orbital coverage model is evaluated separately from the thesis model.

## Selected design

| Parameter | Thesis result |
| --- | --- |
| Constellation | Two satellites, phased by 180° |
| Orbit | Sun-synchronous frozen orbit at 690 km; 06:00 LTAN |
| Payload per satellite | Two refractive telescopes, 38 mm pupil, two H2RG detectors |
| Ground sampling distance | 80 m |
| Acquired swath | 240 km |
| Maximum effective regional revisit | 6.642 days on the sampled latitude grid |
| Mass per satellite | 3.538 kg, including maintenance and deorbit reserves |
| Minimum downlink rate | Approximately 5.024 Mbps |
| Maximum required onboard storage | Approximately 1.327 GB |

These are preliminary design estimates. The revisit result uses the thesis's mean cloud correction, rather than measured weather sequences. The 80 m requirement is a spatial design assumption, not a demonstration of atmospheric CO₂ retrieval performance. The 10 Mbps S-band reference is documented for ASF AS1 equipment; the simulated contacts use NOAA/Fairbanks coordinates, and station service availability remains an assumption.

## Engineering implementation

The MATLAB analysis combines SNR and MTF constraints with detector geometry, telescope field of view, regional revisit, dry mass and maintenance propellant. The selected solution then feeds the frozen-orbit, contact, memory and end-of-life studies. Pipeline selection mass includes maintenance propellant; the final mission mass also includes the deorbit reserve.

The implementation was subsequently reorganized and optimized:

- SNR and MTF calculations evaluate altitude/aperture grids with array operations, separating numerical kernels from plotting and file exports.
- Revisit calculations are cached by their orbital and viewing inputs so equivalent hardware cases reuse the same solution.
- Explicit configurations and result structures replace dependencies on caller workspace variables and intermediate files.
- Independent mission studies share the selected design, and regression tests compare numerical kernels with reference implementations.

Sequential propagation, numerical integration and memory evolution retain their time ordering. The repository includes an optical benchmark; no hardware-independent speedup is claimed.

The later Python extension is **vibecoded**. Its independent search reads the MATLAB input grid, recalculates optical and maintenance-inclusive mass constraints, and evaluates coverage with Orekit. It retains MATLAB’s mean cloud adjustment for selection and reports persistent Markov cloud scenarios separately. The adopted-design simulation is a validation case, rather than the optimization input. Its numerical assumptions, validation and runtime requirements are documented in the [Python guide](python/README.md). Python results do not replace the MATLAB results in the submitted thesis.

## Repository

| Path | Purpose |
| --- | --- |
| [MATLAB_Code/](MATLAB_Code/) | Authoritative thesis analysis, numerical kernels and tests |
| [Compiled data/](Compiled%20data/) | Saved pipeline results, selected mission budgets and figures |
| [Latex_Code/](Latex_Code/) | Thesis sources and figures |
| [python/](python/) | Subsequent Python/Orekit analysis and visualization |
| [Defence/](Defence/) | Supporting defence material |

`MATLAB_Code/legacy/` contains earlier model versions. Local bibliography and tutor revisions are excluded from version control.

## Reproducing the MATLAB analysis

Requires MATLAB, Aerospace Toolbox and the bundled RevisitTime toolkit by Crisp and Livadiotti (University of Manchester).

From the repository root:

```sh
matlab -batch "addpath('MATLAB_Code'); runThesisAnalysis"
matlab -batch "addpath('MATLAB_Code/tests'); results = runChecks();"
```

`MATLAB_Code/master.m` contains the design grid and requirements. The pipeline writes results to `Compiled data/Pipeline` and `Compiled data/Mission`; `Mission/selected_design.csv` records the adopted configuration. `TFG_GSD` overrides the ground sample distance in metres.

To compare the scalar and vectorized optical calculations on the same machine:

```sh
matlab -batch "addpath('MATLAB_Code/tests'); timings = benchmarkOpticalRefactor(); disp(timings)"
```

## Building the thesis

Requires a TeX distribution with `latexmk` and `biber`:

```sh
cd Latex_Code
TFG_LATEX_BUILD_DIR=/tmp/tfg-latex-build latexmk -pdf 00main.tex
```

Intermediate files are written outside the source tree. The compiled document is `Latex_Code/00main.pdf`. To update the linked deliverable after a successful build:

```sh
cp 00main.pdf ../TFG_Miguel_Rosa.pdf
```

The root [TFG_Miguel_Rosa.pdf](TFG_Miguel_Rosa.pdf) is the versioned deliverable.

## License

Code is MIT-licensed. Thesis text and figures are all rights reserved; see [LICENSE.md](LICENSE.md).
