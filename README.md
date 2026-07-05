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
| `Defence/` | Final defence material (simulation videos) |

## Running the MATLAB pipeline

The full sizing pipeline (optical MTF/SNR heatmaps, coverage/revisit simulation, configuration cross-check, dry and total mass) runs from `master.m`:

```bash
TFG_GSD=80 TFG_OUTPUT_ROOT="Compiled data" \
  matlab -batch "cd('MATLAB_Code'); master"
```

`TFG_GSD` sets the required ground sample distance in metres, so the whole study can be repeated for a different resolution requirement. `TFG_OUTPUT_ROOT` sets where results are written.

The revisit computation uses the RevisitTime toolkit (Crisp & Livadiotti, University of Manchester), which must be on the MATLAB path.

## Building the thesis

Requires a TeX distribution with `latexmk` and `biber`:

```bash
cd Latex_Code
latexmk -pdf 00main.tex
```

## License

Code is MIT-licensed; the thesis text and figures are all rights reserved. See [LICENSE.md](LICENSE.md).
