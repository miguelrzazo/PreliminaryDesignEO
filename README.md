# Preliminary Design of an Earth Observation Satellite — TFG

MATLAB-based analytical and numerical tool for satellite mission design and optimization.
Evaluates detector/telescope combinations across SNR, MTF, coverage revisit time, and mass budgets for Sun-Synchronous orbits.

## Repository Structure

```
MATLAB_Code/
├── Analitical model ver 1.0/       # First analytical model iteration
├── Analitical model with wrapper 1.0/  # Wrapper-based version
├── Analitical model with wrapper 1.5/  # Final analytical solution (ver 1.5)
├── Propagation model ver 2.0/      # Numerical propagation model (long runtime)
├── New Folder/                     # Orchestrated pipeline (master.m + refactored functions)
└── *.m                             # Standalone analysis scripts (MTF, SNR, coverage, electrical...)

Latex_Code/     — LaTeX thesis source
Compiled data/  — Reference simulation outputs (CSV/PNG)
Beamer_Slides/  — Presentation slides
Tikz_images/    — TikZ vector graphics for thesis
```

## Model Notes

The **final analytical solution** was obtained with `Analitical model with wrapper 1.5`. The analytical approach was chosen over the propagation model (ver 2.0) because the propagation model has excessively long runtimes; further optimization is possible and would yield a more accurate numerical solution.

`New Folder/master.m` is an orchestrated pipeline that calls all analytical sub-functions in sequence.

## External Dependencies

### RevisitTime-1.0
Coverage revisit-time calculations use the **RevisitTime toolkit** by Crisp & Livadiotti (University of Manchester, 2017–2018), licensed under GPL.

- Source: [University of Manchester RevisitTime](https://github.com/UoM-AAS/RevisitTime) (not included in this repo)
- Functions used: `calcTheta.m`, `listPasses.m`, `ReShape.m`, `RevisitCalc.m`, `RevisitGrid.m`, `RevisitPlots.m`
- Download and place in `MATLAB_Code/RevisitTime-1.0/` before running coverage models.

## Running the Analysis

Open MATLAB, navigate to the desired model directory, and run the `master*.m` script.
Results (PNG heatmaps, CSV tables) are saved to the model's output subdirectories.

---
Date: 11/05/2026
Also a version refactored in Python, improving some functions on MATLAB as well using Claude Code, just for fun
