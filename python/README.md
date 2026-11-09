# Python mission analysis

This later Python extension is **vibecoded**. The [built thesis](../TFG_Miguel_Rosa.pdf)
and its MATLAB analysis remain the primary deliverable.

This later extension can calculate the design from scratch using the same
editable inputs, optical equations, instrument scaling and maintenance budgets
as MATLAB. Orekit supplies propagated ground accesses and station contacts.
The MATLAB pipeline remains the source of the submitted thesis results.

## Independent design search

After installing Orekit, run from the repository root:

```sh
export OREKIT_DATA_PATH=/private/tmp/tfg-orekit/data/orekit-data.zip
/private/tmp/tfg-orekit/env/bin/python python/scripts/run_analysis.py design
```

The `design` command reads the literal inputs directly from
`MATLAB_Code/master.m`; it never reads `selected_design.csv` as a starting
solution. It recalculates all three bands over the 200–1000 km altitude grid,
20–200 mm pupils, 20–500 km swaths and all five constellation configurations.
Detector catalogs, FoV, overlap, thresholds and the eight-year propellant model
come from the same source. Candidate selection retains MATLAB's binding
lambda2 MTF/lambda3 SNR policy and minimum feasible pupil. Ranking includes dry
mass and maintenance propellant; deorbit remains a subsequent budget.

The default `--criterion matlab-mean` retains MATLAB's mean-availability rule:

`effective revisit = worst observed geometric revisit / (1 − cloud fraction)`.

This is a maximum geometric gap with a **mean cloud correction**, not the mean
of all target gaps. The selected cloud fraction is read from `master.m`
(currently `1/6`). Stationary Markov clouds have the same mean clear fraction;
their persistence and extreme-gap distribution are reported as diagnostics.
`--criterion p95` explicitly switches to the stricter empirical 95th percentile
of the worst-target cloudy-access gap. These are different requirements and
must not be compared as if they selected the same quantity.

The orbital/access engine is an explicit model change: Orekit J2/J3 propagated
positions, daylight and timed point accesses on the rectangular study grid
replace RevisitTime's analytical latitude/longitude construction. Its numerical
revisit results need not equal MATLAB, even though the optical/mass equations
and mean cloud adjustment match. Interior geometric gaps are used for selection,
with at least two observations required per target; horizon-censored gaps remain
in the mission diagnostics.

The CLI defaults numerical kernels to one thread to avoid overhead on small
geometry operations; explicit thread environment settings take precedence.

The search visits hardware in increasing maintenance-inclusive mass order.
The widest admissible swath rejects a hardware row when even it fails coverage;
otherwise a monotone search finds the first feasible swath. A failed target
subset is used only as a lower-bound rejection; successful candidates are
checked on the full grid. Propagated orbits are reused between hardware cases.
The winner is verified at a finer ephemeris step, then receives its own
seven-day downlink/memory and eight-year end-of-life budgets.

Outputs are under `python/outputs/design/`: `hardware_candidates.csv`,
`best_design.csv`, `design_summary.json` and the computed mission's
`mission/mission_explorer.html`. A running or failed search is labelled
explicitly. Under the p95 criterion, a regionwide cloudy-spell lower bound may
reject every hardware candidate before orbital evaluations. That is a bound
for the stated Monte Carlo sample, not a calibrated climate impossibility.

## Adopted-design validation case

`config/selected_mission.yaml` is a separate validation experiment. It reads the
adopted thesis design from `Compiled data/Mission/selected_design.csv`; it does
not perform optimization. Use it to compare propagators, sampling steps and
Markov cloud diagnostics at a fixed configuration.

## Run with Orekit

From the repository root, install the JCC wrapper and physical data on macOS
Apple Silicon:

```sh
bash python/scripts/install_orekit.sh
```

The script installs into `/private/tmp/tfg-orekit`; it does not change the system
Python or shell configuration. It needs internet access. On other platforms,
create a conda-forge environment from `python/environment-orekit.yml` and download
the official [Orekit physical data](https://gitlab.orekit.org/orekit/orekit-data).
The conda-forge JCC wrapper supplies `import orekit` and its Java runtime; JPype
is a different binding and is not used here.

```sh
export OREKIT_DATA_PATH=/private/tmp/tfg-orekit/data/orekit-data.zip
/private/tmp/tfg-orekit/env/bin/python python/scripts/run_analysis.py mission
```

This runs an Eckstein–Hechler propagation with explicit J2/J3 zonal coefficients
and mean initial elements, then transforms states into ITRF with IERS 2010
conventions and a WGS84 ellipsoid. Ground-target illumination comes from the
Orekit solar ephemeris. The data path can also be passed with `--orekit-data`.
An unavailable wrapper or data archive causes an error; another propagator is
never substituted automatically.

## Independent numerical reference

A second backend uses SciPy's DOP853 integrator with central, J2 and J3
accelerations. It makes the event-accounting and visualization code usable
without a JVM:

```sh
python python/scripts/run_analysis.py mission --backend scipy-reference
```

The reference uses approximate Greenwich rotation with UTC as UT1, a low-order
solar ephemeris and the thesis mean elements directly as osculating initial
elements. It omits EOP, polar motion, nutation, drag and mean-to-osculating
conversion. Its outputs identify the backend explicitly. They are independent
numerical results, not an executed Orekit run or a replacement for the MATLAB
SGP4 analysis.

Both commands write ignored artifacts into `python/outputs/`:

- `mission_explorer.html`: a standalone interactive globe, target map and
  acquisition/contact timeline. Open it in a browser; no server, CDN or Plotly
  installation is needed.
- `ephemeris.csv`: per-satellite timed Earth-fixed positions and geodetic heights.
- `acquisitions.csv`, `contacts.csv`: interpolated event intervals in seconds
  from the configured UTC epoch.
- `target_metrics.csv`: merged constellation accesses, observed gaps and censored
  horizon-boundary gaps, including targets with no acquisition.
- `mission_summary.json`: configuration, model assumptions, metric definitions,
  input/output provenance and numerical summaries.
- `cloud_summary.json`, `cloud_realizations*.csv`: cloud-adjusted gap statistics
  and the individual Monte Carlo realizations.

Use `--output-dir PATH` to keep multiple experiments separate, or `--no-html` to
export only numerical artifacts.

## What the access metrics mean

The target grid samples the rectangular study-region bounds used in the thesis;
it is not a political CONUS boundary polygon. A nadir-aligned circular sensor
cone approximates the pushbroom footprint. Its effective width is
`240 × (1 − 0.05) = 228 km`, accounting for the adopted overlap. Optical access
requires positive solar elevation at the target. The display's cylindrical
satellite-shadow flag is separate from this illumination condition.

Intervals use cubic interpolation of propagated positions, closest-approach
searches and root refinement to recover grazing passes between samples. The
default ephemeris step is 10 seconds; roots are refined to 1 ms. Compare finer
steps before drawing mission conclusions. The largest sampled no-access gap includes the beginning and end
of the finite horizon. Boundary gaps are censored, so this number neither
certifies a weekly revisit requirement nor establishes a rigorous long-term
bound.

Cloud scenarios apply a seeded two-state Markov chain to the actual acquisition
intervals, splitting intervals at UTC midnight and merging both satellites. The
default fraction is the MATLAB input `1/6`, with persistence 0.5, 200 realizations
and seed 12345. The former Python fraction 0.30 is exported as a separate
sensitivity case. One daily chain is shared across the entire region, assuming
perfect spatial correlation; these are illustrative parameters, not a fitted
weather model. The reported percentile is taken after finding the worst target
in each realization. `cloud_hardware_floor` also evaluates the unavoidable
shared-cloud blackout duration, which can rule out the configured empirical
weekly-gap objective regardless of the hardware.

Contacts use only the thesis NOAA/Fairbanks coordinate (64.84° N, 147.712° W)
above a 10° elevation mask. They do not demonstrate antenna availability,
scheduling, a contracted link rate or an ASF AS1 service.

## Tests and the earlier portfolio model

```sh
PYTHONPATH=python/src /private/tmp/tfg-orekit/env/bin/python -m pytest python/tests
```

Scientific checks cover the adopted design, WGS84 coordinates, the J2/J3
acceleration against an independent potential gradient, retrograde nodal
precession, illumination, event interpolation, constellation interval merging
and censored/no-access targets. Tests that start the real JVM additionally need
the Orekit environment and `OREKIT_DATA_PATH`.

`config/portfolio.yaml` and the `run` subcommand retain an earlier optical and
dry-mass portfolio study. Its coverage uses a latitude-band tiling estimate and
Markov cloud statistics, with a Keplerian Orekit track as a reachability gate.
It does not use the new timed point-access metric, and its grid/catalog differs
from the selected thesis design. Keep its rankings separate from the selected
mission and from the MATLAB optimum.

The offline globe embeds [Natural Earth](https://www.naturalearthdata.com/)
public-domain coastline outlines at 1:110m resolution. These provide geographic
context; they do not change the numerical acquisition region or terrain model.
