# wrfxpy history review

Review date: 2026-11-20.

Repository reviewed: `/data002/projects/dshield/wrfxpy-v3.2`.
Branch: `dshield-v3.2`.
Reviewed HEAD: `1c9f95f839ff8bfbaffb85a11b6abc8574e8909d`.
Jonathon's v3.2 baseline: `ac6fea8ac570da4f8ba98ef6166a9ca5ee816392`.

This is a read-only source-history review, not a numerical comparison of
historical runs. Commit authors identify recorded authorship, not necessarily
who proposed or approved the scientific or architectural decisions.

## Fuel moisture: what changed after Martin Vejmelka

The first subsequent commit touching the FMDA code and Cycler is `fed2edb`,
Jan Mandel, 2016-09-18, "adding prints". The first subsequent change to
`src/fmda/fuel_moisture_model.py` itself is `6248e1a`, 2016-09-20.

Martin's time-lag dynamics and operational Kalman update remain intact.
Comparing `fed2edb^` with the reviewed HEAD, `advance_model` and
`kalman_update_single2` have identical Python syntax trees after removing
logging calls. Later substantive changes concern forcing, spatial estimation,
parameter settings, and additional RNN and blended products.

| Date / commits | Author | Change | Classification |
| --- | --- | --- | --- |
| 2016: `fed2edb`, `6248e1a` | Jan | Logging and diagnostics. | Instrumentation; no method or parameter change. |
| 2018: `bccd9aa`, `c0ce40f` | Jan | Empty observation-service response handling; meteorological fields stored with model state. | Operational fixes and additional output. |
| 2018: `4565083`, `77c9405` | Jan | Rain above `1e10` replaced with zero; rain filtering moved before model advancement. Existing small-rain cutoff was `0.01`. | Input treatment affecting results, not just refactoring. |
| 2019: `de9c822`, `a6bb3cf`, `533116a` | Jan | Geogrid export, longitude/latitude tracer planes, central geographic anchor. | New coupling/output capability and geometry correction. |
| 2020-2021: `341af78`, `ae5abe3`, `0c38fce`, `a296b69`, `2f67754` | Angel | Python modernization, module organization, export generalization, WPS intermediate output. | Refactoring and interfaces; WPS export adds a capability. |
| 2021-2022: including `ce3db6f`, `150b658` | Angel | Observation database integration and missing-data handling. Masked rainfall also became zero (`c60f6cc`, 2021). | Data handling and operations; inputs can change without changing equations. |
| 2022: `003d265` | Angel | Added `fit_tsm_lstsq`; restart reconstruction generalized to the stored number of fuel classes; observation station IDs added. | Numerical method addition, restart generalization, and metadata. |
| 2024-10-15: `903d401`, `21a483e` | Angel | HRRR cycling, four classes, more spatial covariates, different process noise, least-squares solver enabled. RTMA time-constant units corrected. | Method, forcing, and parameter changes. |
| 2024-10-16 to 30: `dc0c5db`, `5b71e9b`, `88d10a4` | Angel | Forecasting without future observation assimilation; observation cutoff and retrieval-window corrections. | Workflow and observation-selection changes. |
| 2025-2026: `b36c851`, `59d6680`, `49e3e17` | Angel | Historical interval runs, NumPy scalar/array compatibility, additional covariate output. | Workflow, compatibility, and output. |
| 2026 June-August: `6a34f91`, `55fce86`, `309e9c0`, `06a29e3` | Jonathon | RNN predictions after generating a meteorological sequence; time-warped predictions for other classes. | New prediction method. |
| 2026 September-October: `73c7a84`, `a588fd3`, `b911b67` | Jonathon | RNN machinery moved into `ml_fmda`; wrfxpy wrapper, recurrent-state persistence, hourly integration. | Package extraction is refactoring; state continuation changes execution/restart behavior. |
| 2026-10-14: `b911b67` | Jonathon | Arithmetic mean of FMDA and RNN moisture predictions. | New combination method with fixed 50/50 weights. |
| 2026-11-04: `146d6bd` | Jonathon | RNN interval operation; WPS intermediate writes disabled in the Cycler. | Workflow/output changes. |
| `1c9f95f` | Geometry correction cherry-picked onto v3.2 | HRRR spacing adjusted for the geogrid reference sphere. | Spatial-coordinate correction, not moisture-model tuning. |

### Parameters and units

| Quantity | Earlier implementation | Current HRRR implementation |
| --- | --- | --- |
| Fuel classes | RTMA: 1, 10, 100 hours | 1, 10, 100, 1000 hours |
| Process-noise diagonal, `Q` | RTMA: `[1e-4, 5e-5, 1e-5, 1e-6, 1e-6]` | `[1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 1e-5]` |
| Initial covariance, `P0` | Three moisture variances of `.01`, two correction variances of `.001` | Four moisture variances of `.01`, same correction variances |
| HRRR accepted observations | Initially `0 < FMC < 50` | Corrected to `0 < FMC < 0.5` in `5b71e9b`; FMC is fractional |
| HRRR minimum count | Introduced as `count > 10` | Still strictly more than 10 |
| Blend weights | No blend | FMDA `.5`, RNN `.5` |

The state order for `Q` is moisture classes followed by `dE` and `dS`.
For shared classes and correction states, HRRR uses ten times the RTMA
process-noise variance. That changes covariance evolution and assimilation
strength. These settings were already in Angel's initial HRRR Cycler.

Four classes already existed in Martin's direct WRF assimilation path.
The HRRR Cycler added the fourth class relative to RTMA, not to the repository
for the first time. The direct WRF path retains its own process-noise settings.

The RTMA initializer originally passed `[1, 10, 100]` to a model expecting
seconds. `903d401` added `* 3600`. NetCDF restart loading already reconstructed
time constants in seconds, so the defect particularly affected continuation
of a freshly initialized model in memory, such as recursive spin-up. It did
not necessarily affect every restarted hourly run.

Unchanged physical constants are the rain threshold `.05 mm/h`, rain-intensity
scale `8 mm/h`, rain-wetting time constant `14 hours`, and saturation moisture
`2.5`. Equilibrium formulas, the observation-variance formula, and trend-surface
clipping to `[0, 2.5]` also remain inherited from Martin.

### Spatial estimation and observation selection

HRRR adds latitude, longitude, temperature, humidity, pressure, soil variables,
snow, wind, and radiation to the spatial covariates. This changes the statistical
model. The model's FM10 field remains a covariate; the number of covariates used
is limited by the number of observations.

`fit_tsm_lstsq` was added in `003d265` in 2022 and wired into assimilation in
`903d401` in 2024. It replaces selected direct solves with least-squares solves,
changing behavior for singular or poorly conditioned systems. Its
equal-observation-variance branch still uses `np.linalg.solve`, so it is not
uniformly protected against singular matrices. The original `fit_tsm` remains.

`5b71e9b` temporarily changed HRRR observation retrieval from a centered
plus/minus 30-minute window to the preceding hour. `88d10a4` restored the
centered window. The actual assimilation routine still selects observations
strictly within 30 minutes of model time. The Cycler's preliminary count/filter
uses 60 minutes, so these are distinct checks.

### Jonathon's RNN and blend

Jonathon added a learned model and operational integration while retaining
Martin's FMDA equations. `src/fmda/fuel_moisture_rnn.py` wraps external
`ml_fmda` code with trained weights, a scaler, and time-warp parameters.
Their numerical evolution cannot be established from this wrfxpy branch alone.

The blend averages four moisture fields after FMDA assimilation, retains
FMDA's `dE` and `dS`, and exports no blended covariance. A code TODO explicitly
questions the credibility of retaining those correction terms with the blended
moisture state. The blend is an output product; the Cycler continues the
separate FMDA and RNN states.

`1c9f95f` changes geogrid spacing from `3000` to approximately `2999.4213 m`
using `3000 * 6370000 / 6371229`. This accounts for different reference spheres
and corrects where moisture is mapped into WPS, without retuning its dynamics.

## Python 3 transition

Angel Farguell Caus introduced Python 2/3 compatibility in `341af78`,
2020-01-07, "transition to compatibility with python 2 and 3 using
python-modernize routine". This was compatibility with both versions, not an
immediate Python-3-only conversion. Later cleanup included `8fc7556`,
2023-07-17, removing Python 2 compatibility code from `forecast.py`.

## JSON workflow-state history

The persisted job description and explicit workflow state appeared in stages
in Jan Mandel's commits to `src/forecast.py`:

| Date | Commit | Change |
| --- | --- | --- |
| 2016-10-30 | `16186e7` | "separating process_output": wrote `job.json` with `job_id`, `postproc`, and `grid_code`; separate postprocessing loaded it. |
| 2016-10-30 | `e1e5bfd` | Used the name `job_process_output.json` while adding the separate processing script. |
| 2016-11-11 | `fa079b4` | Stored scheduler job number and Python process ID. |
| 2016-11-11 | `24c014c` | Renamed `job_process_output.json` back to `job.json`. |
| 2016-11-13 | `919c019` | Explicit `Running`, `Processing`, and `Completed` states, updated during execution. |
| 2016-11-18 | `aef5951` | Process creation time, actual process/scheduler checks, and states such as `Stopped` and `WRF stopped`. |
| 2016-11-25 | `7912d1a` | Extracted initialization into `make_job_file`, with initial state `Preparing`. |
| 2016-11-26 | `4e588f6` | Stored original forecast arguments separately as `input.json`. |

Thus `16186e7` introduced the persisted job description; `919c019` introduced
explicit workflow state in that JSON. Git records Jan as author but does not
establish who originally proposed the design.

The basic design has remained a mutable `job.json` storing coarse status,
process IDs, scheduler job number, and postprocessing information. Later
extensions include satellite-only processing, `Postprocessing failed` status,
restart information, domain count, and station-output settings.

There was no major redesign into a general workflow engine. The file does not
track every preprocessing step or provide automatic resumption from an
arbitrary completed step. Current `make_job_file` remains recognizably the
November 2016 implementation with additional fields.

## Source locations

All paths below are relative to the reviewed wrfxpy repository, not this
notes repository:

- `src/fmda/fuel_moisture_model.py`: time-lag model, Kalman updates, exports, blend.
- `src/fmda/trend_surface_model.py`: original and least-squares spatial solvers.
- `src/fmda/fuel_moisture_da.py`: observation ingestion and assimilation.
- `src/fmda/fm10_observation.py`: observation values and variance model.
- `src/fmda/fuel_moisture_rnn.py`: RNN wrapper and state/output persistence.
- `src/rtma_cycler.py`, `src/hrrr_cycler.py`: forcing, parameters, cycling policy.
- `src/forecast.py`, `src/cleanup.py`: job JSON creation, updates, and status checks.
