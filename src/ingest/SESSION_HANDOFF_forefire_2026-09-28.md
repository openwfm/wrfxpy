# ForeFire session handoff — 2026-09-28

The ForeFire pipeline now runs a forecast from a wrfxpy job file with no WRF-SFIRE
anywhere in the loop. Winds and fuel moisture both come from HRRR forecasts.

## 1. TL;DR

- **`src/ingest/forefire_forecast.py` is the missing orchestration** (`089fb97`). Job
  file in, perimeters out. Verified end to end on two fires.
- **All three inputs are HRRR now.** Winds from the forecast cache via WindNinja, fuel
  moisture from FMDA's HRRR region, fuels from LANDFIRE.
- **Moisture varies per step and per fuel category.** On GALENA_RX all nine steps got
  their own fuel table, with the 1-hour class recovering 0.081 -> 0.117 overnight.
- **`etc/forefire.json` moved from `wksp_fmda/CONUS` to `wksp_fmda/CONUS_HRRR`.** CONUS
  is RTMA-driven, and RTMA is an analysis with no hours past the present, so a forward
  forecast silently reverted to the default fuel table partway through its own window.
- **A genuinely forward run works**: ignition at T+30, every wind step from HRRR leads
  valid in the future, every moisture step from FMDA `-fNN` forecast files.
- **Piece 4 (the pSAF re-check) is still open** and is now more urgent, because the
  moisture source has changed under the calibration.

## 2. What was built

`forefire_forecast.py`:

    PYTHONPATH=src:src/ingest python src/ingest/forefire_forecast.py <job.json> \
        [--steps 8] [--issue-utc ...] [--domain-km N] [--overwrite]

It reads the ignition, builds the per-step netcdfs through `forefire_grib`, writes the
`.ff` scripts through `make_script_set` (which applies the per-step moisture), runs
them, and collects the outputs. Every one of those pieces already existed; the only
thing joining them was a hand-written scratch script per fire.

**Job time formats.** Three spellings are in circulation and all three appear in real
files: `2026-09-27_21:08:00` in job files and workspace `input.json`, ISO with a `Z`,
and the relative `T+30` that the templates under `tests/` use. Accepting only the first
is how a template silently becomes a 1970 ignition.

**Domain.** Derived from the job's own grid — `cell_size / subgrid_ratio` is the fire
mesh WRF-SFIRE would have used, `cell_size * domain_size` its extent — so a ForeFire run
and a WRF-SFIRE run of the same job are comparable rather than accidentally different.
When that exceeds WindNinja's measured cell ceiling the **domain** is trimmed, not the
mesh: coarsening the solve degrades exactly the ridge-top winds that drive the fastest
spread, and an 8 h fire does not reach the edge of a 50 km box. A 61 km job at 25 m
(5.95 M cells) trims to 50 km.

## 3. Verified end to end

**GALENA_RX** (Oregon, 2026-09-27 21:08, 31 km at 25 m, 8 steps):

- winds from the 18z cycle leads 3-5, then the 00z cycle leads 0-4 — `find_gribs`
  switching to the fresher cycle as soon as it exists, which is the piece-2 rule working
- nine distinct fuel tables, 1-hour class 0.0807 -> 0.1166 overnight
- **3.1 -> 76.0 ha over 8 hours**, steady growth, plausible for a prescribed burn

**A forward run** from `tests/default_simple_forecast.json`, whose `T+30` resolved to 30
minutes ahead of the clock:

| | |
|---|---|
| wind steps | 06z cycle, leads 6-12, valid 12:00-18:00 — **steps 3-7 ahead of the clock** |
| moisture steps | 10z cycle, leads 2-8 — **every hour from a `-fNN` forecast file** |

That is the whole point of the exercise: nothing in the run waited on a WRF-SFIRE
forecast, and nothing was read from a file that only exists in hindsight.

**The forward fire does not spread, and that is correct.** Fuel 2 sits at Md
0.164-0.171 against its moisture of extinction of 0.150 for the first five steps — too
wet. The sixth step reaches 0.1426, below extinction but only barely, and Rothermel's
rate of spread goes to zero as Md approaches `me`, so no measurable growth in one hour
is the right answer rather than a stalled front. Colorado front range, late September,
1862-2190 m. **I first read this as the front failing to restart after a stall; the
measurements do not support that** — the front is present and valid at every step
through 18:00, and only one step was ever below extinction.

## 4. A correction carried over

Last session's §7 traced a bug where an out-of-range-dry 1-hour moisture class was
replaced by roughly the 100-hour value. The clamp fixed it. Worth keeping straight:
**the fix prevents the stall; it was never shown that a stalled front can restart.**
That question is still open and needs a case where moisture crosses `me` decisively
mid-run, not marginally.

## 5. Open items

1. **Piece 4, the pSAF re-check** — now more urgent. `fmda_dir` has moved to
   `CONUS_HRRR`, so the moisture feeding every fire has changed under a calibration that
   was fitted on RTMA-driven output. Recipe in the 09-25 handoff §8. Baseline to beat:
   geo mean 1.56, median 1.53, 25th percentile 1.11, n=58.
2. **FMDA forecast mode has no schedule.** It is run by hand. For the pipeline to be
   usable on demand the `-fNN` files have to be fresh, which means a cron entry after
   each 00/06/12/18 cycle. **Use a separate flock** — `hrrr_cycler_fmda.lock` is held by
   the hourly analysis run. ~90 s per lead; `forecast_length` is 12.
3. **The last step of a run loses its front** — writes an empty geometry with the
   timestamp reset to the anchor. Seen on Eagle_Springs and again here. Harmless to the
   forecast, but it is the same signature as a genuinely broken chain, so it makes real
   breakage harder to spot.
4. **`ff_growth.track` collapsed seven step geojsons into two rows** on the forward run.
   Did not affect the result — the areas were read from ForeFire's own logs — but it
   would corrupt scoring, so it needs a look before piece 4.
5. Items 3-6 of the 09-24 handoff §9 remain untouched: the fuel 2/9 transposition
   upstream, the false "unknown moisture source" log line, scoring against observations
   rather than only WRF-SFIRE, and the two fires with no overlapping valid time.

## 6. NEXT SESSION

Piece 4. Everything it needs is in place, and the moisture source has changed under it.

Then item 2 above, because an on-demand pipeline that depends on a hand-run FMDA job is
not yet operational.

`etc/forefire.json` is gitignored; `etc/forefire.json.initial` carries the required keys
(`09594ab`) and is the only version-controlled record of what a working install needs.

Commits remain unpushed on `james_ngfs` and `release-fmda-fixes`; JH pushes them
manually (SSH password).
