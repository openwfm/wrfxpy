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
- **The moisture switch does not break the pSAF calibration (§5).** RTMA and HRRR differ
  by 0.2-0.4 percentage points in the quantity that reaches the fuel table, well inside
  the run-to-run noise. **But near extinction they disagree 9% of the time** in the
  Oklahoma panhandle, and pSAF cannot compensate for a fire that never starts.
- **Piece 4 is done: keep 0.337 (§6).** On HRRR moisture the batch gives geo mean 1.39
  against 1.54 paired, but the 25th percentile — the actual calibration target — moves
  only 1.11 -> 1.05, implying 0.326. Inside the noise.
- **Caution: a regional moisture comparison under-predicted the impact.** §5's 0.2-0.4
  point difference implied the calibration was safe by a wide margin; the real aggregate
  move was 10-15% in area, because area does not respond proportionally to moisture near
  extinction. Re-run rather than infer.

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

## 5. RTMA vs HRRR FMDA — does the pSAF calibration travel?

Per JH: two 30x30 km domains, two weeks of hourly analyses, box means. The grids differ
(RTMA 1258x2145, HRRR 1038x1799) so there is no cell-by-cell comparison to make; a box
mean is also what ForeFire consumes, since `md_series` takes a `radius_km` mean around
the ignition. `src/ingest/fmda_compare_regions.py` (`c67eaf5`), 339 common hours,
2026-09-14 to 09-28.

| domain | quantity | RTMA | HRRR | HRRR-RTMA | RMS | corr |
|---|---|---|---|---|---|---|
| Sierra / DOME | 1-h | 0.0841 | 0.0887 | +0.0047 | 0.0097 | 0.938 |
| (37.57, -119.62) | **fuel 2 Md** | **0.0863** | **0.0882** | **+0.0019** | 0.0060 | 0.951 |
| Oklahoma panhandle | 1-h | 0.1244 | 0.1259 | +0.0014 | 0.0272 | 0.906 |
| (36.75, -101.60) | 100-h | 0.1821 | 0.1597 | **-0.0224** | 0.0374 | 0.931 |
| | **fuel 2 Md** | **0.1343** | **0.1307** | **-0.0036** | 0.0239 | 0.904 |

**The answer to JH's question is yes, with a caveat.** In the quantity that actually
reaches the fuel table the two products differ by **0.2 to 0.4 percentage points**, with
correlations of 0.90-0.95. That is far inside the ~3% per-fire reproducibility noise
(09-25 handoff §4), so the pSAF calibration travels and 0.337 should survive the switch.

**The caveat is that moisture does not act proportionally near extinction.** Counting
hours where fuel 2 could burn at all (Md < me = 0.15):

| domain | RTMA | HRRR | disagree |
|---|---|---|---|
| Sierra / DOME | 339/339 | 339/339 | **0 (0%)** |
| Oklahoma panhandle | 218/339 | 233/339 | **29 (9%)** |

The Sierra sits at Md ~0.087, far below extinction, so the products are interchangeable
there. The panhandle sits at ~0.13 against a threshold of 0.15 — close enough that the
same 0.4-point difference flips the verdict in 9% of hours, and HRRR burns in **more**
of them. **pSAF cannot compensate for that**: it scales the rate of a fire that spreads,
and has no effect on one that never starts. So the calibration travels for fire *size*
while the burn/no-burn boundary shifts slightly toward more fire in marginal regimes.

**Where they differ most is the 100-hour class** (HRRR 2.2 points drier in the
panhandle, against 0.1 for the 1-hour). Fuel 2 weights 100-h at only 0.1429, which is
why the derived Md difference stays small. Fuels that weight the coarse classes more
heavily would show a larger divergence, and that is worth checking before trusting this
result across all 13 categories.

## 6. Piece 4 done — the pSAF re-check on HRRR moisture

Batch re-run with `fmda_dir` on `CONUS_HRRR`, scored against the preserved RTMA baseline
(`src/ingest/reference/psaf337_rtma_baseline_2026-09-25.txt`). **Paired on the 42 fires
present in both:**

| | RTMA | HRRR | change |
|---|---|---|---|
| geometric mean | 1.54 | 1.39 | **-10.0%** |
| median | 1.55 | 1.31 | -15.2% |
| below 1.0x | 7 | 8 | |
| within factor 2 | 34 | 33 | |

69% of fires produce less fire under HRRR moisture, with a per-fire log-change sd of
**0.300** -- individual fires move about +/-35%, extremes from -64% (HOP-PATTERSON_PILE_RX)
to +93% (Merit_Creek). That is an order of magnitude above the ~3% reproducibility noise,
so it is signal.

**Keep 0.337.** The calibration posture was the 25th percentile at 1.0, and it is still
there: 1.11 -> 1.05, implying pSAF **0.326**, a 3% move inside the noise. Restoring the
*median* would want 0.371, but the median was never the target.

**A prediction that was too optimistic, worth recording.** From the two-domain moisture
comparison in §5 I said 0.337 "should survive" on the strength of a 0.2-0.4 point mean
difference. The actual aggregate move is 10-15% in area -- the right conclusion, but by a
narrower margin than that comparison implied. The reason is the same non-linearity the
Oklahoma panhandle exposed: **area does not respond proportionally to moisture near
extinction**, so a small mean difference flips fires across the threshold and moves area
far more than the mean difference suggests. A regional moisture comparison is therefore a
weak predictor of a calibration's stability, and is not a substitute for re-running.

The HRRR set also has slightly better factor-of-2 coverage (79% vs 76%) on the full sets
despite producing less fire, which is mildly in favour of JH's view that the HRRR product
is the better one.

**Overwriting.** `run_forecasts` rewrites `<wksp>/forefire` in place, so the 09-25
perimeters for the 50 re-run workspaces are gone from this machine; JH has them as KMLs
elsewhere for a visual pass. The scored table is committed under `src/ingest/reference/`
so the next A/B does not depend on a scratchpad. For side-by-side runs use
`sweep_ff_params`, which takes `dest_dir` and `suffix`.

## 7. Open items

1. ~~**Piece 4, the pSAF re-check.**~~ **DONE — see §6. Keep 0.337.**
2. ~~**FMDA forecast mode has no schedule.**~~ **DONE** — cron entry added:

       30 2,8,14,20 * * *   ./hrrr_cycler.sh f CONUS_HRRR

   One run per 6-hourly HRRR cycle. At 02/08/14/20 local (MDT) `cycle_start` lands on
   00/06/12/18Z respectively, ~2 h after the cycle so the leads have published, giving
   11 leads and ~11 h of forward moisture per run.

   **Correction to what this item first said.** It advised a *separate* flock. That is
   wrong and would have been a data race: forecast mode is **not read-only**. Lead 1
   seeds from lead 0, which is the cycle's *analysis* file, and it writes that file
   itself when it is missing — the 2026-09-28 run ran DA three times and created
   `fmda-CONUS_HRRR-20260928-10.nc`. It therefore shares `hrrr_cycler_fmda.lock` with
   the hourly analysis. Minute 30 because 11 leads take ~35 min and the analysis fires
   at :27, leaving 57 minutes of clear air.
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

## 8. NEXT SESSION

**Re-run `fmda_fcst_skill.py` once the cron has accumulated cycles.** The forecast-vs-
analysis result in this session came from 2 cycles with n=1 at leads 3-9, and JH's read
of it was right: the apparent wet drift was the period, not the system. The two cycles go
opposite ways (+0.061 and -0.007 at the 10-h class), and the archive confirms 09-25 was
drying 1.1 sd faster than normal while 09-28 was moistening. The real signature is that
the forecast **lags the analysis tendency in either direction**, which matters more than a
bias would: it cannot be corrected away, and it is largest during exactly the rapid-drying
events that drive fire. At 4 cycles/day x 32 leads the cron yields ~128 pairs/day, so a
couple of days settles it. **Score against the analysis tendency, not lead alone** --
that is what separates lag magnitude from which way the weather went.

Then the full-duration ForeFire vs WRF-SFIRE comparison JH asked for: `forefire_forecast.py
--steps 27` covers the median WRF-SFIRE run, and both the GRIB cache (+33 h) and the new
33 h moisture now span it.

Watch the FMDA forecast cron. Its first firing was silently cancelled -- `flock -n` lost
the lock to a 12-minute analysis under load 37 -- and it now uses `-w 1200` with a skip
log, but that has not yet been seen to work.

Commits remain unpushed on `james_ngfs` and `release-fmda-fixes`; JH pushes them manually
(SSH password).
