# ForeFire session handoff — 2026-10-01

A second forward batch on fresh fires, a measurement bug that had been quietly
subtracting fire fronts since before the pSAF calibration, and an assessment of the
RAP/RRFS grib sources.

**Section numbers continue from the 2026-09-29 handoff on purpose.** §5i-§5k were
written there first and are referenced by number from commit messages and from that
file's open items; renumbering them here would break every one of those references.
Read this file as the 10-01 continuation of `SESSION_HANDOFF_forefire_2026-09-29.md`,
not as a standalone.

## 1. TL;DR

- **13 fresh fires re-run as genuine forecasts** at 50 m, pSAF 0.337, HRRR winds and
  HRRR FMDA, everything cut off at issue time. Out of sample against the 20 of §5d/§5f,
  and the first batch with continuous FMDA forecast coverage. 8 h: geo mean 1.63,
  median 1.87, **36% within 2x** (§5i).
- **§5d's 8-hour prediction fails out of sample** -- it forecast 78% within 2x after
  rescaling to 0.337, against 36% measured -- and its "degraded moisture is suppressing
  these numbers" caveat is refuted with it (§5i).
- **`ff_growth.track()` was subtracting ForeFire's fronts instead of adding them**
  (`0d5ab57`). ForeFire packs separate fronts as rings of one polygon and the GeoJSON
  hole convention removed them. One cron run read 132.4 ha against a true 1621.8. Found
  by JH reading the kml and saying the gap did not look real (§5j).
- **The calibration survives it**: a duplicate ring appears in 1 of the 120 most recent
  cron runs, so pSAF 0.337 does not need refitting (§5j).
- **RAP and RRFS assessed, not ported** -- JH is discussing the choice with others.
  Merging the branch is the wrong tool; the port is six files and two edits, and RAP has
  a cycle/lead defect that would 404 on three cycles in four (§5k).
- **Two new tools**: `ff_fwd_score.py` scores forward runs at both horizons
  (`faef20b`, `3f1a37b`); `ff_moisture_bias.py` carried over from 09-30.
- **Forward kml names are tagged `_fwd`** (`9107510`), so a forward run no longer
  overwrites the cron run of the same fire in the shared `perim_copy_dir`.

## 5i. Second forward batch, 13 fresh fires, pSAF 0.337 (2026-10-01)

Every WRF-SFIRE run that finished after the 09-30 batch was built, re-run as a genuine
forecast: `forefire_forecast.py` from the job file, 50 m with the barrier mask, pSAF
0.337, `--steps 27`, winds from the HRRR forecast cache and moisture from `CONUS_HRRR`,
both cut off at `--issue-utc` = the first whole hour after ignition. No overlap with the
20 fires of §5d/§5f, so this is **out of sample**.

**The FMDA confound of §5d is gone.** The forecast cycler has filled in -- four cycles a
day (01/07/13/19), leads to f32, continuous from 09-28 -- and every fire here ignites on
09-29 or later, so all 13 ran on genuine forecast moisture for their whole window.

| fire | full valid | FF ha | WRF ha | ratio | 8 h | FF ha | WRF ha | ratio |
|---|---|---|---|---|---|---|---|---|
| IA_2 | 10-01 21:00 | 2981.7 | 6570.8 | 0.45 | 10-01 05:00 | 445.3 | 1006.9 | 0.44 |
| LAC-352011 | 10-01 18:00 | 57.3 | 124.1 | 0.46 | 10-01 02:00 | 15.8 | 49.2 | 0.32 |
| Wayne_-_Davis_Road | 10-01 21:00 | 105.2 | 92.6 | 1.14 | 10-01 03:00 | 13.6 | 12.2 | 1.11 |
| STARR_33_RX | 10-01 21:00 | 77.3 | 51.9 | 1.49 | 10-01 05:00 | 12.1 | 17.5 | 0.69 |
| Alcorn-_Blue_Room_Rd | 09-30 09:00 | 337.1 | 203.5 | 1.66 | 09-29 16:00 | 51.2 | 27.4 | 1.87 |
| SOUTH_PARK_RX | 10-02 00:00 | 31.1 | 17.1 | 1.82 | 10-01 05:00 | 6.2 | 3.1 | 2.01 |
| SAN_FRANCISCO_2 | 10-01 09:00 | 1580.5 | 682.9 | 2.31 | 09-30 17:00 | 180.0 | 18.2 | 9.90 |
| IRIS | 10-02 00:00 | 19814.2 | 7684.1 | 2.58 | 10-01 07:00 | 2208.4 | 910.9 | 2.42 |
| Dorran_Road | 10-01 00:00 | 1321.8 | 414.9 | 3.19 | 09-30 07:00 | 416.2 | 222.0 | 1.87 |
| Hatch_Grande | 10-01 00:00 | 1483.8 | 456.8 | 3.25 | 09-30 07:00 | 485.0 | 173.1 | 2.80 |
| GRAPEVINE | 10-02 06:00 | 3292.7 | 686.8 | 4.79 | 10-01 11:00 | 278.2 | 65.4 | 4.25 |

| | n | geo mean | median | below 1x | within 2x |
|---|---|---|---|---|---|
| full window | 11 | 1.69 | 1.82 | 18% | **36%** |
| 8 h | 11 | 1.63 | 1.87 | 27% | **36%** |

### §5d's 8-hour prediction does not hold out of sample

§5d rescaled its pSAF 0.41 batch to 0.337 and predicted geo mean 1.42, median 1.30,
**78% within 2x**, 17% below 1x at the 8 h horizon. Measured on fresh fires: 1.63, 1.87,
**36%**, 27%. The central tendency survives; **the fraction landing within a factor of
two does not** -- it is less than half the predicted value.

**And §5d's hopeful caveat is refuted.** It argued that seven fires with degraded
moisture were dragging the numbers down and "true 8 h skill is probably better than
shown". Coverage is now continuous and skill is **not** better. Degraded FMDA was not
what was holding those numbers down.

What limits the product is the same thing §5h found: **per-fire spread**. The 8 h ratios
run 0.32 to 9.90. No global multiplier reaches that, and this is now the third
independent line -- §5f's sweep, §5h's moisture comparison, and this batch -- pointing at
per-fire treatment rather than a better global constant.

### Two fires extinguished, and FMDA says they should have

`Walker_Hill_Piles_RX` and `Walker_Hill_2_Piles_RX` produced 28 step files of which 26
hold no front at all: ForeFire reached ~0.96 ha and the fire went out. The truncation
guard flagged both, correctly, but the cause is physical rather than a failure. These are
high-elevation Colorado pile burns (terrain 1834-3836 m) and **every fuel sat above its
moisture of extinction** for the whole window:

| fuel | Md used | `me` |
|---|---|---|
| 2 | 0.246-0.270 | 0.15 |
| 5 | 0.212-0.238 | 0.20 |
| 8 | 0.336-0.354 | 0.30 |
| 10 | 0.336-0.354 | 0.25 |

with the 100-h class clamped at the 0.5 ceiling for hours. WRF-SFIRE still burned 2.5 and
11.6 ha. **Both models are arguably right**: a pile burns because it is a pile, not
because the surrounding fuel carries fire, and ForeFire has no pile -- it is being asked
to spread a surface fire through fuel the moisture model says will not carry one. Worth
remembering before treating a prescribed pile burn as a forecast verification case at
all.

### SAN_FRANCISCO_2 is Finley with the sign reversed

9.90 at 8 h against 2.31 at the full window: WRF-SFIRE holds 18.2 ha eight hours in and
reaches 682.9 by the end. §5g's Finley had ForeFire late and WRF early; here WRF is the
slow one. The timing disagreement is **not one-signed**, which matches §5h.

### Two tooling fixes this batch needed

- **`ff_fwd_score.py` (`faef20b`, `3f1a37b`).** Nothing in the repo could score a forward
  run: `ff_batch_ratio.py` reads `<wksp>/forefire`, these land in
  `wksp_forefire/<grid_code>`. It reports both horizons from one run and carries the
  truncation guard. The first version matched timestamps exactly, which is right for a
  cron run but wrong here -- ForeFire stamps a forward perimeter at the second its front
  step landed, so times drift to :01 and those hours were discarded. SOUTH_PARK_RX
  matched 10 hours of 28 and had **no 8 h point at all**. Fixed by snapping to the hour;
  the ratios barely moved (1.68 -> 1.69 full), so the loss was sparse rather than biased.
- **Forward kml names are tagged (`9107510`).** `run_job` now passes `suffix='fwd'`, so a
  forward run writes `<grid_code>_steps_fwd.kml`. Untagged, it wrote the same name as the
  cron run into the shared `perim_copy_dir` and silently replaced it. The 13 runs of this
  batch were retagged after the fact and their untagged duplicates removed, per JH.

## 5j. `track()` was subtracting fronts, and JH caught it by eye (2026-10-01)

JH compared the forward runs against the cron runs of the same fires and read the gap
as "a bit random". Measured, it is not random -- it is tight agreement with outliers:

| | fwd/cron |
|---|---|
| the 8 that agree | geo mean **0.90**, geo sd **1.14**, range 0.68-1.07 |
| SAN_FRANCISCO_2 | 11.94 |
| LAC-352011 | 0.21 |

Adding those two takes the scatter from 1.14 to 2.56. JH then said the SAN_FRANCISCO_2
gap was **not visible in the kml** and was worth re-checking. It was not visible because
it was not there.

### The bug

`ff_growth.track()` measured each perimeter with plain `shape()` + `unary_union`.
ForeFire's `dumpMode=geojson` packs every separate front as a **ring of one polygon**,
and the GeoJSON convention makes ring 0 the shell and the rest holes -- so the second
ring is **subtracted**. `forefire.read_ff_geojson` exists for exactly this (its
docstring records a 106,869 ha step collapsing to 5,857) and `track()` never used it.

On SAN_FRANCISCO_2's cron run a second ring opens at 09-30 17:30:

| step | valid | measured | true |
|---|---|---|---|
| 18 | 09-30 17:00 | 64.2 | 64.2 |
| 19 | 09-30 17:30 | **36.3** | 102.1 |
| 53 | 10-01 09:00 | **132.4** | 1621.8 |

**Under by 12.3x.** A fire cannot shrink, so the track was physically impossible from
step 19 on, and nothing checked that.

### Not every nested ring is an island

The naive repair -- difference the islands out -- is also wrong, and gave 146.9 ha. Over
all 822 step perimeters of this batch the ratio of the largest nested ring to its
container is **sharply bimodal with an empty middle**:

| ratio | steps |
|---|---|
| no nested ring | 720 |
| < 0.2, real islands | 70 |
| 0.2 - 0.8 | **0** |
| > 0.8, duplicates | 32 |

The high group is a front dumped as a thin band: on the last step the two rings are
1621.8 and 1607.4 ha, a mean width of **7.25 m against a 25 m fire mesh**, so it cannot
be unburned ground. `_perimeter()` now differences out only rings below half the burning
area. The empty middle is what makes that threshold evidence rather than tuning.

### What the rescore moves

- **Forward batch: nothing.** All 13 forward runs were unaffected at their scored times,
  so §5i's numbers stand exactly as recorded.
- **Cron batch: one fire in 19.** SAN_FRANCISCO_2 132.4 -> 1621.8 ha, ratio 0.19 -> 2.38.
  Batch geo mean **1.35 -> 1.54**; median 1.81 and within-2x 9 of 19 both unchanged.
- **SAN_FRANCISCO_2 is no longer an outlier.** fwd/cron is **1.10**, in the pack with the
  other eight. The 11.94 was measurement, and the "threshold sensitivity" explanation
  offered for it was wrong.
- **The historical calibration is safe.** A duplicate ring appears in **1 of the 120 most
  recent cron runs**. One fire in 65 moves a geometric mean by a few percent, well inside
  noise, so pSAF 0.337 does not need refitting. (An earlier guess in this session that
  pSAF might be calibrated high followed from assuming the bug was common. It is not.)

### LAC-352011 survives, and suggests a second metric

0.23 after the fix, so that one is real. Per JH it sits in odd terrain near Edwards AFB
with ponds constraining growth in places. The cron run reaches WRF-SFIRE's **final** area
13 h early (step 24) and ends at 268 ha against 124, yet JH reports its **step-33
perimeter is extremely close to WRF-SFIRE's final perimeter**. That is a rate error on a
path both models agree about, and area-at-a-fixed-time cannot distinguish it from a fire
spreading the wrong way. **Time-to-area belongs beside area-at-time** in the scorer.

### The lesson worth keeping

Every area number in this project rests on `track()`, and the defect survived a pSAF
calibration over 65 fires, three handoffs and two scoring tools, because an
under-measured fire looks exactly like a fire that did not spread. It was caught by
someone looking at the map. **A monotonicity check is one line and would have caught it
on the step it happened** -- worth adding to the scorers.

## 5k. RAP and RRFS: what a port would cost, and what is not yet verified (2026-10-01)

**Status: assessed, not done.** JH wants to discuss RAP vs RRFS with others before
committing. Nothing here has been merged or run; this is what that discussion needs.

`RAP.py` and `RRFS.py` live in `/data/jhaley/clean_wrfxpy/wrfxpy/src/ingest` on
**`origin/release`**, not on master, and are absent from `james_ngfs`. They matter
because **NAM is being retired** and ForeFire currently has HRRR as its only validated
wind source.

### Do not merge the branch

`james_ngfs` against `origin/release`: merge base **2023-10-09**, 97 files, +5,079 /
-27,210 lines, ~19 files changed on both sides. That is three years of unrelated
divergence to obtain two modules. **Port the files instead.**

### The port is four copies and two edits

| | |
|---|---|
| copy | `src/ingest/RAP.py`, `src/ingest/RRFS.py` |
| copy | `etc/vtables/{Vtable,METGRID.TBL}.{RAP,RRFS}` -- **absent** from james_ngfs |
| edit | `src/forecast.py` -- dispatch branches |
| edit | `src/ingest/retrieve_gribs.py` -- imports, elif chain, usage string |

Nothing else: james_ngfs has no `available_grib_sources` list and no
`create_job_json.py`.

**The base class already fits.** `grib_forecast.py` in james_ngfs expects exactly the
attributes these set -- `cycle_hours`, `period_hours`, `hours_behind_real_time`,
`grib_forecast_hours_periods`, `remote_url`, `browse_aws` -- the same contract
`HRRR.py` satisfies, and `grib_source.py` already handles `s3://`. Neither module needs
editing to run.

### Verified live against the buckets (2026-10-01)

| path | result |
|---|---|
| `noaa-rap-pds/rap.20260930/rap.t12z.wrfprsf01.grib2` | 200, 226 MB |
| `noaa-rrfs-ops-pds/rrfs.20260930/12/rrfs.t12z.prslev.3km.f001.conus.grib2` | 200, 563 MB |
| `noaa-rrfs-ops-pds/.../rrfs.t12z.2dfld.3km.f001.conus.grib2` | 200, 339 MB |

**Both filename templates are correct.** That is the part most likely to be wrong in an
unused module, and it is not.

### Three defects or unknowns, in order of how early they bite

1. **RAP's cycle/lead pairing is wrong.** `cycle_hours = 1` with
   `grib_forecast_hours_periods = [{'hours': 51}]`, but only the **03/09/15/21Z** cycles
   reach FH51. Measured: `t15z/f51` 200, `t14z/f51` 404, `t14z/f21` 200, `t14z/f22` and
   `f23` 404. So three cycles in four would request files that do not exist. `HRRR.py`
   already solves this exact problem -- `cycle_hours = 6`, "only the 00/06/12/18Z cycles
   reach f48". **Fix to `cycle_hours = 6` for 51 h runs, or `hours: 21` for hourly
   cycles, before any test.**
2. **All four vtables are byte-identical copies of the HRRR ones.** `cmp` confirms
   `Vtable.RAP`, `Vtable.RRFS`, `METGRID.TBL.RAP`, `METGRID.TBL.RRFS` all match
   `Vtable.HRRR` / `METGRID.TBL.HRRR`. `RRFS.py`'s docstring says so and gives the
   `strings` / `ncdump` checks; `RAP.py`'s docstring does not. **Neither is validated.**
3. **The metgrid level counts are guesses.** RAP's docstring says 38 levels while its
   code returns 41; RRFS returns 46 marked "confirm from met_em files after the first
   successful metgrid run". A wrong value fails in **real.exe**, not in ingest, so it
   surfaces late.

### Volume

RAP ~226 MB per forecast hour; RRFS ~902 MB (both products). An 84 h RRFS cycle is
**~75 GB**. `/data` has 61 TB free and `ingest/` already holds 5.3 TB, so a few cycles a
day adds roughly a terabyte a fortnight on top of the HRRR and NAM caches. Affordable,
not free.

### For ForeFire specifically

**RRFS 3 km is the natural NAM replacement** and should feed WindNinja as well as HRRR
does; **RAP at 13 km is coarser than what we already have**, so it is a Canada/backup
option rather than an upgrade. One catch: `forefire_grib.find_gribs` parses
`hrrr.tHHz.wrfprsfNN.grib2` filenames, so ForeFire use needs a parser change beyond the
port. **The WRF-SFIRE path needs nothing extra.**

### Suggested order if it goes ahead

Copy the six files, fix RAP's `cycle_hours`, wire both dispatch points, then pull a
single RAP file end to end through `retrieve_gribs.sh` to prove the plumbing **before
anything touches WPS**. Vtable validation is a separate and slower job needing a real
ungrib/metgrid run.

## 6. What today added to the open list

The full list stays in `SESSION_HANDOFF_forefire_2026-09-29.md` §6 and §7, which this
file does not duplicate. Today added three items to it:

1. **Add a monotonicity check to both scorers.** A fire cannot shrink. One line would
   have caught §5j's defect on the step it happened, instead of after a 65-fire pSAF
   calibration, three handoffs and two scoring tools.
2. **Score time-to-area beside area-at-time** (§5j). LAC-352011's cron run reaches
   WRF-SFIRE's final area 13 h early while tracing what JH reports is nearly its final
   outline. That is a rate error on an agreed path, and the current metric cannot
   distinguish it from a fire spreading the wrong way.
3. **The RAP/RRFS decision** (§5k), pending JH's discussion with others.

Still open from 09-29 and untouched today: scoring the ensemble members rather than
only the chained track, isolating the mesh at 50 m with analysis moisture, the roads
layer for fire breaks, and re-running the moisture-bias test on a larger batch.
