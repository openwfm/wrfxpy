# ForeFire session handoff — 2026-09-25

Started as the scoring pass the 09-24 handoff's §10 asked for. It closed the pSAF
question (§2), then turned up a moisture bug that had been silently extinguishing fires
(§7), and ended with three of the four pieces of HRRR-forecast-driven forecasting built
and running (§8). Five bugs fixed in total.

## 1. TL;DR

- **FMDA moisture was silently extinguishing fires (§7).** An out-of-range-*dry* 1-hour
  class was replaced by roughly the 100-hour value, so the driest fuels came out the
  wettest. Eagle_Springs held 0.946 ha for 24 h where WRF-SFIRE burned 1297. Fixed
  (`3610523`); it now scores 0.84-0.89. **This is the most important thing in this handoff.**
- **pSAF 0.337 is verified and survives the moisture fix.** 65 contaminated fires gave
  geo mean 1.55 / median 1.53; the re-run on clean output (58 fires) gives **1.56 / 1.53**.
  The tail moved where it should — the 25th percentile went from 1.00 to 1.11, implying
  pSAF 0.315 against the 0.337 in use, a 7% difference that sits inside the noise.
  **Keep 0.337.** The bug was severe per fire and immaterial to the global constant.
- **The `pSAF^1.54` exponent holds over a 1.9x move**, as it did over the 0.6 -> 0.175
  move. It can be trusted for future adjustments.
- **The templates now explain the number.** Both `.ff` templates carry a comment block
  recording the calibration, so 0.337 does not read as a magic constant.
- **Correction: the size dependence has not vanished.** It is `area^-0.198` over fires of
  10 ha or more — half the old `-0.40`, still real. An all-fires fit says `-0.050`, and
  that is an artifact of two tiny fires.
- **`ff_batch_ratio.py` could silently score the whole archive** when given a bare number
  of hours that matched a filename. Fixed (§6). It did not affect any number here.
- **Reruns are not reproducible to better than ~3%**, so differences below about 5%
  between ForeFire configurations are noise (§4).
- **Real forecasting is now unblocked (§8).** HRRR and HRRRA are the same product, so
  forecast cycles need no WindNinja change. `retrieve_gribs.sh HRRR` was dead on an
  `if`/`elif` slip and is fixed; `find_gribs` now picks a cycle deterministically and
  honours an `issue_utc` so hindcast scores cannot use cycles issued after the fire
  started. **FMDA forecast mode now runs** -- 9 leads, 25 min, no errors -- and ForeFire
  can read its `-fNN` output. Pieces 1, 2 and 3 done; only the pSAF re-check is left.
- **The 0.02 moisture floor is hit routinely, not exceptionally (§8).** 2.18% of CONUS
  by late afternoon, 14.75% of the Great Basin -- and part of that field is *negative*,
  so the floor is rejecting filter undershoot rather than implausible dryness. JH
  confirms the regions match the web server's moisture maps.
- **Two `retrieve_gribs.py` bugs, one per repo.** The `if`/`elif` slip above, and in the
  release repo `== 'NAM' or 'NAM218'` -- a truthy string that silently routed NAM227,
  CFSR, NARR, GFSA, GFSF, RAP and every RRFS_* to NAM218 (`e1444cc`). RRFS is a candidate
  NAM replacement, so that one would have bitten exactly when it mattered.

## 2. The result

Scored with `src/ingest/ff_batch_ratio.py`, ForeFire against WRF-SFIRE burned area at the
latest valid time both models reach:

    cd /data/jhaley/wrfxpy && PYTHONPATH=src:src/ingest \
      python src/ingest/ff_batch_ratio.py 18

The argument is hours. `ff_batch_ratio.py` normally takes the batch script and uses its
mtime to bound the window, but `/tmp/psaf_work/psaf337_batch.py` had been cleaned up. The
templates' own mtime is the same boundary — the moment pSAF became 0.337 — so hours from
that works and is not a guess. **If a batch script still exists, pass it**; the reason
that option exists is that a guessed window let two stale pSAF-0.6 fires into the 0.175
table and moved its geometric mean from 0.531 to 0.577.

| | predicted (09-24 §2) | measured |
|---|---|---|
| geometric mean | 1.46 | **1.55** |
| median | 1.54 | **1.53** |
| below 1.0x | 26% | **17%** |
| within a factor of 2 | — | **75%** |
| n | 57 | **65** |

The median is essentially exact. The batch runs slightly more conservative than aimed
for, which is the harmless direction given JH's preference for overprediction.

### Why stop at 0.337

Rescaling the scored ratios by the exponent shows it sitting at the peak of the
within-2x curve:

| pSAF | geo mean | median | below 1.0x | within 2x |
|------|----------|--------|-----------|-----------|
| 0.300 | 1.30 | 1.29 | 26% | 72% |
| **0.337** | **1.55** | **1.53** | **15%** | **75%** |
| 0.380 | 1.85 | 1.83 | 11% | 60% |
| 0.420 | 2.15 | 2.13 |  9% | 43% |

Moving up to 0.38 buys 6 points fewer underpredictions and costs 15 points of accuracy.
There is no argument for further global tuning.

## 3. Correction — the size dependence is halved, not gone

I reported mid-session that it had collapsed (`ratio ~ area^0.02`). **That was wrong.**
The all-fires fit is flat because two sub-10-ha fires (0.7 and 1.3 ha) sit at the far
left of a log axis and act as leverage points. Over the 63 fires of 10 ha or more:

    ratio ~ area^-0.198     log-log r = -0.478

against `area^-0.40` at pSAF 0.6. Real, and about half as steep.

**The residual is not a clean power law**, which matters before anyone fits one:

| WRF area | n | geo mean ratio |
|---|---|---|
| <200 ha | 20 | 1.71 |
| 200-600 ha | 25 | 1.77 |
| >600 ha | 20 | 1.18 |

Non-monotonic. Small and mid-size fires behave alike; it is the **large** fires that are
already close to right. A single exponent describes that shape badly.

**Method note worth keeping.** Both errors this session — this one and the 0.175/0.337
mix-up in §5 — came from trusting a summary statistic without checking what was in the
sample. A log-log slope over a range that spans 0.7 ha to 6000 ha is dominated by its
endpoints; bin the data before believing the fit.

## 4. Template comment

`/home/jhaley/forefire/tests/ignition_template.ff` and `restart_template.ff` — **outside
the git repo**, referenced from `etc/forefire.json` under `templates.ignition` and
`templates.restart`. This handoff is therefore the only version-controlled record of what
they say. Both now carry, directly above the parameter:

    # Calibrated against WRF-SFIRE burned area over a 65-fire CONUS batch, 2026-09-24/25.
    # 0.6 (ForeFire default) overpredicted 3.85x; 0.175 overcorrected to 0.58x.  0.337 puts
    # the 25th percentile at 1.0x: geo mean 1.55, median 1.53, 17% of fires below observed,
    # 75% within a factor of two.  Area responds as pSAF^1.5, so changes here are steep.
    # Deliberately biased high -- per JH, overprediction is the safer error for a forecast.
    # Caveat: this is model-vs-model.  WRF-SFIRE is not truth; on Dome it ran 0.68x an IR
    # perimeter.  Re-derive against observations before treating 0.337 as physical.

Inserted with a `sed` address matched on the parameter *name*, not a line number, so it
stays put if the file is reordered. The two templates still differ only in
`minSpeed=0.001` and the ignition-vs-restart block, which is how they should differ.

Verified by re-running `RX_Stray_Creek__LCRD` (27 ha) under the cron `flock`: 51 steps,
all "Success: Process completed", no "wrong input file, check your settings" — the
generated `.ff` scripts parse and the step chain completes.

**The rerun did not reproduce exactly, and that is worth knowing.** Same pSAF, same
inputs, same valid time (09-24 21:00): 270.8 ha against the 261.8 ha scored a few hours
earlier, a **3.4%** difference, ratio 9.96 against 9.63. Not diagnosed. The likely
candidate is the hourly FMDA files — `md_series` reads them per step and the FMDA cycler
rewrites them, so a rerun can pick up revised moisture for hours already simulated;
front-tracking nondeterminism is the other possibility.

Either way it sets a noise floor: **differences below about 5% between two ForeFire
configurations are not evidence of anything.** The mesh-insensitivity result from 09-22
(<2% across 25/50/100 m) sits below that floor and should be re-read as "no detectable
effect", not "a 2% effect".

## 5. Correction — which value was actually running

Worth recording because it cost time. JH's instruction this morning was "leave it set at
0.175", but the templates had been at **0.337** since 11:04 on 09-24 — 0.175 was scored,
found to underpredict (geo mean 0.577, 90% of fires below 1.0x) and raised the same day.
The cron forecasts JH looked at and called "a big improvement" were the 0.337 ones, so
the instruction and the observation agreed; only the label was off.

I also mis-scored the batch twice before getting it right, both times by reading the
wrong files:

- the per-step geojsons in the **workspace root** are bare GeoJSON Features, not
  FeatureCollections, so a reader expecting `features` finds nothing;
- the combined `<name>.geojson` in the root is **WRF-SFIRE's own** perimeter export, not
  ForeFire's. Scoring against it compares WRF-SFIRE with itself and returns a suspiciously
  tidy 1.01-1.19 — the ~5% gap between TIGN_G and the polygonised perimeter.

ForeFire's perimeters are the geojsons in `<wksp>/forefire/`. `ff_batch_ratio.py` already
reads the right ones; the lesson is to use it rather than writing another ad-hoc scorer.

## 6. Bug fixed — `ff_batch_ratio.py` could silently score the whole archive

Its argument is "a batch script, or a number of hours", and it decided which by calling
`osp.exists()` **first**. A stray executable named `1` has been sitting in the repo root
since 2025-12-04, so `ff_batch_ratio.py 1` resolved to that file, took its mtime as the
window start, and scored **125 fires spanning three different pSAF values** — reported,
with a straight face, as the last hour. Geometric mean 1.52, median 1.26: plausible
enough to be believed.

Fixed by trying `float()` first and only falling back to a path, with an explicit error
if the argument is neither. A bare number can no longer be captured by a file that
happens to share its name.

The `18.01`-hour run that produced every number in §2 was never affected — nothing is
named `18.01` — and re-running the 1 h check after the fix returns exactly the one fire
that was re-run. **Worth deleting `/data/jhaley/wrfxpy/1`**; it is an old scratch script,
but it is not this tool's job to be robust against it and now it is.

## 7. The big one — FMDA moisture was silently extinguishing fires

JH flagged Eagle_Springs as a run whose ForeFire steps "didn't complete properly". They
completed: 54 steps, no crash, no parse error. **The fire was extinguished by its own
moisture input** and sat at 0.946 ha for 24 hours. WRF-SFIRE burned 1297 ha.

### The chain

FMDA reported the 1-hour dead fuel moisture at the fire as very dry — 0.0113, 0.0076,
0.0051, 0.0042 over 21:00-00:00 UTC. All are under the `valid_range` floor of 0.02, so
`md_series` set that class to None. `write_fuel_table` then filled the gap with **the
unweighted mean of the classes that survived**:

    (0.0798 + 0.3463) / 2 = 0.2130      <- matches the written table filename exactly

Fuel 2 — 86% of the domain and the fuel under the fire — carries cawfe weights
(0.5714, 0.2857, 0.1429), so:

    0.5714*0.2130 + 0.2857*0.0798 + 0.1429*0.3463 = 0.1940   against me = 0.15

Above extinction. No spread, anywhere, for ten steps. By the time FMDA rose back over the
floor at 02:00 the front had been frozen for five hours and never recovered.

| treatment | fuel 2 Md | outcome |
|---|---|---|
| all three classes, as FMDA gave them | 0.0787 | burns |
| **1-h dropped, filled from the mean** | **0.1940** | **dead** — what ran |
| 1-h dropped, SAV weights renormalised | 0.1686 | still dead |
| **1-h clamped to the 0.02 floor** | **0.0837** | burns — the fix |

Wind and fuel were both exonerated first: fuel under the fire is category 2, and the winds
(0.5-1.2 m/s) were *higher* than RX_Stray_Creek's, which burned 270 ha. The controlled
test was re-running the identical fire with moisture off — **1100 ha**, against
WRF-SFIRE's 1297. Nothing was wrong except the moisture.

### Why this is the dangerous direction

It fires precisely when the fine fuels are driest. **Per JH the pattern here — 1-h and
10-h dry, 100-h wet — is the signature of a sharp temperature rise and humidity drop
shortly before ignition.** The 100-hour class lags by design and is still carrying
pre-drying moisture, so borrowing from it during a rapid-drying event imports stale wet
data into the class that matters most for spread. WRF-SFIRE's own `FMC_GC` agrees: its
1-hour class sits at 0.103 and never exceeds 0.116.

This also sharpens JH's standing rule that moisture in the 5-6% range predicts big
growth: the bug converted exactly that regime into no fire at all.

### The fix

`md_series` now **clamps** to `valid_range` and logs every clamp. `write_fuel_table` still
has to fill a genuinely absent class, but takes the **nearest** class rather than the mean
of all of them — neighbouring size classes track each other, the 1-h and 100-h do not.
Its comment claimed "the SAV mix" while the code took a plain unweighted mean; both now
agree. Commit `3610523`.

Eagle_Springs re-run with the fix: **1151.9 ha against WRF-SFIRE's 1297.2 — ratio 0.89**,
where it previously could not be scored at all. It scored 0.84 again in the full batch
re-run below; the 0.84/0.89 gap between two runs of the same fire is the reproducibility
noise from §4, not a difference in treatment.

### It contaminated the calibration

Only 9 batch fires still have their `fueltables` directory, and **4 of those 9 are hit** —
detectable because the substituted value equals the mean of the other two:

| fire | tables affected | ratio in the 0.337 table |
|---|---|---|
| Sheep_Station_Rx | 4 of 11 | **0.48** — lowest in the batch |
| Ranger_Academy_Burn_1_RX | 59 of 1532 | **0.59** — fourth lowest |
| GRADE | 4 of 11 | 0.94 |
| Eagle_Springs | 5 of 9 | unscorable |

**pSAF 0.337 was chosen by putting the 25th percentile at 1.0 — keyed on exactly the tail
this bug creates**, so the batch was re-run with the fix and re-scored.

### Result of the re-run: the calibration holds

Seven fires clamped across the 62 workspaces — Eagle_Springs, GRADE, Sheep_Station_Rx,
Anderson_Butte_Rd_MM4, HOP-PATTERSON_PILE_RX, Horsefly_605_RX and STEEL_PASS — three more
than the surviving fuel tables had revealed.

| | contaminated (n=65) | **clean (n=58)** |
|---|---|---|
| geometric mean | 1.55 | **1.56** |
| median | 1.53 | **1.53** |
| below 1.0x | 17% | **21%** |
| within a factor of 2 | 75% | **76%** |
| 25th percentile | 1.00 | **1.11** |
| size dependence (>=10 ha) | area^-0.198 | **area^-0.217** |

Per fire the fix is dramatic; in aggregate it is not, because only 7 of 62 fires were
touched:

| fire | before | after |
|---|---|---|
| Eagle_Springs | 0.946 ha, unscorable | **0.84** |
| Ranger_Academy_Burn_1_RX | 0.59 | **0.82** |
| Sheep_Station_Rx | 0.48 | **0.58** |
| GRADE | 0.94 | **0.98** |

The 25th percentile moved from 1.00 to 1.11, which implies **pSAF 0.315** where 0.337 is
in use — a 7% change, inside both the ~3% per-fire reproducibility noise (§4) and the
sampling noise on a 58-fire percentile. **0.337 stands**, and it errs on the
overprediction side of 0.315, which is the side JH asked for.

The honest reading: this bug badly damages individual forecasts and barely moves a
constant fitted across dozens of them. Both facts matter — the second is why the
calibration did not have to be redone, the first is why the bug was worth finding.

### Two smaller things in the same run

- **Per-step moisture stops partway.** From step 27 of 54 the fuel table reverted to the
  base `ff_fuels_behave13.csv` at Md=0.1, because FMDA is an *analysis* and has no files
  past the present while the forecast runs to 09-26. Inherent, not a bug, but it means a
  forward forecast loses its moisture signal partway — worth thinking about for the
  8-hour low-latency target, where most of the window is in the future.
- **The last step writes an empty front** with the timestamp reset to the anchor.
  Cosmetic here, but it is the same signature as a genuinely broken chain, so it makes
  real breakage harder to spot.

## 8. Driving ForeFire from HRRR *forecasts* — pieces 1, 2 and 3 done

Per JH: `hrrr_cycler.py` can drive FMDA from HRRR forecasts, and the same cycles should
drive the WindNinja interpolation. NAM218 is deprecated in a few weeks and caching moves
to HRRR forecast cycles at 00/06/12/18. That is what turns this into real forecasting
rather than hindcasting.

**The wind side was almost already built.** `HRRR` and `HRRRA` retrieve the *same*
product — `hrrr.tHHz.wrfprsfNN.grib2`. HRRRA is just "always lead 3, from the cycle
three hours back". Same grid, same variables, so WindNinja needs nothing. And
`_grib_valid_time` already computed **cycle + lead**, `grib_dir` was already a parameter,
the cache layout `hrrr.YYYYMMDD/conus/` already matches the glob, and `hrrr_cycler.py`
never deletes GRIBs.

### Two bugs found on the way in (`16161c7`)

`retrieve_gribs.py` had `if grib_src_name == 'HRRR'` followed by a second **`if`** rather
than `elif`, so the HRRR branch set the source and then fell through the HRRR_AK/NAM
chain, matched nothing and hit its `else`. **`./retrieve_gribs.sh HRRR ...` always died
with "Invalid GRIB source HRRR"** despite HRRR being fully supported. Since
`cache_grib_files.py` drives caching through that script, switching its source list to
HRRR would have failed on the first call.

Its success report also iterated the returned *manifest dict*, printing keys like
`ingest/grib_files` and `ingest/colmet_prefix` — which read exactly like paths and made
an empty download look like a successful one.

Verified after the fix: it retrieved `hrrr.t12z.wrfprsf09/f10`, valid 21:00 and 22:00 UTC
with the clock at 15:40 UTC — real future data — and `find_gribs` mapped them correctly
with no change.

### Piece 2: deterministic cycle choice and an issue time (`42fa193`)

`find_gribs` kept whichever file `glob` returned first, commented "a later lead wins
nothing, so keep the first". True of the analysis cache, where each valid time has
exactly one file. **False as soon as forecast cycles are cached**: 12z+f09 and 18z+f03
are both valid at 21:00 and differ by six hours of forecast error.

For a fixed valid time, newest cycle and shortest lead are the same ordering (valid =
cycle + lead), so one rule covers both: **smallest lead wins**, ties broken on the path
so repeated calls agree.

**`issue_utc` is the part that matters for verification.** Without it, "forecasting" a
fire from last week would quietly use cycles issued *after* the fire started, and every
hindcast score would be optimistic — which matters because the 8-hour product will be
judged on exactly those scores. With it, only cycles published by that moment are
eligible: `cycle + issue_lag_h <= issue_utc`, the lag defaulting to 1 h to match
`HRRR.hours_behind_real_time`. A valid time whose every candidate post-dates the issue
time is **dropped, not back-filled** — it genuinely could not have been forecast.

Default `issue_utc=None` keeps the old behaviour, which is what a retrospective analysis
run wants. Exposed as `--issue-utc` on `forefire_grib.py` and as a `build_step_ncs`
keyword.

Tested against a synthetic three-cycle cache: the publication boundary is exact (18:59 ->
12z f09, 19:00 -> 18z f03), `issue_lag_h` works as a knob, selection is stable across
calls, the real HRRRA cache selects byte-identically to before, and a full
`build_step_ncs` run through WindNinja is unchanged.

### Piece 1: HRRR cycles cached, and the cache finally has retention (`d979f2e`, `207eb8c`)

**The pinning trap was real.** `cache_grib_files.py` computes 00/06/12/18 and hands
`retrieve_gribs` the *range*; NAM218 lands on those cycles by itself because its
`cycle_hours` is 6. **HRRR's is 1**, so the same call takes whatever hourly cycle is
newest -- 14z at 15:40Z -- and the whole 00/06/12/18 scheme would have been lost with
nothing failing. `retrieve_gribs` has always accepted `cycle_start` (hrrr_cycler passes
it) but the CLI did not expose it, and `retrieve_gribs.sh` forwarded only `$1..$4`, so a
fifth argument vanished silently. Both fixed; `PINNED_CYCLE = {'HRRR'}` pins HRRR and
leaves NAM's behaviour untouched.

Verified against the live archive: pinned to 06z for a window valid 18:00-19:00Z it
fetches `t06z.wrfprsf12/f13`, where unpinned it took ~14z.

**Retention, and the asymmetry that drives it.** `ingest/` had never had any. There was
an untracked `clean_grib_cache.py` whose idea -- never delete a GRIB a workspace still
symlinks from `wps/<SRC>/GRIBFILE.???` -- is right but is a trap as a *default*:
workspaces are never pruned either, so every old run pins its GRIBs forever and the
sweep reclaims steadily less until it reclaims nothing. It survives as `--keep-linked`.

The rewrite turns on a per-source `replaceable` flag, because the two halves of the
cache have opposite risk:

| | status | default |
|---|---|---|
| HRRR, HRRRA, HRRR_AK | archived on AWS, cheap to refetch | **swept** |
| NAM218/198/196/227 | **being retired -- this cache becomes the only copy** | **protected** |

Per JH, keep the NAM gribs of all varieties. A bare run therefore cannot touch them
however the cache is laid out; naming one still works for a dry run, with a warning, but
`--apply` on a protected source refuses unless `--force-irreplaceable` is also given
(exit 1). Dry runs stay unrestricted -- only the irreversible step gains friction.

Default sweep at 30 days is 6.9 GB of HRRR_AK. Naming NAM218 would show 2.5 TB and
NAM198 1.9 TB; **those were deliberately not applied.** The HRRRA cache ForeFire reads
has nothing older than 30 days, so it is unaffected.

**The two ingest trees are now one.** `cache_grib_files.py` fills
`/data/jhaley/wrfxpy/ingest/HRRR` while the FMDA cycler reads
`/data/jhaley/clean_wrfxpy/wrfxpy/ingest/HRRR` -- separate directories, so every GRIB was
being fetched twice (~56 GB/day duplicated). Per JH the second is now a **symlink** to
the first. The merge kept both sides: 7 GRIBs existed only in the cycler tree, plus the
`.size` sidecars the downloader uses for cache bookkeeping, and the 2 collisions were
byte-identical before removal. Verified that writes through the link land in the shared
tree, so future cycler downloads populate the one cache.

**A second, worse `retrieve_gribs` bug, in the *other* repo** (`e1444cc` on
`release-fmda-fixes`). The two repos carry diverged copies of this file; the clean one is
newer and had `elif grib_src_name == 'NAM' or 'NAM218':`, which parses as
`(name == 'NAM') or ('NAM218')`. A non-empty string is truthy, so the branch matched
**every name that reached it**: NAM227, CFSR_P/S, NARR, GFSA, GFSF_P/S, RAP and every
RRFS_* silently resolved to NAM218, downloaded NAM data under the requested name, and
printed NAM's Vtable as the one to use. Only HRRR_S and HRRR escaped, being earlier in
the chain. It matters now rather than eventually: **RRFS is a candidate NAM replacement**,
and every attempt to retrieve it would have quietly returned NAM218 right up to the point
NAM stops being published.

**Still to do for piece 1:** a cron entry for the sweep, and a decision on when NAM comes
out of `sats`. Per JH keep both for now and compare NAM- vs HRRR-driven forecasts as the
deprecation date approaches. Neither was done unasked.

### Piece 3: FMDA forecast mode runs, and ForeFire can read it (`0c88972`, `1f39148`)

`./hrrr_cycler.sh f CONUS_HRRR` had never been run. It works: it seeds from an
assimilated analysis at the cycle, then advances without DA, pinning to the 12z cycle.
9 leads in 25 minutes, no errors, ~90 s per lead. `forecast_length` set to 12 (via
`etc/fmda_cycler.json`, which is a **symlink** to `fmda_cycler_conus.json`).

The naming confirms the reader: `fcst_hour` counts from `cycle`, not `cycle_start`, so
valid time = the cycle in the name + NN. `read_fmda_hourly` now parses with `_FMDA_RE`
and orders candidates by nearest valid time, then analysis over forecast, then shortest
lead. Verified byte-identical output on the analysis-only path.

**Not done: repointing `fmda_dir` at CONUS_HRRR.** That changes every nightly forecast,
so it belongs with the piece 4 pSAF check, not ahead of it.

### Where the 1-hour fuel goes under the 0.02 floor, and why it matters

The floor is not an edge case. Over the forecast window the fraction of CONUS under it
climbs from 0.00% at 15:00 UTC to **2.18% at 23:00** -- ordinary afternoon drying.
Concentrated overwhelmingly in the interior West, at 23:00 UTC = late afternoon there:

| region | cells < 0.02 | % of box |
|---|---|---|
| Great Basin / Intermountain (36-45N, 122-110W) | 16682 | **14.75%** |
| N Rockies / Plains (44-49N, 115-100W) | 5071 | **7.00%** |
| S Texas / W Gulf (25-31N, 101-94W) | 294 | 0.56% |
| East of 95W | 754 | 0.10% |

Not a localized artifact -- a coherent regional signal exactly where and when fires run.
Eagle_Springs (45.7N, -107.7) sits in the N Rockies box. **JH confirms these regions
match the moisture maps the cycler pushes to the web server**, which is an independent
check on the whole chain: the reader, the grid indexing and the regional breakdown all
agree with a product produced by a different code path.

**Part of the field is negative, and that reframes the floor.** The 1-hour minimum is
**-0.0594**; 2991 cells (0.16%) are below zero. Fuel moisture cannot be negative, so
`valid_range` is not rejecting implausibly dry fuel -- it is rejecting **filter
undershoot**. Two checks say that is exactly what it is:

- The histogram decays smoothly through zero (340835 cells in 0.08-0.12, then 235758,
  122715, 44700, 27930, 9730, then 2811 negative, 168, 12). No spike, no second mode --
  the tail of a distribution whose mode is well above zero, continuing past it.
- Negatives sit where the fuel is genuinely driest: the neighbourhood mean around a
  negative cell is **0.0168** against **0.1411** domain-wide, eight times drier.

So clamping is the right treatment and the old substitution was the wrong one in the
worst possible place: at peak burning time, in the driest 15% of the Great Basin, it
replaced a near-zero fine-fuel moisture with roughly the 100-hour value and stopped any
fire there from spreading. Had forecast moisture been switched on before the clamp
landed, this would have fired routinely rather than occasionally.

### Piece 4, the only one left: the pSAF re-check

Changing the moisture source changes Md at every step of every fire, so this gets
measured, not assumed.

1. **Repoint `etc/forefire.json`**: `fmda_dir` -> `wksp_fmda/CONUS_HRRR`,
   `fmda_geo_file` -> `CONUS_HRRR-geo.nc`. Both exist and are current. This is the one
   change; nothing else moves, so the comparison isolates the moisture source.
2. **Re-run and score** with `ff_batch_ratio.py`, passing the batch script so the window
   is bounded by its mtime (§6 -- a bare number can be captured by a filename).
3. **Compare against the clean baseline in §2**: geo mean 1.56, median 1.53, 25th
   percentile 1.11, n=58.
4. **Decide**: keep 0.337 if the 25th percentile stays within ~10% of 1.11 -- that is
   inside the ~3% per-fire reproducibility noise (§4) plus sampling noise on 58 fires.
   Otherwise `pSAF_new = 0.337 / p25^(1/1.54)`.

Per JH the RTMA/HRRR difference in resulting spread is small, so the constant will
probably survive. **Watch the per-fire scatter rather than the mean anyway**: §7 is the
cautionary tale, where a bug that moved Eagle_Springs from unscorable to 0.84 moved the
58-fire aggregate from 1.55 to 1.56.

Once forecast FMDA is in the loop, a *forecast-driven* run should be scored separately
from an analysis-driven hindcast -- and with `--issue-utc` set (§8, piece 2), or the
score silently credits cycles issued after the fire started.

## 9. Open items

Items 3, 4, 5 and 6 of the 09-24 handoff's §9 are untouched and still open: the fuel 2/9
transposition upstream, the false "unknown moisture source" log message, scoring against
observations rather than only WRF-SFIRE, and the fires with no overlapping valid time --
now **two**, Nethery and STEEL_PASS: Eagle_Springs turned out to be the moisture bug in
§7, not a timing problem, and scores 0.89 once it can burn.

Item 1 is closed. Item 2 — per-fire pSAF from NGFS detection fields — is the next real
move, now with a concrete target (§3).

## 10. NEXT SESSION

**Piece 4 of §8: the pSAF re-check.** Repoint `fmda_dir` at `CONUS_HRRR`, re-run, score,
compare against geo mean 1.56 / median 1.53 / 25th pct 1.11 (n=58). Everything it needs
is in place -- pieces 1, 2 and 3 are done and committed.

Two loose ends from piece 1, both left deliberately for JH:

- a cron entry for `clean_grib_cache.py` (nothing prunes `ingest/` yet);
- when NAM comes out of `sats` in `cache_grib_files.py`. Per JH, keep both for now and
  run NAM-vs-HRRR forecast comparisons as the deprecation date approaches. **Scope that
  to WRF-SFIRE**: `resolve_grib_source` takes either, but WindNinja *silently garbles*
  NAM218 -- exits 0, plausible valid time, speeds of 25.8 to 45,632 m/s -- so any
  ForeFire-side number from NAM winds would be garbage that does not announce itself.

Longer term, the per-fire pSAF question (§9 item 2) is the larger scientific one, but
weigh §9 item 5 first. The whole calibration is model-vs-model and WRF-SFIRE is not
truth; on Dome it ran 0.68x an IR perimeter while ForeFire ran 6.0x. A **global**
constant fitted against WRF-SFIRE inherits that bias on average, which is easy to revise
later. A **per-fire** correction fitted the same way would bake it in fire by fire, which
is not. If there are enough usable IR perimeters to fit against instead, that ordering is
worth the delay.

Commits remain unpushed on `james_ngfs` and `release-fmda-fixes`; JH pushes them manually
(SSH password).
