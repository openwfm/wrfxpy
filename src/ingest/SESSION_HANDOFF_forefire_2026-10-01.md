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

> **Corrected 2026-10-05.** On 31 fires at pSAF 0.337 the 8 h measurement is 1.40 / 1.36
> / 58%, against §5d's predicted 1.42 / 1.30 / 78%. The centre was predicted almost
> exactly; only the spread was over-optimistic. The 11-fire reading below was too small.
> See `SESSION_HANDOFF_forefire_2026-10-05.md` §5p.

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

### A second front read as a hole (2026-10-02, per JH)

JH found two contours at steps 27 and 28 of Hatch_Grande in both the cron and forward
kml, "almost like a second front was ignited partially within the first", and asked
whether subtracting fronts inside the outer perimeter explained the shrink. It did.

| step | rings | |
|---|---|---|
| 26 | 3 | 1493.3 BURNING, 5.95 + 3.56 islands inside |
| 27 | 5 | **1601.1 BURNING** (grew), **147.7 "island"** overlapping the outer ring, 57.5 BURNING, 5.95 + 2.03 islands |

The outer ring grew by 108 ha and the measured area **fell** 1483.8 -> 1457.8, because
the second front's 147.7 ha boundary was differenced out. It cannot be a hole: most of
that ground was already enclosed at step 26, and **burned ground does not become
unburned**.

### The discriminator is containment, not area

A hole lies entirely inside the ring it punctures. `read_ff_geojson` decides containment
from a representative point, so a ring that merely overlaps is taken for a hole.
Measured on both artefacts and the genuine islands beside them:

| ring | % of its area outside the burn |
|---|---|
| genuine islands (5.95, 5.39, 2.03 ha) | **0.000** |
| Hatch second front (147.7 ha) | 8.37 |
| SAN_FRANCISCO thin band (1607.4 ha) | 8.24 |

So a nested ring is differenced out only when enclosed, and unioned in otherwise. **This
replaced the area-ratio threshold below**, which caught the thin-band case and not the
second-front case; containment is an invariant rather than a tuned constant and covers
both. Shrink events across 63 runs fell from 2 to 1, and the survivor
(Lauderdale, -1.41%) carries no ring artefact and is the restart re-mesh loss the
warning threshold exists for. Final areas move where the artefact was:
SAN_FRANCISCO_2 cron 1621.8 -> 1754.3 ha, Hatch_Grande cron 1387.5 -> 1488.8. The
forward batch's scored times precede both artefacts, so §5i is unchanged again.

### Superseded: the area-ratio threshold

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

## 5l. RAP and RRFS ported, and verified through real.exe (2026-10-02)

Per JH, brought across from `origin/release` into `james_ngfs` as a **port, not a merge**
-- merge base 2023-10-09, 97 files, ~19 conflicts for two modules. Commits `f018aab`,
`c0918e0`, `aa0358c`, `713bb66`, `f7e4221`.

### What is in

| source | job name | products | grid |
|---|---|---|---|
| RAP | `RAP` | wrfprs | 13 km CONUS/NA |
| RRFS CONUS | `RRFS` | prslev + 2dfld | 3 km |
| RRFS NA | `RRFSNA` | prslev + 2dfld | 13 km |
| RRFS Alaska | `RRFSAK` | prslev + 2dfld | 3 km |
| RRFS Hawaii | `RRFSHI` | prslev + 2dfld | 2.5 km |
| RRFS Puerto Rico | `RRFSPR` | prslev + 2dfld | 2.5 km |

`_P` and `_S` are **halves of one input**, not alternatives: metgrid's `fg_name` lists
both, so fetching only `_P` passes ungrib and fails in metgrid. AK was added because
`HRRR_AK` was the only other Alaska source; PR because those fires currently run off
GFS and 2.5 km convection-allowing output should beat it.

### Three defects found, all the same defect

Every one was a cycle paired with a lead that cycle does not serve, and each would have
404'd silently on most cycles:

| | written as | measured | fixed to |
|---|---|---|---|
| RAP | `cycle_hours=1`, 51 h | only 03/09/15/21 reach f51; standard cycles stop at f21 | `cycle_hours=6`, `cycle_offset=3` |
| RRFS CONUS | `cycle_hours=3`, 84 h | only 00/06/12/18 reach f084; t15z stops at f018 | `cycle_hours=6` |
| RRFS AK/HI/PR | -- | same rule, checked per domain before writing the classes | `cycle_hours=6` |

The first fix was itself wrong and a live run caught it: `cycle_hours=6` selects
multiples of 6, i.e. **00/06/12/18** -- precisely RAP's non-extended cycles. The
selector had no way to express a cycle grid that is offset from midnight, so
`grib_forecast.py` gained **`cycle_offset`** (default 0, every existing source
untouched) and subtracts a timedelta rather than replacing `.hour`, so 01z rolls back
to the previous day's 21z instead of going negative.

### The one that actually blocked JH

`forecast.py` carries a **second** list of grib sources at job validation (line ~1495),
separate from `resolve_grib_source`. Every source above was fully implemented and still
rejected before it could run. I had grepped for `available_grib_sources` -- the name it
has on the release branch -- found nothing, and concluded no such list existed. It is
now a named `valid_grib_sources` that the error message formats from, cross-checked to
carry the same 17 entries as the dispatch.

### Verified end to end on the PIUTE fire

Run in `/data/jhaley/new_wrfxpy/wrfxpy`, which shares this tree through
`/home/jhaley/work -> /data/jhaley`, so its `forecast.py` and `retrieve_gribs.py`
symlinks already resolved here. Both variants reached `wrf.exe`:

| | ungrib -> metgrid | levels declared | **levels in met_em** | real.exe |
|---|---|---|---|---|
| RRFS | 07:22 -> 07:49 | 46 | **46** | pass, wrfinput 116.4 MB |
| RAP | 08:15 -> 08:37 | 41 | **41** | pass, wrfinput 116.4 MB |

**So the HRRR-copied vtables work for both**, which was the port's stated risk, and
both level counts were right -- `RAP.py`'s docstring claiming 38 was the thing that was
wrong, now corrected. RAP reaches WRF in 23 minutes against RRFS's 41, being one
product per lead rather than two.

**What this does NOT establish** is that the fields are physically sensible. A wrong
vtable entry can map a plausible-looking field to the wrong variable and still satisfy
metgrid and real. Only comparing the WRF output against the HRRR or NAM run of the same
fire shows that, which is what the PIUTE variants are staged for.

### Two notes for whoever runs these next

- **Volume.** RAP ~230 MB per lead; RRFS CONUS ~900 MB for both products, so an 84 h
  cycle is ~75 GB. Hawaii and Puerto Rico are ~30 MB per lead -- the domains are tiny.
- **The staged RRFS job had the RAP `grid_code`**, so both variants would have written
  to one workspace. Worth checking when cloning a job to change only `grib_source`.

### The physical check: RAWS says NAM was wrong and the new sources are right

Passing `real.exe` showed only that the vtables are structurally sound. The open worry
was that a mis-mapped field could look plausible and still satisfy metgrid and real.
The PIUTE runs answer it.

**Burned area**, same fire, same 26 h window, TIGN_G:

| | final ha | vs NAM |
|---|---|---|
| NAM | 62.9 | -- |
| RAP | 82.8 | **+32%** |
| RRFS | 82.2 | **+31%** |

RAP and RRFS land **within 0.7% of each other** on different native grids (13 km vs
3 km) with different file structures, which already argues against a field being
mis-read.

**Why**, from the near-surface fields: NAM is far wetter. Checked against the 14 RAWS
within 40 miles and 1000 m of the fire's ~2300 m elevation (Synoptic, network 2):

| valid UTC | local | RAWS | NAM | RAP | RRFS |
|---|---|---|---|---|---|
| 10-02 00:00 | 17:00 | 18% | 42% (+24) | 23% (+6) | 27% (+10) |
| 10-02 03:00 | 20:00 | 26% | 53% (+28) | 30% (+4) | 36% (+10) |
| 10-02 06:00 | 23:00 | 28% | 58% (+30) | 24% (-4) | 34% (+7) |
| 10-02 09:00 | 02:00 | 28% | 61% (+32) | 23% (-6) | 32% (+3) |
| 10-02 12:00 | 05:00 | 26% | 61% (+36) | 30% (+4) | 34% (+9) |
| 10-02 15:00 | 08:00 | 24% | 52% (+28) | 30% (+6) | 30% (+6) |

| | mean bias | RMSE |
|---|---|---|
| NAM | **+29.6 pts** | 29.9 |
| RAP | **+1.8 pts** | 5.1 |
| RRFS | +7.5 pts | 7.9 |

**NAM sits outside the entire station spread at every hour** -- above the wettest of the
14, not merely above the median -- holding RH near 60% overnight in the Sierra in
October while every nearby gauge read 20-28%. RAP is effectively unbiased; RRFS is
modestly wet.

**So the 32% larger fire is in the direction the observations support**, and a
mis-mapped humidity field would not land within 2 points of 14 independent gauges. If
anything a 30-point RH error should cost more than 32% in area, so the burned-area gap
probably *understates* how much NAM was damping this fire.

Three limits worth keeping with the result:

- RAP and RRFS share the same HRRR-derived vtable, so this validates **that vtable**,
  not each source independently.
- The 14 stations span 40 miles and 1500 m of elevation; the median carries real spatial
  uncertainty, though NAM's error is far outside it.
- The fire is PDT and these runs were scored the morning of 10-02, so the observations
  cover **10-01 afternoon through 10-02 morning**. The local afternoon of 10-02 was not
  yet observed; `10-02 00:00Z` is the one true afternoon point.

### ForeFire on the same three runs: an exceptional fire, and not for the reason expected

With three WRF-SFIRE runs of one fire differing only in driving weather, ForeFire was
run on each by the analysis path (`run_forecasts`, chained to each run's own wrfout
winds, pSAF 0.337). NAM's was already on disk from the 10-02 cron, so only RAP and RRFS
needed running. All three are clean -- 50 perimeters each, **zero shrink events**, which
is also the first independent exercise of §5j's monotonicity guard on fresh output.

| source | ForeFire ha | WRF-SFIRE ha | ratio |
|---|---|---|---|
| NAM | 183.9 | 62.9 | **2.92** |
| RAP | 167.6 | 82.8 | **2.02** |
| RRFS | 290.5 | 82.2 | **3.53** |

For scale, the 25-fire cron batch the same morning had a geometric mean of 1.54 and a
max of 4.55, so PIUTE sits at the top of the range under every driver.

**ForeFire does not amplify the driver spread -- it reorders it.**

| | NAM | RAP | RRFS |
|---|---|---|---|
| WRF-SFIRE | 62.9 | 82.8 (**+32%**) | 82.2 (**+31%**) |
| ForeFire | 183.9 | 167.6 (**-9%**) | 290.5 (**+58%**) |

In WRF-SFIRE, RAP and RRFS agree to within 1% and both beat NAM by about a third. Given
the same three winds, ForeFire puts NAM in the middle and splits RAP and RRFS by 73
points. **An earlier reading of this as "ForeFire roughly doubles the sensitivity to the
driver" was made on the first two points and is wrong**; the third removes the pattern.

### It is not the fuel map, which is what JH's granite observation predicted testing

Per JH, PIUTE "mainly burns in forested areas surrounded by granite" -- an exceptional
case. If ForeFire were spreading across rock that WRF-SFIRE respects, that alone would
explain everything. It is not:

| | non-burnable (cat 14) | cat 5 | cat 8 | cat 10 |
|---|---|---|---|---|
| ForeFire | **37.07%** | 26.24% | 22.39% | 11.93% |
| WRF-SFIRE | **37.26%** | 25.59% | 22.78% | 12.14% |

The same map to within 0.2 points. **37% of the domain is bare rock in both models**, so
ForeFire is over-spreading *within* the forest pockets rather than crossing the granite.
That also weakens the terrain explanation, since both see identical elevation and fuel
on the same 25 m mesh -- for the record the fire mesh is far steeper than d01 suggests
(mean slope 16.4°, 90th percentile 30.7°, 32% of cells over 20°, against 5.1° mean at
1 km), but WRF-SFIRE integrates that same terrain.

Moisture is ruled out too, and points the wrong way as it did in §5h and §5j: FMDA gave
ForeFire 1-h 0.020-0.059 while WRF-SFIRE's own `FMC_GC` ran 0.006-0.049, mean 0.030.
**WRF-SFIRE burns a third as much on drier fine fuel.**

### What is left is connectivity

With 37% of the domain blocked in scattered granite, the burnable ground is a maze of
narrow necks, and spread becomes a question of whether the front threads a gap -- which
is exactly where `perimeterResolution 25` and ForeFire's front tracking differ most from
WRF-SFIRE's level set on identical cells. One neck taken that should not be compounds
for the rest of the run. It also explains the volatility across drivers: with a third of
the domain blocked, small wind-direction differences decide **which** necks open, so
±1 m/s reorders the outcome instead of scaling it. `ff_domain_connectivity.py` has been
open since 09-28 for this, and PIUTE is the case that justifies running it.

**At 8 hours the spread is much tighter** -- 2.93 / 2.19 / 2.07, a range of 0.86 against
1.51 at the endpoint -- so the low-latency product is less exposed to this than the
full-duration one.

### Scope limit, per JH

**PIUTE is exceptional and must not be allowed to pull a recalibration.** 37%
non-burnable with the burnable fraction in disconnected pockets is nothing like the
batch pSAF 0.337 was fitted on, which is dominated by flat, continuous prescribed burns
in the southeast. It is a good diagnostic case and a bad calibration case -- the same
judgement §5f recorded for Pulaski, for a different reason.

### `ff_domain_connectivity.py` confirms it, and PIUTE is at the ceiling

Run against all 391 retained domains, and against the RAP/RRFS pair separately.

| | PIUTE | min | 10th | median | 90th | max |
|---|---|---|---|---|---|---|
| burnable | **62.4%** | 0.347 | 0.880 | 0.982 | 0.998 | 1.000 |
| reach | **79.8%** | 0.000 | 0.787 | 0.999 | 1.000 | 1.000 |
| fragility | **1.00** | 0.008 | 0.015 | 0.053 | 0.350 | **1.000** |

**Fragility 1.00 is the maximum the measure can take**: eroding the burnable mask by two
cells collapses the reachable area completely, against a batch median of 0.053. PIUTE is
not in the tail, it is at the ceiling, and it trips all three tests when only 73 of 391
domains trip any. RAP and RRFS give 62.9% / 81.3% / 1.00, identical as expected on a
shared fuel map.

That is the mechanism. Where removing two cells of fuel severs the domain, spread is
decided by whether the front threads a neck rather than how fast it travels along one --
exactly where front tracking at `perimeterResolution 25` and a level set on identical
cells diverge -- and it explains the reordering across drivers, since a small
wind-direction change decides **which** necks open.

**Independent corroboration from LAC-352011.** It is flagged too, at 56.4% burnable,
**0.6% reach**, fragility 1.00 -- almost none of its burnable ground connects to the
ignition. That is the fire whose `fwd/cron` ratio of 0.21 survived the `track()` fix at
0.23 and which JH identified as odd terrain with ponds near Edwards AFB. Three
independent routes -- JH's read of the terrain, the area outlier, and this measure --
reach the same conclusion without reference to each other.

**This sharpens the per-fire pSAF item** (09-24 handoff §9 item 2). **73 of 391 domains
should be excluded from any rate calibration**, because their area response to a
parameter change is binary rather than proportional: the answer is "did it get through
the gap", and no exponent absorbs that. Fitting on them adds noise that looks like
unexplained variance.

## 6. What today added to the open list

The full list stays in `SESSION_HANDOFF_forefire_2026-09-29.md` §6 and §7, which this
file does not duplicate. Today added three items to it:

1. ~~**Add a monotonicity check to both scorers.**~~ **DONE.** `ff_growth.shrink_events`
   finds every step where the burned area falls, and both `ff_batch_ratio.py` and
   `ff_fwd_score.py` now call it: below 10% the run is reported and still scored, at or
   above 10% it is disqualified and named. Regression-checked against the pre-fix
   measurement of SAN_FRANCISCO_2, which it flags at **-43% on step 19** and would have
   skipped; the fixed measurement raises nothing. The thresholds are measured rather
   than picked -- over 57 runs the only real losses left after the ring fix are
   Hatch_Grande at 1.75% and Lauderdale at 1.41%, both from a chained restart re-meshing
   the perimeter it was handed, so 1% warns and 10% disqualifies.
2. **Score time-to-area beside area-at-time** (§5j). LAC-352011's cron run reaches
   WRF-SFIRE's final area 13 h early while tracing what JH reports is nearly its final
   outline. That is a rate error on an agreed path, and the current metric cannot
   distinguish it from a fire spreading the wrong way.
3. ~~**The RAP/RRFS decision**~~ -- **ported and validated** (§5l). What remains:
   re-check against the 10-02 **local afternoon** RAWS once those observations land,
   since the runs were scored before that window; chase RRFS's **+7.5 point wet bias**,
   which is small but one-signed across all six hours; and validate the vtables for the
   NA, AK, HI and PR domains, which have their own vertical structures and have only
   been checked as far as filename and cycle correctness.

4. ~~**Run `ff_domain_connectivity.py` on PIUTE**~~ **DONE (§5l).** Flagged on all three
   measures, fragility at the ceiling.
5. **Exclude the 73 connectivity-flagged domains before any per-fire pSAF fit**
   (§5l, 09-24 handoff §9 item 2). Their area response to a parameter change is binary,
   not proportional, so they add noise a rate exponent cannot absorb. This is now a
   precondition for that fit rather than a caveat on it.

## 7. NEXT SESSION

**Try the RRFS_NA forecast of PIUTE.** Per JH. **DONE 2026-10-05 and it cannot be run** --
RRFS NA is on a rotated lat-lon grid that WPS 4.4's ungrib does not support, so it fails
before metgrid. See `SESSION_HANDOFF_forefire_2026-10-05.md` §5m. The original reasoning
follows. The CONUS 3 km run is done (§5l); `RRFSNA`
is the 13 km North America product, so it tests the one RRFS domain whose vertical
structure has not been through metgrid -- `num_metgrid_levels` is still the module's
guess for NA, as it is for AK, HI and PR, and a wrong value fails in `real.exe` rather
than in ingest. It also gives a resolution pair on an identical fire: 3 km against 13 km
with everything else held, which is the cleanest read available on how much the extra
resolution is worth. Note `RRFS_NA` runs `period_hours = 3` while the bucket carries
hourly files, so its boundary updates are 3-hourly by choice, not by availability.

The job is `jobs/PIUTE_..._RRFS.json` in `/data/jhaley/new_wrfxpy/wrfxpy` with
`grib_source` changed to `RRFSNA` -- and **change `grid_code` too**, or it overwrites the
CONUS workspace, which is the trap the staged RRFS job already fell into once.

Still open from 09-29 and untouched today: scoring the ensemble members rather than
only the chained track, isolating the mesh at 50 m with analysis moisture, the roads
layer for fire breaks, and re-running the moisture-bias test on a larger batch.
