# ForeFire session handoff — 2026-09-21

Branch `james_ngfs`. **Four commits this session**, all local:

```
2360f33  Select HRRR cycles that actually reach f48, and step back on retry
54aa288  Record the 2026-09-18 ForeFire session: wind is not the shape gap
89151ac  Build ForeFire input netcdfs from GRIB winds, without WRF-SFIRE
0a79f15  Anchor every ForeFire step's clock to the ignition step
```

Plus one deliberate change outside git: `cron_forefire.sh` (§4). Follows
`SESSION_HANDOFF_forefire_2026-09-18.md`, whose §11 set the goal — run ForeFire
from an ignition point, an ignition time and GRIB files alone. **That is built and
verified.**

**Headline:** the GRIB-driven path works end to end and reproduces the coupled-WRF
result on a real fire — Union_400257, ForeFire from HRRR+WindNinja gave a 10-member
ensemble mean of **176.1 ha against WRF-SFIRE's 178.1 ha** at a common valid time,
a 1.1% difference. Getting there turned up **two bugs that were corrupting every
ForeFire ensemble ever run here**, neither related to winds: output times were
mislabelled (§3), and more than half of every ensemble was being silently emptied
by overlapping cron runs (§4).

---

## 1. `forefire_grib.py` — ForeFire inputs without WRF-SFIRE

`src/ingest/forefire_grib.py`, commit `89151ac`. Builds the same netcdfs
`make_FF_nc` builds, from an ignition alone:

    netcdf var    was                        now
    fuel          wrfinput_d01:NFUEL_CAT     LFmosaic_FBFM13_CONUS.tif
    altitude      wrfinput_d01:ZSF           LF2020_Elev_220_CONUS
    windU/windV   wrfout:UF/VF               WindNinja downscaling a raw HRRR wrfprs
    domain        WRF's LCC metadata         synthesised for the grid around the ignition

The fuel raster is the one `etc/vtables/geo_vars.json` already points geogrid at, so
a GRIB-driven run and a geogrid run see identical fuels.

```bash
PYTHONPATH=src:src/ingest python src/ingest/forefire_grib.py \
  <lat> <lon> <ign_utc> <run_dir> --grid-code <GC> \
  --steps 8 --domain-km 30 --fire-dx 30 [--wn-mesh 50] [--windrf ...]
```

**It needs no change to `forefire.py`.** Files are named exactly as `ff_nc_name`
expects, and `make_script_set` skips its own builder for any step whose netcdf
already exists, so `make_script_set`, `run_timing_table` and the `.ff` templates
drive these unmodified. The 09-18 §11 worry about `make_timing_table` was
misplaced — the ensemble machinery needed nothing.

**Verified, not assumed:** ForeFire ran on a generated netcdf and produced a
perimeter centred on the ignition, which confirms `startFire[lonlat=...]` resolves
through the synthesised `BBoxWSEN` — the part of the domain maths most likely to be
wrong. The structure was diffed against a real `make_FF_nc` output from the Dry
River arm A run: identical dimension names, variable names, dtypes, dimension
tuples, `type` attributes and `domain` attribute set.

**Cost**, 1200x1200 at 25 m with a 50 m WindNinja mesh: **~2 min per step**,
dominated by the warps and the netcdf write, not the solve. A 10 km/100 m test was
~11 s per step.

## 2. Decisions inside it that are choices, not derivations

Do not "fix" these without knowing they were deliberate.

- **`windrf` is a parameter, not a constant.** `--windrf` takes all 14
  per-category factors, defaulting to `etc/nlists/default.fire`. Per JH it should
  be set **per fire**, from observations and the general fire location. 09-18 §5
  is why it cannot simply be omitted.
- **Category 14's `windrf` of `1e-7` is a sentinel, not a reduction.** Applying it
  would punch near-zero wind holes into the field, so non-burnable cells take the
  modal burnable factor (0.36) to keep the field continuous.
- **Fuel remapping follows `geo/var_wisdom.py`**, including the `nearest` fill for
  codes 0, 91 (urban) and 93 (agriculture), so a fire is not stopped by a road or a
  field edge. Verified: 91 -> 2, 93 -> 1, 99 and -9999 -> 14.
- **AAIGrid writes row 0 at the north edge; ForeFire is indexed south-to-north**,
  so rasters are flipped on read. Getting this wrong mirrors the fire north-south
  and **is invisible on a symmetric domain**.
- **Naive/aware, again.** An ignition string carrying `Z` parses tz-aware while
  GRIB valid times are built naive from filenames. Normalised through one
  `as_naive_utc` helper at the boundaries rather than a guard per comparison —
  the same instruction the FMDA handoff gives for its own copy of this.
- **The GRIB search starts an hour before the ignition's own hour**, so the
  ignition falls *inside* a step. Starting at the ignition drops the ignition step
  and the fire silently begins at the next whole hour.

## 3. Every ForeFire output time was mislabelled

**This affects every run in the archive, not just GRIB-driven ones.**

ForeFire runs **one continuous clock across a restart chain, anchored at the
ignition**. The state files say so:

    FireDomain[sw=(0,0,0);ne=(30000,30000,0);t=10806]
        FireFront[id=2;domain=0;t=1640]        <- 1640 is the ignition offset

and `include[]` re-seats the front at the `FireDomain` time. So `START_STEP`,
`END_STEP` and `END_TIME`, all `ign_seconds` from that anchor, were **correct**.

What was wrong is the date handed to `loadData`. Each step got **its own**
timestamp, and that date is what ForeFire stamps outputs with — `valid_at` comes
out as *(that date + t)*. Every step after the ignition was mislabelled by its own
offset from the ignition step, growing one step at a time. On Union, step 05
printed t=18000 — truly 01:00Z — stamped 05:00Z; step 08 printed t=28800, truly
04:00Z, stamped 11:00Z. Proof that `valid_at = own_nc + t`:

    step 04  nc 20:00  t=1800   -> 20:30
    step 10  nc 23:00  t=12600  -> 02:30
    step 24  nc 06:00  t=37800  -> 16:30

Fixed in `0a79f15`: every step's `loadData` gets the **ignition step's** timestamp.
The simulation was always right; only the labels were wrong — which made a
time-lagged ensemble look like its members ended at different times when they did
not.

**A dead end worth not repeating:** the first attempt made `END_TIME` per-step,
computed from timestamps. It fixed the members where `END_TIME > END_STEP` and
broke the rest, because a smaller `END_TIME` makes the second `goTo` run backwards
and become a no-op, leaving the member stopped at `END_STEP`. It was treating the
symptom. Reverted.

## 4. Overlapping cron runs were emptying most of every ensemble

`cron_forefire.sh` ran every 20 minutes (`7,27,47`) with **no `flock`**, and
`run_forecasts` **deletes every file in the shared run directory** on entry:

```python
for stale in glob.glob(f"{forefire_dir}/*"):
    if os.path.isfile(stale) or os.path.islink(stale):
        os.remove(stale)
```

`run_dir` is one fixed path for every fire. A 47-step fire takes longer than 20
minutes, so a second invocation routinely wiped the first one's netcdfs mid-step.
ForeFire then fails with `wrong input file, check your settings...`, writes a
restart holding `FireDomain t=0` with no nodes, and **every later step in the chain
inherits no fire**. On Union the run-log mtimes show it exactly:

    09-20 15:47   21 steps ran   <- last good restart 37800
    09-20 16:07   14 steps ran   <- step 39600 fails, chain dies here
    09-20 16:27   12 steps ran

That left **26 of 47 members empty**. Step 25 re-run in isolation completes fine,
so nothing was wrong with the data.

**Fixed in `cron_forefire.sh`** (untracked, matching every other cron script here):
`flock -n /tmp/cron_forefire.lock`. Two details that matter — the python is no
longer backgrounded, because with `&` the shell exited and flock released the lock
immediately, making it a no-op; and `forefire.log` is truncated *inside* the lock so
a skipped tick cannot destroy the log of the run still writing it. Skips are
recorded in `forefire_skipped.log`.

After the fix, Union re-ran **47 of 47 steps in 1.9 minutes**. The original took
over 20 minutes for 21 — most of that was three invocations thrashing one directory.

## 5. Union_400257 — the result

North-east Arkansas, ignition 2026-09-20 20:27:20Z at 33.3499, -92.64239. 30 km at
25 m. WRF-SFIRE crashed on process 19 with the last wrfout at 2026-09-21 19:30, per
JH; that sets the timing table's end and is unrelated to §4.

**At a common valid time of 2026-09-21T06:00Z:**

    WRF-SFIRE chained track        178.1 ha
    HRRR+WindNinja, 10 members     169.9 - 193.9 ha, mean 176.1

**1.1% apart, with the WRF value inside the ensemble range.** Driving ForeFire from
the cached hourly HRRRA files as if they were forecast GRIBs reproduces the
coupled-WRF result on this fire. Full ensembles, each internally consistent after
§3 and §4:

    WRF-SFIRE (30-min steps)  47/47 members  all 19:30Z   869-1107 ha  spread 1.27
    HRRR+WindNinja (hourly)   10/10 members  all 06:00Z   170-194 ha   spread 1.14

They end at different instants because the f03 cache could not reach 19:30Z — that
was 3.5 h in the future at run time (§8). The ensembles are therefore **not**
comparable to each other on area or spread; only the 06:00Z figures above are.

**A number to distrust if you find it elsewhere:** an earlier reading of this fire
gave "21 members, 190-215 ha, spread 1.13". That was the §4-truncated run, only
10.5 h in. It is wrong.

## 6. Red Bank — IN FLIGHT AT HANDOFF

`Red_Bank_2026-09-15_17_00_00_57A034B6-C5DA-4780-BA17-014940D059CC`, north Texas,
ignition 2026-09-15 18:27:21Z at 32.56616, -97.90548. Same geometry as Union
(30 km, 25 m, 1200x1200), 61 wrfouts spanning 15:00 09-15 to 21:00 09-16.

Chosen because it is **entirely historical**: 30 GRIBs, no gaps, covering the whole
WRF window, so unlike Union the HRRR ensemble can run to the **same end time** and
give a true like-for-like comparison of both areas and spreads.

Two jobs were running when this was written:

- HRRR netcdf build, 28 steps into `/home/jhaley/forefire/tests/ffwksp_redbank_hrrr`
  (`--wn-mesh 50`, justified: DEM std **31.5 m**, next to Dry River's 26.7 m which
  09-18 §4 measured as converged at 50 m).
- WRF ensemble re-run with `overwrite=True` under the cron lock, to relabel it
  with §3. It is **far slower per step than Union** — 33 of ~56 steps in about 6
  minutes, against Union's whole 47-step run in 1.9 minutes — but it is progressing,
  not stuck. Red Bank's fire is much larger, and each step's `goTo[t=END_TIME]` runs
  to 97200 s here against Union's 84600. Watch for the 240 s `apptainer.timeout_s`
  on the early steps, whose second `goTo` is the longest single run in the ensemble;
  a step killed there would write no final geojson.

**A background job does not survive the session ending** — the same warning as
09-15 §F. Check both before assuming results exist.

## 7. WindNinja usage — two traps to add to 09-18 §2

Both cost a run and neither is in that section:

1. **`--time_zone` is mandatory** whenever `--forecast_filename` is used:
   `Exception caught: Option 'forecast_filename' requires option 'time_zone'.`
   The timestamp WindNinja puts in its output filename is in that zone, **not UTC**.
2. **The env must be `conda activate`d, not just put on PATH.** Its activation
   scripts set `PROJ_DATA`/`GDAL_DATA`; without them the DEM `gdalwarp` dies with
   `PROJ: proj_create_from_database: Open of .../share/proj failed` even though
   `proj.db` is right there. This bites the DEM prep, not WindNinja, so it looks
   unrelated.

Confirmed again this session: WindNinja ingests the raw cached `wrfprs` files
directly and resolves valid time correctly (t00z f03 -> 03Z).

## 8. The f03 ceiling

The FMDA HRRR cycler's cache (`clean_wrfxpy/ingest/HRRRA`, 71 GB, 8 days) is
**f03 only, hourly**. That makes this path a hindcast/nowcast: it cannot forecast
past roughly now. Union hit exactly this — its ensemble stops at 06:00Z because
19:30Z had not happened yet.

`find_gribs` / `_grib_valid_time` in `forefire_grib.py` is the **single seam** to
widen when HRRR forecast cycles are cached, which JH says is planned. The
`cycle_hours = 6` selection committed this session as `2360f33` is already the
correct selector for those, since only the 00/06/12/18Z cycles reach f48.

## 9. Open items

1. **Nothing is pushed.** Four commits this session; `james_ngfs` is local only.
2. **Red Bank is unfinished** (§6), and the WRF side of it may be wedged.
3. **Every archived ForeFire run has wrong `valid_at` stamps** (§3) and most have
   mostly-empty ensembles (§4). Geometry is fine; labels and completeness are not.
   Re-running needs `overwrite=True` and the lock. **How far back to go is
   undecided** — only Union has been redone.
4. **Workspaces touched by cron between the §3 fix and their re-run are mixed** —
   `.ff` scripts carry the new anchor while geojsons are stale. The 09:47 cron
   regenerated 192 files in the Union workspace this way before it was re-run.
5. **Why members stop before `END_TIME`.** Union's WRF members stopped at t=37800
   with `END_TIME=84600` before the fixes; afterwards they reached it. Whether any
   residual early stopping remains is unverified on other fires.
6. **`forefire.log` has no rotation** and now also receives skip lines.
7. `windrf` per-fire calibration is still an open question (§2, 09-18 §5).
8. **Southern Great Plains agreement**, per JH: WRF-SFIRE and ForeFire agreed
   surprisingly well there over the preceding few days, worth investigating. It
   sits against the Dry River result where the two failed in *opposite* directions
   (09-18 §7). Nothing scored yet.
9. 09-18 §10 items 3 (`.ff` templates hardcode `fuelstrans.csv`) and 5 (Ranger Road
   three-cell experiment) are untouched.

## 10. NEXT SESSION

1. **Finish Red Bank** (§6) and get the first like-for-like ensemble-vs-ensemble
   comparison at one valid time. That is the test Union could not provide.
2. **Decide the re-run scope** for §9 item 3. Every archived ensemble is affected.
3. **Wire NGFS ignitions in.** `forefire_grib.py` takes lat/lon and time as
   arguments by design; per JH the ignition should come from the same logic that
   drives the NGFS forecasts, whose points he considers better than NIFC's.
4. **Set expectations from §5**: on this evidence a GRIB-driven forecast lands on
   top of the coupled one, so the limit is still ForeFire's shape deficit
   (09-18 §7), not the winds.
