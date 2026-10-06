# ForeFire session handoff — 2026-10-05

Continues `SESSION_HANDOFF_forefire_2026-10-01.md`. Section numbering carries on from
there so existing references keep resolving.

## 5m. RRFS_NA cannot be used with this WPS: it is on a rotated lat-lon grid

The 10-01 handoff's NEXT SESSION item was the RRFS_NA forecast of PIUTE, to test the one
RRFS domain whose vertical structure had never been through metgrid. **It never reached
that test.** It died in ungrib:

```
unknown model from NCEP 134
GRIB2 Unknown Projection: 1
ungrib understands projections 0, 10, 20, 30, and 40
```

`num_metgrid_levels` was never the issue. The NA product is on **GRIB2 grid definition
template 1, rotated lat-lon**, which WPS 4.4's ungrib does not implement:

| domain | grid | template | ungrib |
|---|---|---|---|
| RRFS CONUS 3 km | Lambert Conformal (1799 x 1059) | 30 | **works** (§5l) |
| **RRFS NA 13 km** | **rotated lat-lon (1127 x 683)** | **1** | **unsupported** |
| RRFS AK 3 km | polar stereographic (1649 x 1105) | 20 | should work |
| RRFS HI 2.5 km | Mercator (321 x 225) | 10 | should work |
| RRFS PR 2.5 km | Mercator (544 x 310) | 10 | should work |

**NA is the only one of the five on an unsupported grid**, checked with `wgrib2 -grid`
against the files each class had already downloaded. That matters for the Puerto Rico
plan in §5l: PR is Mercator and is clear.

**One `wgrib2 -d 1 -grid` call per domain would have predicted this before the port.**
It is now in each class's docstring.

### A trap for anyone regridding it

The NA header reads `lat -36.93 to 36.93, lon 299.0 to 60.9`, but **those are rotated
coordinates**. The real geographic extent, from `wgrib2 -ijlat` on the corners, is

```
lon 135.9 - 358.1,  lat -1.6 - 41.5      south pole lat=-35.0 lon=247.0
```

Anyone sizing a target grid from the header values will get a grid in the wrong place.
The file also reports `winds(grid)`, so the winds are grid-relative and any regrid must
rotate them to earth-relative for WPS.

### Option 2 (regrid before ungrib): the command is right, the build is not

Per JH, via Gemini:

```
wgrib2 in.grib2 -set_grib_type same -new_grid_winds earth \
       -new_grid latlon lon_min:num_lons:lon_spacing lat_min:num_lats:lat_spacing out.grib2
```

That is **structurally correct**. `-new_grid_winds earth` is not optional here -- it is
exactly what the grid-relative winds need -- and a regular lat-lon output is template 0,
which ungrib does support. Two caveats on using it as written: the bounds must be
**geographic, not the header's rotated values**, and the full NA domain at 0.1 degrees is
enormous, so subset to a box around the fire (235:101:0.1 35:61:0.1 covers PIUTE at about
the native 13 km).

**It fails here for an unrelated reason.** The installed wgrib2 is v0.2.0.8 (2/2019) and
was built without the interpolation library:

```
$ wgrib2 ... -new_grid latlon 235:101:0.1 35:61:0.1 out.grib2
IPOLATES package is not installed
exit=8, 0 records written
```

So `-new_grid` is unavailable, not wrong. Nothing else on the machine can regrid either
-- no cdo, ncremap, ESMF, NCL or xesmf; `pygrib` 2.1.4 and `pyproj` are present, which is
a possible but laborious route since the hard part is writing valid GRIB2 back out.

**Feasibility: rebuilding wgrib2 with `USE_IPOLATES=1` is the cheap path.** gcc 9.2.1,
gfortran, make and cmake are all present; only the source is missing. Build it to a
local prefix rather than over `/usr/bin/wgrib2`.

### Recommendation

**Do not pursue RRFS_NA for its own sake.** For CONUS fires the 3 km CONUS product is
better anyway; NA was only ever interesting as a 3 km vs 13 km resolution comparison, and
that comparison was already impure because `RRFS_NA` uses `period_hours = 3` while CONUS
uses 1. Per JH, WRF and WPS are expected to be deprecated next year, so WPS gaining
template 1 is unlikely and should not be waited on.

The reason to build wgrib2 with IPOLATES anyway is **general capability** -- it unlocks
any future source on an unsupported grid, not just this one.

## 5n. GRIB retention windows, measured

Per JH: retention governs what can be re-tested and reproduced, so it belongs in the
record rather than in anyone's head. Measured 2026-10-05 by HEAD request against each
host the source classes actually use.

| source | host the class uses | retention | reproducible back to |
|---|---|---|---|
| HRRR | S3 `noaa-hrrr-bdp-pds`, NOMADS fallback | **years** | 730 d verified, archive goes further |
| **HRRR_AK** | **NOMADS only** | **~1 day** | **yesterday** |
| RAP | S3 `noaa-rap-pds` | **years** | 730 d verified |
| **RRFS** (all domains) | S3 `noaa-rrfs-ops-pds` | **56 days** | **2026-08-12** |
| NAM218, NAM227 | S3 `noaa-nam-pds` | **years** | 730 d verified |
| **NAM196, NAM198** | **NOMADS only** | **~7 days** | 7 d |
| GFS | S3 `noaa-gfs-bdp-pds` | **years** | 730 d verified |

Lookback probes, HTTP status by days back:

```
                   2    7   30   60  120  365  730
HRRR_s3           200  200  200  200  200  200  200
HRRR_AK_s3        200  200  200  200  200  200  200
RAP_s3            200  200  200  200  200  200  200
RRFS_s3           200  200  200  404  404  404  404
NAM227_s3         200  200  200  200  200  200  200
GFS_s3            200  200  200  200  200  200  200

NOMADS             1    2    4    7   10   14
HRRR_AK           200  404  404  404  404  403
NAM196 hawaii     200  200  200  200  404  403
NAM198 alaska     200  200  200  200  404  403
```

### Two consequences

**RRFS is the binding constraint, at 56 days.** Every RRFS result has an expiry after
which it cannot be reproduced from the source data. **The PIUTE comparison of §5l rests
on 2026-10-01 cycles and becomes irreproducible around 2026-11-26.** Anything to be re-run
against those GRIBs must be re-run before then, or the files kept deliberately -- they
are in `new_wrfxpy/wrfxpy/ingest/RRFS_CONUS_[PS]` and `ingest/RAP` now, and nothing
protects them from a cleanup. This is also why the two Puerto Rico and Hawaii test fires
JH staged could not be run at their own dates: Bike_Trail is 2026-07-09 and SADDLE_ROAD
2026-08-11, one day before the window opens.

**Three sources are needlessly limited, and it looks like an oversight.** `HRRR_AK`,
`NAM196` and `NAM198` list only NOMADS, where retention is about **1 day** for HRRR
Alaska and **7 days** for the NAM nests. Their S3 equivalents exist and carry the same
paths -- `nam.<date>/nam.t00z.hawaiinest.hiresf00.tm00.grib2`,
`nam.<date>/nam.t00z.alaskanest...`, `hrrr.<date>/alaska/hrrr.t00z.wrfprsf00.ak.grib2`
all return 200 at 30 days back, and HRRR Alaska at 730. In each class the `browse_aws`
line is present but **commented out**, and `git log -S` puts that in `835a893`
(2024-05-26), the commit that first added Hawaii and Alaska ingest. So it was never
enabled rather than switched off after a problem.

**HRRR_AK at one day is effectively non-reproducible**, which matters now that it is the
Alaska replacement for NAM198 (§5m table). Adding `"s3://noaa-hrrr-bdp-pds/"` ahead of
the NOMADS url, and uncommenting `browse_aws`, is a three-line change per class -- left
for JH to decide, since it changes an operational path.

### Practice this suggests

Cache the GRIBs that back any published comparison, rather than relying on the host. The
`ingest/` tree already holds 5.3 TB and nothing distinguishes files that back a result
from files that were merely downloaded once.

### S3 enabled on the three NOMADS-only sources (`45f1a17`)

JH's instinct was that the lines had been commented for a reason, so this was checked
before changing anything.

**Availability was not the reason.** All three products were on S3 on **2024-05-26**, the
exact day `835a893` commented them out, and on 2024-05-20. Lead coverage is complete
today: f00-f60 for both NAM nests, f00-f48 for HRRR Alaska.

**The likely reason is a code coupling.** `grib_source` tests `s3://` availability through
`browse_aws`, and that line is commented out in the same three classes. A `remote_url`
added on its own dereferences a missing attribute, so the pair was probably switched off
together rather than shipped half-wired. Both lines are therefore restored together,
matching what `HRRR` and `NAM218` already do.

**Verified by retrieval, twice.** First with the attributes overridden at runtime for all
three, then with the change as shipped for `NAM196` -- each for a date **30 days back**,
where NOMADS holds nothing, so the files can only have come from S3.

```
HRRR_AK  with s3 -> OK, 3 files for 2026-09-05   hrrr.t00z.wrfprsf01.ak.grib2
NAM196   with s3 -> OK, 3 files for 2026-09-05   nam.t00z.hawaiinest.hiresf01.tm00.grib2
NAM198   with s3 -> OK, 3 files for 2026-09-05   nam.t00z.alaskanest.hiresf01.tm00.grib2
```

**Not verified:** S3-sourced files have not been through ungrib. The paths are identical
to NOMADS so the content should be, but that is an inference, not a measurement. Sizes
differ sharply by product -- 233 MB per lead for HRRR Alaska, 19 MB for hawaiinest,
**775 MB** for alaskanest, which is worth knowing before turning deep-archive access on
for it.

## 5o. Hawaii and Puerto Rico verified through real.exe (2026-10-05)

Using JH's staged fires, re-dated to 2026-10-03 because both sit outside RRFS's 56-day
window (§5n) -- Bike_Trail is 2026-07-09, SADDLE_ROAD 2026-08-11, one day short. These
are **plumbing tests, not forecasts**: the weather no longer matches the fire, so they
prove the WPS chain and nothing about skill.

| | ungrib -> metgrid | levels declared | levels in met_em | real.exe |
|---|---|---|---|---|
| RRFS_PR, Bike Trail | 05:26 -> 05:32 | 46 | **46** | pass, wrfinput 116.4 MB |
| RRFS_HI, Saddle Road | 05:34 -> 05:40 | 46 | **46** | pass, wrfinput 116.4 MB |

Both submitted to the queue and integrated. **So the inherited 46 metgrid levels are
right on the 2.5 km Mercator domains as well as 3 km Lambert, and the HRRR-copied vtable
reads both.** Every domain in the §5m table now has a verified replacement path, and the
Hawaii row -- the only one with a real retirement deadline -- is closed.

**One failure worth recording, because it was mine.** Saddle Road first died with
`Unsatisfiable: RRFS_HI_P forecast is only available for 84 hours`. The job carries a
`cycle_start_utc` pinning the GRIB cycle, and when re-dating the window I shifted
`start_utc`, `end_utc` and the ignition times but not that, so it asked for a 2026-08-11
cycle to serve a 2026-10-04 forecast -- 53 days of lead from an 84-hour product. Bike
Trail has no `cycle_start_utc`, which is why it sailed through. **Re-dating a job means
shifting `cycle_start_utc` too**, and the guard caught it immediately rather than
downloading nothing for an hour.

### RRFS_AK closes the set (2026-10-05)

JH's Kopshesut fire, 67.07 N 157.91 W, above the Arctic Circle, re-dated to 2026-10-03
like the others. `cycle_start_utc` **was** present in this job and was shifted with the
rest -- the §5o lesson applied.

| | ungrib -> metgrid | levels in met_em | real.exe |
|---|---|---|---|
| RRFS_AK, Kopshesut | 06:02 -> 06:18 | **46** | pass, wrfinput 451.3 MB |

**All four runnable RRFS domains are now verified through `real.exe`**, and 46 metgrid
levels held on every one: 3 km Lambert (CONUS), 3 km polar stereographic (AK), 2.5 km
Mercator (HI, PR). Only NA, on the unsupported rotated lat-lon grid, cannot run.

**A check worth keeping: the obvious coverage test would have said no.** The AK grid's
corners read 40.5 N and 64.0 N, which looks like it excludes a fire at 67.07 N. On a
polar stereographic grid the corners do not bound the domain -- the top edge reaches
**75.4 N** at centre, and `wgrib2 -lon` put the fire on a real cell. Checking corner
latitudes on a projected grid is a trap.

Wall clock to `wrf.exe`, which scales with domain size and product count rather than
anything interesting: Puerto Rico 6 min, Hawaii 7, Alaska 16 (61x61 km domain, ~825 MB
per lead), CONUS 41.

## 5p. The forward batch at one pSAF: 31 fires, and a correction to §5i

The 20 forward runs still at pSAF 0.41 were re-run at 0.337, so the forward set is one
batch. **20 of 20 succeeded in 9 min 35 s**, 1.3-4.3 min a fire at 4 concurrent, against
~36 min for a full run -- about **19x cheaper**, because clearing the perimeters while
keeping the wind netcdfs makes it ForeFire-only. Script:
`scratchpad/rerun_psaf337.sh`, unattended, with the 0.41 output preserved in
`psaf410/` inside each fire directory.

### The baseline, now worth quoting

| batch | n | geo mean | median | within 2x |
|---|---|---|---|---|
| 09-30, pSAF 0.41 (§5d) | 20 | 1.87 | 2.01 | 45% |
| 10-01, pSAF 0.337 (§5i) | 11 | 1.69 | 1.82 | 36% |
| **one batch, pSAF 0.337, full window** | **31** | **1.42** | **1.49** | **58%** |
| **one batch, pSAF 0.337, 8 hours** | **31** | **1.40** | **1.36** | **58%** |

Reverting 0.41 to 0.337 did what §5d intended: 1.87 to 1.42, within-2x from 45% to 58%.

### Correction to §5i

§5i reported that §5d's 8-hour prediction for 0.337 -- geo mean 1.42, median 1.30, **78%
within 2x** -- "does not hold out of sample", on the strength of 11 fires giving 1.63 /
1.87 / 36%. On 31 fires at the same setting the measurement is **1.40 / 1.36 / 58%**.

**§5d predicted the centre almost exactly and the spread badly.** That is a different and
milder conclusion than §5i's. The 11-fire batch was too small and happened to hold
several extremes; the lesson is the one §5h already recorded about n=14, applied again
and ignored once.

### Three fires worth watching

- **Finley_RX collapsed to 0.06** (12.5 ha against 203.8). At 0.41 it was moisture-limited
  and barely moved; at 0.337 it dies almost entirely. §5g's extinction threshold, now on
  the wrong side -- and the clearest case yet that pSAF and moisture are not independent
  knobs on a moisture-limited fire.
- **Pulaski 6.57**, still the worst, as §5f expected: industrial Little Rock, fuel maps
  unlike the ground. Excluded from calibration by JH.
- **Lauderdale tripped the shrink guard** at 4.2% on step 26 -- above the 1% warning,
  below the 10% skip, so it scored and was flagged. First time §5j's guard has fired on
  genuinely new output rather than on a regression test.

### Timing, measured, for sizing an ensemble

Split by file mtimes over runs whose mtimes are intact (the retagged ones are useless for
this -- their geojson times are the retag):

| | total | WindNinja + netcdf | ForeFire | per step |
|---|---|---|---|---|
| alone (Bigfoot) | **19.4 min** | 18.5 (96%) | 0.9 (4%) | 0.69 min |
| 4 concurrent (n=20) | 35.7 min | 32.7 (92%) | 2.4 (7%) | 1.27 min |

28 steps, 27 h, 50 m, 31 km domain, GRIBs and FMDA already cached. **A full-duration
forecast is about 20 minutes and WindNinja is essentially all of it**; an 8-hour forecast
is about 6, since cost scales linearly with steps.

### What that means for a perturbed ensemble

Moisture, ignition location and pSAF **leave the wind field untouched**, so the netcdfs
are built once and shared. Measured resource use per ForeFire member: **32 MB RSS and
3-28% of one core** on a machine with 240 GB free and 32 cores, so 50 members is ~1.6 GB
and nowhere near saturating the ~25 cores JH expects to have.

| | 27 h | 8 h |
|---|---|---|
| netcdf + WindNinja, once | 18.5 min | ~5.5 min |
| 50 members, ForeFire only, serial | ~25 min | ~7 min |
| 50 members at 25-way | ~1 min | ~0.3 min |
| **total at 25-way** | **~20 min** | **~6 min** |

**A 50-member ensemble costs about what one deterministic forecast costs**, provided the
expensive part is shared. Three conditions, none of them yet met:

1. **The domain must be fixed, not recentred per member.** `forefire_grib` sizes the
   domain around the ignition; if each member recentres on its perturbed ignition, every
   member needs its own WindNinja pass and the cost becomes ~15 hours. One design choice,
   40x.
2. **Moisture is not a `setParameter`.** `sweep_job` perturbs anything that is; moisture
   reaches ForeFire through the per-step `fuelsTableFile`, so it needs `md_series` offset
   before `write_fuel_table` (§5h). Ignition location is likewise in the `.ff` script.
3. **`sweep_job` is serial.** Members are independent, so this is where the remaining
   gain is -- and the 25-way figure above is an **estimate, not a measurement**. Container
   launch and filesystem contention over 28 steps x 50 members is 1,400 apptainer
   invocations, and nothing has tested whether that stays linear. One fire, 50 variants,
   timed, would settle it.

## 5q-5r. Moved to the 2026-10-06 handoff

The 10-06 work -- the perturbed ensemble and the fuel-table bug behind it, the container
overhead, and the finding that every fire validated so far is a known incident -- lives in
`SESSION_HANDOFF_forefire_2026-10-06.md`, keeping its §5q/§5r numbering so references from
the commit messages still resolve.

## 6. Open items

Carried from `SESSION_HANDOFF_forefire_2026-10-01.md` §6, unchanged, plus:

1. ~~**Validate AK, HI and PR through metgrid**~~ **DONE (§5o).** All four runnable RRFS
   domains pass `real.exe` with 46 metgrid levels confirmed on each.
2. **Build wgrib2 with IPOLATES** only if the regrid path is ever wanted, to a local
   prefix. Per JH, RRFS_NA is being skipped, so this is optional and general-capability
   rather than needed. The one-file test in §5m is the proof if it is built.
3. **Cache the GRIBs behind published comparisons** (§5n). RRFS keeps 56 days; the PIUTE
   result ages out around **2026-11-26** and nothing in `ingest/` marks which files back
   a result.
4. **S3-sourced NAM nest files have not been through ungrib** (§5n). The paths match
   NOMADS so the content should, but that is an inference.

## 7. NEXT SESSION

Three, in the order they unblock each other.

**1. Time a real ensemble: one fire, 50 variants, measured.** §5p puts a 50-member
ensemble at about the cost of one deterministic forecast -- but the 25-way figure is an
**estimate**. 28 steps x 50 members is 1,400 apptainer invocations and nothing has tested
whether container launch and filesystem contention stay linear. The work needed to
measure it is the work needed to build it: member-level parallelism in `sweep_job`, which
is serial today. Do this first, because the answer sizes everything else.

Two conditions come with it, from §5p:

- **Fix the domain; do not recentre per member.** `forefire_grib` sizes the domain around
  the ignition. If each perturbed ignition recentres it, every member needs its own
  WindNinja pass and ~20 minutes becomes ~15 hours. One flag, 40x.
- **Moisture is not a `setParameter`.** It reaches ForeFire through the per-step
  `fuelsTableFile`, so a moisture ensemble needs `md_series` offset before
  `write_fuel_table` (§5h). Ignition location is likewise in the `.ff` script.

**2. Finley_RX at 0.06 is a signal, not an outlier.** At pSAF 0.41 it was moisture-limited
and barely moved; at 0.337 it dies almost entirely, 12.5 ha against WRF-SFIRE's 203.8.
The batch now contains a fire the global pSAF kills outright, which says plainly that
**pSAF and moisture are not independent knobs** on a moisture-limited fire. §5g has the
mechanism -- fuel 2 above extinction for 18 h -- and §5f showed pSAF has no leverage
there. This is the strongest argument yet for the moisture-perturbed ensemble of §5h,
since it is the one axis that moves such a fire at all.

**3. Isolate the mesh.** Open since the 09-29 handoff §7 and now the **only** uncontrolled
variable between the two ForeFire modes: the forward set runs at 50 m, the cron set at
25 m. Until that is separated, the forward-vs-cron ratio mixes forecast-vs-analysis with
mesh. Re-run the forward set at 25 m, or the cron set at 50 m, and the comparison in the
last session's three-way table becomes attributable.

### The forward baseline to quote

```
forward mode, pSAF 0.337, 31 fires, 50 m + barrier mask
  full window   geo mean 1.42   median 1.49   58% within 2x   19% below 1x
  8 hours       geo mean 1.40   median 1.36   58% within 2x   26% below 1x
```

34 forward runs exist on disk; 31 score. Bigfoot has no WRF-SFIRE workspace, being the
VIIRS-only fire the forward path exists for, and the two Walker Hill pile burns
extinguished (§5i).
