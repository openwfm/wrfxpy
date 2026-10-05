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

## 6. Open items

Carried from `SESSION_HANDOFF_forefire_2026-10-01.md` §6, unchanged, plus:

1. **Build wgrib2 with IPOLATES** if the regrid path is wanted, to a local prefix.
   Then the one-file test above is the proof: regrid, `wgrib2 -grid` to confirm
   template 0, and run ungrib on the result.
2. ~~**Validate AK, HI and PR through metgrid**~~ **DONE.** All four runnable RRFS
   domains pass `real.exe` with 46 metgrid levels confirmed on each.
