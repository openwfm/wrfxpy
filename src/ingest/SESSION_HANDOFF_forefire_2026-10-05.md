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

## 6. Open items

Carried from `SESSION_HANDOFF_forefire_2026-10-01.md` §6, unchanged, plus:

1. **Build wgrib2 with IPOLATES** if the regrid path is wanted, to a local prefix.
   Then the one-file test above is the proof: regrid, `wgrib2 -grid` to confirm
   template 0, and run ungrib on the result.
2. **Validate AK, HI and PR through metgrid.** Their grids are supported, but their
   `num_metgrid_levels` are still the module's guesses, as NA's was. §5l confirmed 46
   for CONUS only.
