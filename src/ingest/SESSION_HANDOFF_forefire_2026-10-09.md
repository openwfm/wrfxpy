# ForeFire session handoff — 2026-10-09

HAWK, a fire that crossed I-5 at the Grapevine, run as a 2x2 of fuel map against ignition
point. The remap verdict from 10-07 survives and strengthens, but the session's real
finding is that hard barriers make a forecast hypersensitive to an ignition error that the
detection cannot resolve. Separately, GFS winds now reach WindNinja.

**Section numbers continue from the 2026-10-07 handoff**, which continues 10-06, 10-05 and
10-01. Read this as the 10-09 continuation.

## 5u. HAWK: the 2x2 that separates fuel map from ignition side

### Why a 2x2 and not a pair

The pair was set up as PICO and Pulaski were: original raster against `noburn91`, winds and
moisture identical. JH then observed that **both NGFS jobs for this fire put the ignition
on the wrong side of the interstate**, which makes a pair unreadable -- a smaller no-burn
fire could be the barrier holding, or the fire simply confined to the wrong side. Two more
members were queued at a corrected ignition before the first pair finished.

### The interstate, and a correction

A first attempt located I-5 by searching the fuel raster for a dense north-south corridor
of code 91 and found one at lon -118.848. **That was wrong** -- near the ignition latitude
I-5 runs NW-SE, and the corridor at -118.848 is a different road 4.5 km away. Isolating
corridors by morphological opening (interstates are 3+ cells wide where local roads are
one) found the real one: a 10 km chain at bearing 322 deg passing ~0.6 km east of the
ignition. JH then supplied an actual centreline, which is what the scoring now uses
(`i5_line.py` in the session scratchpad, 10 vertices, lat 34.689 to 34.792).

Both methods agree on every side assignment and on distances to within ~90 m, so they
corroborate each other; the supplied line removed the need to extrapolate past the ends of
the detected segment.

### Everything points east except the forecasts

| point | side of I-5 | dist to centreline |
|---|---|---|
| VIIRS 2-instrument centre (JH) | **EAST** | 392 m |
| `forecast_ign` HAWK | WEST | 613 m |
| `forecast_ign` DANNY | WEST | 125 m |
| DANNY InitialFcast | WEST | 17 m |
| terrain-corrected pixel centroid | EAST | 1568 m |
| SWIR centroid | EAST | 75 m |

`ign_latlon` is the mean of SWIR pixel centroids over the first three hours
(`src/ngfs/ngfs_incident.py` ~line 779), used whenever it agrees with the nominal mean to
within 0.04 deg -- 4.4 km, far too loose to catch this. But no choice among these points is
defensible: **`satzen` is 61.43 deg and `pixel_area` 12.50 km2**, a GOES-19 pixel stretched
to roughly 5.5 x 3.2 km, with `viirs_pixel_ign` False. The side call is being made an order
of magnitude inside the sensor's own footprint.

JH's best estimate, the centre of detections from two VIIRS instruments, is
**34.735403, -118.793339**. VIIRS resolved it but too late to act on -- which is the
low-latency argument in its sharpest form: the estimate that arrives in time cannot place
the fire, and the one that can place it arrives too late.

### HAWK and DANNY are one detection

`HAWK {9F712D58-...}` and `DANNY {22AECF2E-...}` carry the **identical
`feature_tracking_id`** (`I6-2026-10-08T23:11:17.000Z_0069`), the identical `uuid`, and both
read `local_id_and_known_incident_name = "362527 - HAWK"`. Only the incident ID and name
differ. Two full forecasts were run for one hotspot. **De-duplicate on `feature_tracking_id`
or `uuid`, never on `incident_name` or `incident_id_string`.**

### WRF-SFIRE's ignition seed decides two of the three jobs

`fire_ignition_radius1 = 60.0` m with `fire_ignition_ros1 = 1.0` m/s over a 1200 s window.
Measured from TIGN_G, the seed reaches 76 m equivalent radius in five minutes, 98 m by the
end of the ignition window and 124 m by 00:00 -- independent of fuel. Against the distances
above:

- **DANNY InitialFcast, 17 m** -- inside the 60 m seed at t=0. On both sides of I-5 the
  instant it ignites, under any fuel map. The crossing carries no information.
- **DANNY operational, 125 m** -- crosses about 48 minutes after ignition from seed growth
  alone.
- **HAWK, 613 m** and the **VIIRS centre, 392 m** -- outside the seed, so their crossings
  are real spread.

ForeFire's template seeds a bare point (`startFire[lonlat=...]`, no radius, no forced ROS).
**Any ForeFire-vs-WRF-SFIRE comparison on a fire ignited within ~150 m of a barrier is
measuring the ignition treatment, not the models.** That is a systematic bias toward
WRF-SFIRE crossing barriers, and it will show up in the batch wherever small fires sit near
roads.

### The crossing, and the 09:00 discrepancy

JH read a crossing near 09:00 UTC off a PNG. Measured against the centreline, WRF-SFIRE
with the original map **arrives** at the highway at 09:00 and stalls at 115-124 m for five
hours, then breaks through at 16:00-17:00, finishing 1774 m past it. Arrival at 09:00,
breakthrough later; the PNG explains the difference. The five-hour stall is the barrier
working *partially* even in the original map, which makes the no-burn comparison a
hardening of a partial barrier rather than the creation of one.

### The 2x2, against an observed outcome

JH saw a perimeter for this fire (provenance unknown, file not available): the real fire was
**about 150 acres = 60.7 ha, heavily suppressed**, constrained by the highway, with a pocket
cut out around a building cluster at 34.735439, -118.796673.

| member | final ha | west ha | east ha | x observed | abs err |
|---|---|---|---|---|---|
| **FF GOES noburn91** | 13.7 | 12.2 | 0.0 | **0.23x** | **47 ha** |
| WRF-SFIRE noburn91 | 259.7 | 259.7 | 0.0 | 4.3x | 199 ha |
| WRF-SFIRE orig | 1822.8 | 1233.4 | 589.4 | 30.0x | 1762 ha |
| FF GOES orig | 4965.9 | 1402.1 | 3513.4 | 81.8x | 4905 ha |
| FF VIIRS noburn91 | 4795.5 | 18.9 | 4729.0 | 79.0x | 4735 ha |
| FF VIIRS orig | 9638.0 | 1559.2 | 7998.8 | 158.8x | 9577 ha |

**Both models agree the remap stops the crossing**: WRF-SFIRE 589 -> 0 ha east, ForeFire
3513 -> 0 at the GOES seed, 1559 -> 19 ha on the far side at the VIIRS seed. Reproducible
across two independent models.

**Every original-map member is 30-159x over; both no-burn members land within ~4.4x.**

### Three things this got wrong at first, recorded because they recur

1. **The 13.7 ha member was called a collapse artefact.** It is the best of the six by
   absolute error. This is the Pulaski lesson -- against a small suppressed fire, near-zero
   growth is a good estimate -- and it was missed again.
2. **But it is right for the wrong reason.** The GOES seed is west, so that member is
   12.2 ha west and 0 east while the real fire was east: correct area, essentially zero
   overlap. **Area alone rewarded the member that got the geography backwards.** Score
   shape as well as area on fires like this.
3. **The ignition cell was not the cause.** It is FBFM5 in both rasters. The 13.7 ha is the
   fire sealed into a small road-bounded parcel.

### The buildings pocket is reproduced, by the right mechanism

The cluster at 34.735439, -118.796673 sits 89 m east of the centreline. The coordinate
itself lands on an FBFM5 pixel, but the 7x7 neighbourhood is dense code 91. In `noburn91`
that block becomes 14, which is in `BARRIER_CODES`; in the original map `remap_fuel` fills
91 from the nearest burnable neighbour. The perimeters follow:

- FF GOES orig -- buildings **inside** the burn, 206 m in from the edge
- FF VIIRS noburn91 -- buildings **spared**, fire edge stopping **59 m** short

That is shape skill, not just area. JH's suppression explanation and the fuel-map
explanation are not exclusive: the map puts a break there and crews would work a break that
already exists.

### The finding that matters operationally

Moving the ignition 1124 m changes the answer by **1.94x** under the original map
(4966 -> 9638 ha) and by **350x** under the remapped map (13.7 -> 4795 ha). Hard barriers
put the seed in a road-bounded parcel, and which parcel it lands in decides the forecast.
The 600 m error that selects the parcel is far below the pixel that produced it.

**Per JH: the ignition point should always be an ensemble axis -- it accounts for the
largest part of forecast variability.** This fire is the strongest case yet, because the
barrier map converts a sub-pixel placement error into a 350x spread.

The remap is still preferred -- on this fire it is the difference between a ~4x error and a
30-159x error, and it reproduced an observed structure-sparing pocket -- but never from a
single ignition point.

## 5v. GFS winds through WindNinja

### The silent failure

WindNinja has **no GRIB2 reader for GFS**. `ncepGfsSurfInitialization` is registered as
UCAR-GFS-GLOBAL-0.5-DEG and reads netCDF; only HRRR (`NCEP-HRRR-3km-SURFACE`) and NAM 12 km
have direct GRIB2 readers. Handed a cached `gfs.tHHz.sfluxgrbfNNN.grib2` (Gaussian grid,
template 40) WindNinja parses the valid time correctly, then:

```
ERROR 1: Point outside of projection domain
Run 0: Solver time was 0.000000 seconds.
=== exit: 0 ===
```

It writes an **all-zero wind field** and a saturated 100% cloud field, exit code 0. The data
is fine -- the GRIB has U=2.09, V=2.59 m/s at the fire.

`_check_wind` in `src/ingest/forefire_grib.py` only rejected speeds above 100 m/s (the
NAM218 garbage signature), so this passed. **A GFS-driven ForeFire run would have produced
a no-wind forecast that looked successful.** The guard now also rejects a maximum of zero
outright -- a mass-conserving solver cannot return an identically flat field over terrain.
HRRR fields verified still passing. Docstring updated to describe both modes.

### The conversion, solved

`src/ingest/gfs_windninja.py` (new) turns cached sflux GRIB2 into a file WindNinja reads.
Four things all had to be right, each found by bisecting against a file WindNinja downloaded
itself:

1. **The path selects the reader.** The file must sit at
   `<root>/UCAR-GFS-GLOBAL-0.5-DEG-<dem_basename>/<stamp>/<stamp>.nc`, stamp `YYYYMMDDTHHMM`.
   The same bytes elsewhere are matched to a different reader and fail with
   `Failed to find height for height_above_ground4`.
2. **Structure** -- lowercase `u-`/`v-` names, `latitude`/`longitude` (latitude descending,
   longitude 0-360), convective cloud on its own time dimension, and a `latitude_longitude`
   grid-mapping variable referenced by every data variable.
3. **The `Grib2_*` variable attributes** and the `Analysis_or_forecast...` / `history` /
   `History` globals. Structurally identical without them, it is still rejected.
4. **Full level sets** -- `height_above_ground2` = 10/20/30/40/50/80/100 and
   `height_above_ground3` = 2/80/100, or it crashes on `GetRasterBand(7) - Illegal band #`.

**One valid time per file.** A multi-time file reads back inconsistently: writing new data
into every time slice changed only the first WindNinja run, runs 2 and 3 reproducing the
untouched file byte for byte. Per-time files give four distinct outputs tracking the source:

| valid time | GFS domain mean | WindNinja 6.1 m |
|---|---|---|
| 15:00Z | 2.64 | 2.71 |
| 18:00Z | 2.97 | 2.88 |
| 21:00Z | 2.83 | 2.42 |
| 00:00Z | 2.11 | 1.58 |

This matches how `forefire_grib.windninja_uv` already drives WindNinja -- one GRIB, one
valid time, one call -- so wiring it in should be small.

WindNinja's own downloader (`--wx_model_type "UCAR-GFS-GLOBAL-0.5-DEG" --forecast_duration N`)
also works and was verified, but the UCAR THREDDS `best` dataset is a rolling aggregate, so
a fetch is **not reproducible after the fact**. The cached GRIB2 is the reproducible
artifact, which is the argument for the converter.

### Why this is worth having, per JH

Not resolution -- GFS is 0.25-0.5 deg against HRRR's 3 km. **Lead time.** `GFSF` is already
configured for `grib_forecast_hours_periods = [{'hours': 384, 'period': 3}]`, hourly to 120 h
then 3-hourly to 384. At day 5-16 the comparison is not GFS against HRRR, it is GFS against
nothing. JH's target use is **ensemble forecasts for prescribed fires**, where location and
ignition time are both known, so the two error sources that dominated HAWK collapse and
weather plus fuel moisture carry the spread.

### JH's alternative, probably the stronger path

Run WRF-SFIRE off GFS forcing and take winds *and* fuel moisture from the wrfout, starting
from equilibrium moisture estimated from GFS variables. That is the best moisture estimate
available at long range, FMDA being CONUS-only by choice, and it gives physically consistent
wind and moisture from one solution rather than two fields that never saw each other.
WindNinja also ships a `wrfSurfInitialization` reader, so it can initialise from wrfout
directly if both are wanted.

The tradeoff is cost per member, and a prescribed-fire ensemble is where member count
matters most -- which suggests the two are complementary: WRF-SFIRE off GFS for a few
scenario members, direct GFS -> WindNinja for the large cheap ensemble around them.
**Worth deciding before building either out, since it changes what the converter is for.**

## 6. Open items

Carried forward from 10-07 and still open: blanket code 91 may be too aggressive (NLCD 23/24
is the ready refinement); rural roads invisible to both FBFM13 and NLCD; the fragility
measure needs re-checking; agriculture code 93 still filled burnable.

New from this session:

1. **Wire ignition location in as a standing ensemble axis**, sampled over the detection
   footprint rather than a point. This is now the top item -- per JH it is the largest
   source of forecast variability, and 5u quantifies a 350x spread under the barrier map.
   `satzen` and `pixel_area` give the sampling width directly.
2. **De-duplicate NGFS jobs on `feature_tracking_id` / `uuid`** before launching or scoring
   a batch. HAWK/DANNY doubled the compute for one hotspot and would double-count it in any
   batch score.
3. **Decide GFS -> WindNinja versus GFS -> WRF-SFIRE -> wrfout** before building further.
   See 5v.
4. **Score shape, not only area**, on fires near barriers. On HAWK the best member by area
   had essentially zero overlap with the real fire.
5. **The ForeFire point seed versus WRF-SFIRE's 60 m seed** is a systematic bias wherever a
   small fire sits near a barrier. Quantify it across the batch before trusting those
   comparisons.
6. **ForeFire ran 2.7x WRF-SFIRE** on the matched HAWK configuration and 5.3x at the VIIRS
   seed, well above the 1.42 batch geometric mean. HAWK is not typical of the pSAF
   calibration set; worth knowing why before adding it to any fit.
7. **`_md_offset`, mesh isolation (forward 50 m vs cron 25 m), Finley_RX at 0.06, and
   RRFS_AK through metgrid** remain untouched from earlier handoffs.

## 7. State at end of session

**Committed 2026-10-12, on `james_ngfs` in `/data/jhaley/wrfxpy`:**

- `61417b1` `src/ingest/forefire_grib.py` -- `_check_wind` all-zero guard plus docstring
- `689894f` `src/ngfs/ngfs_incident.py` -- JH's `base_folder` change at line 90, pointing
  behave FMDA at `/data/jhaley/clean_wrfxpy/wrfxpy/wksp_fmda/CONUS/` (the old path stops
  at 202609). Still hardcoded; meant to move into a configuration file.
- this handoff, in the commit following `689894f`

**Deliberately left uncommitted:**

- `src/ingest/gfs_windninja.py` -- the GFS converter, verified end to end on the HAWK case
  but held until it has been exercised on long-range forecasts. Known nits before it is
  committed: the usage line names `gfs_to_windninja.py`, the docstring lists two
  requirements where 5v lists four, and `EPOCH` is hardcoded to 2026-09-09.

**ForeFire outputs**, 50 geojsons each, 25 m mesh, issue 2026-10-09T00:00:00Z, 24 steps:
`wksp_forefire/HAWK_{fuelcmp25,viirs25}_{orig,noburn91}`.

**Scoring scripts**, copied out of the session scratchpad so they survive it, are in
`wksp_forefire/`: `i5_line.py` (JH's I-5 centreline and the side/distance test) and
`hawk_compare.py` (the 2x2 scorer). Move `i5_line.py` into `src/ingest/` if barrier scoring
continues; it is currently sourced by path. The scratchpad also held the converter's
development history and the reference THREDDS netCDF WindNinja downloaded -- that reference
file is gone with the session, but `gfs_windninja.py` no longer needs it.

**Two process notes.** A watcher reported "pipeline died" falsely because its `pgrep`
pattern named a script whose contents had been renamed but whose filename had not. And
`pkill -f <pattern>` self-matched and killed the running shell again -- the hazard is
recorded and was still repeated; kill by PID.
