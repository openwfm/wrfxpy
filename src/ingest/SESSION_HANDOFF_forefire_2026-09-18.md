# ForeFire session handoff — 2026-09-18

Branch `james_ngfs` at `4ef50da`. **No commits made this session.** Two source
files were edited on explicit instruction and are left uncommitted with backups
(§8). Follows `SESSION_HANDOFF_forefire_2026-09-11.md`, whose §9 item 14 asked
where the ForeFire gap lives. This session answers part of that, and it is not
the wind.

**Headline:** ForeFire was run twice on the same fire with the only difference
being the wind field — coupled WRF `UF`/`VF` against WindNinja/HRRR. Scored
against a real IR perimeter, the two arms produce **identical elongation, 1.73
against an observed 2.90**. Swapping the entire wind field moved area 15% and
shape not at all. The shape deficit is structural.

---

## 1. WindNinja is installed and verified

Conda env `windninja`, conda-forge `windninja-3.13.0.1` (WindNinja 3.13.0, SCM
382cd19, released 2026/03/09, OpenMP 32 threads). Separate env on purpose;
`wrf_test` and `wrfx` untouched. This closes the 09-08 §9 "investigated, nothing
installed" item — the prebuilt CLI removed the compiler problem exactly as that
section predicted.

**It ingests the raw HRRR `wrfprs` GRIBs directly**, which 09-08 §9 listed as
genuinely unresolved:

    --initialization_method wxModelInitialization
    --forecast_filename ingest/HRRRA/hrrr.20260918/conus/hrrr.t00z.wrfprsf03.grib2

It auto-identifies them as `NCEP-HRRR-3km-SURFACE` and resolves the valid time
correctly (t00z + f03 -> 03Z). No subsetting, no fetcher, no reformatting.
`NOMADS-HRRR-CONUS-3-KM` and `PASTCAST-GCP-HRRR-CONUS-3-KM` are both in its model
list.

### What the cached FMDA GRIBs contain

`ingest/HRRRA/hrrr.<YYYYMMDD>/conus/hrrr.tNNz.wrfprsf03.grib2`, ~415 MB each,
711 messages, hourly cycles, **f03 only**, 5 days retained, ~9.6 GB/day, nothing
prunes it.

- `UGRD`/`VGRD` on 40 pressure levels (50-1013.2 mb) plus **10 m** and **80 m** AGL
- `GUST`, `WIND` (10 m max), storm motion, shear, helicity
- everything WindNinja's surface init wants: `TMP`/`RH`/`DPT`/`SPFH` at 2 m,
  `PRES:surface`, `TCDC`/`LCDC`/`MCDC`/`HCDC`, `DSWRF`, `HPBL`, `SNOWC`
- plus `SFCR` (roughness), `FRICV` (friction velocity), `HGT:surface` (msg 608,
  the same message `static/hrrr.terrainh.nc` was built from)

Grid: Lambert 1799 x 1059, 3 km, LoV 262.5, Latin1 = Latin2 = 38.5.

### DEM

`landfire/LF2020_Elev_220_CONUS` (30 m, EPSG:5070) was already staged and is the
right source. Warped to the local UTM zone it matches WRF's own `ZSF` to
**RMS 2.2 m, correlation 1.0000** over 1.44 M points — so DEM choice is not a
confound in any of the comparisons below.

---

## 2. Three traps, all of the silent kind

1. **`--output_speed_units` defaults to `mph`.** Feeding that into ForeFire's
   `windU`/`windV` inflates winds by 2.24x with nothing looking wrong. Always
   pass `mps`.

2. **WindNinja does NOT read raw NAM218 `awphys` files.** It exits 0, prints a
   plausible valid time, and produces **speeds of 25.8 to 45,632 m/s**, mean
   9985. Checking the exit code is not enough; check the values. This matters
   because every operational wrfxpy run is NAM218-driven, so a same-driver
   controlled comparison cannot be built by pointing WindNinja at the job's own
   GRIB. It needs an HRRR-driven WRF run (§6).

3. **HRRR winds are grid-relative** (`uvRelativeToGrid = 1`), rotation
   `sin(38.5 deg) * (lon - 262.5)`, up to ~16 deg at CONUS edges. **WindNinja
   applies this internally** — verified on a coherent 6.7 m/s Nevada case:
   WindNinja 307.2 deg against earth-relative 306.9 (0.3 deg) and grid-relative
   320.4 (13.2 deg). It reads the **10 m** field, not 80 m. So the rotation
   caveat applies only to our own extraction code, never to the WindNinja path.
   `src/ingest/HRRR.py`'s own consumers are unaffected; the helper written this
   session lives in the session scratchpad only.

**Method note worth not re-learning.** The first attempt to verify the rotation
used a light-wind evening case (TELEPHONE) and appeared to confirm it —
WindNinja 248.9 deg against earth-relative 248.3. That was luck: the field had
circular resultant length 0.149, directions scattered across the whole compass.
**A direction comparison on a field with resultant length below ~0.3 cannot
distinguish anything**, and a single-cell check there will look like agreement.
Use a spatially coherent field.

---

## 3. Cost and the cell ceiling

Cost depends **only on total horizontal cell count**, not on domain size or mesh
resolution separately. Matched pairs agree to 0.2 GB and ~7%:

    domain   mesh   horiz cells   wall    peak RSS
    120 km    75 m     2.56 M     268 s    54.4 GB
     80 km    50 m     2.56 M     294 s    54.2 GB
    120 km    60 m     4.00 M     423 s    84.7 GB
    100 km    50 m     4.00 M     407 s    84.6 GB

Both `wall` and `RSS` scale as `cells^0.93`. **The budget is between 4.0 M
(runs) and 5.76 M (fails)** cells, to spend however you like between extent and
resolution. True 50 m works today on domains up to ~100 km; the 120 km fire
subgrid at 50 m is the one combination just over the line.

**The 5.76 M failure is an int index overflow, not memory** —
`std::bad_array_new_length` at 6.8 GB peak RSS with 232 GB free, which WindNinja
reports as "appears to have run out of memory". More RAM will not fix it.
`--ascii_out_resolution` clears it by solving on a mesh it can index and writing
output at the requested resolution.

---

## 4. Convergence is terrain-dependent, and that is the real constraint

`--ascii_out_resolution 50` produces a correctly-shaped 2400x2400 file from any
solve mesh, but **interpolated output is not resolved detail**.

Complex terrain (Wind River, ZSF std 472 m), RMS difference against a 75 m solve:

    300 m  0.526 m/s     150 m  0.297 m/s
    200 m  0.396 m/s     100 m  0.217 m/s

Roughly halving with resolution and **no plateau**. The field keeps gaining
structure rather than settling: std 1.302 -> 1.446 m/s and peak speed
**11.3 -> 20.0 m/s**. The ridge-top accelerations that drive fast spread are the
least converged part.

Flat terrain (Dry River, ZSF std 26.7 m), 50 m against 25 m: **0.086 m/s RMS
(1.3%) and 0.5 deg**. Converged.

**Set the solve mesh from terrain roughness, not a fixed rule.** And the honest
description of WindNinja in mountains is terrain-aware interpolation, not
prediction — it is mass-consistent, not momentum-solving, so it has no
representation of separation, lee rotors or hydraulic jumps, which is what
dominates mountain winds. JH's standing view, and the measurements agree.

---

## 5. `UF`/`VF` already carry `windrf`

Measured on the fire grid: `UF`/`VF` are **not** simply the wind at
`fire_wind_height` (6.096 m). They carry the per-fuel-category `windrf`
reduction from `namelist.fire` (0.36 for most categories here):

    WRF |UF,VF|           mean 0.939 m/s
    WRF |UF,VF| / windrf  mean 2.483 m/s
    WindNinja @ 6.096 m   mean 3.090 m/s

**Any external wind field written into the slot `UF`/`VF` fills must have
`windrf` applied per fuel cell, or the fire is driven ~3.3x too fast.** ForeFire
then applies its own `windReductionFactor` on top of whatever it is given. This
is the quantified version of the 09-10 §18 "double reduction". Compare like with
like by dividing `|UF,VF|` by `windrf` before checking against a diagnostic model.

Category 14 (no fuel) has `windrf = 1e-7`; mask those cells rather than dividing.

---

## 6. WindNinja against WRF-SFIRE — the controlled comparison

Dry River, 2026-09-15 19:00 UTC, flat plains, 1200x1200 fire grid at 25 m.
**Two WRF-SFIRE runs of the same fire on an exactly identical grid and terrain**
(max |dZSF| = 0.000 m), one NAM218-driven and one HRRR-driven, against WindNinja
from HRRR. 1,332,871 cells, burnable, fire-free and above 0.5 m/s in both runs,
`windrf` divided out.

    comparison                     B/A   dir bias  dir RMS  anomaly corr
    WRF(HRRR) vs WindNinja(HRRR)  1.04      +2.4      3.7      0.129
    WRF(NAM)  vs WindNinja(HRRR)  1.13      +2.7      4.2      0.126
    WRF(NAM)  vs WRF(HRRR)        1.09      +0.3      3.0      0.874

Controlling for the driver cuts WindNinja's speed excess from 13% to **4%**, so
most of the gap was HRRR-vs-NAM, not the downscaling. Its direction error
(3.7 deg RMS) is barely above the **3.0 deg RMS between two WRF runs**, i.e.
inside the inter-model spread.

**But the spatial anomaly structure is uncorrelated** — 0.13 where two WRF runs
reproduce each other at 0.87 — while the variance amplitudes are comparable
(CV 3.1-4.3%). WindNinja makes variability of the right size in the wrong
places. Its texture comes from mass conservation over terrain; WRF's comes from
boundary-layer structures.

A light-wind evening case (WOODRIDGE) gave +64 deg bias and 83 deg RMS. That was
the regime, not an incompatibility — see the method note in §2.

---

## 7. The two-arm ForeFire test — the result that matters

Dry River. Both arms identical except `windU`/`windV`: arm A from the wrfout's
`UF`/`VF`, arm B from hourly WindNinja/HRRR with `windrf` applied. Same wrfouts,
fuels, terrain, ignition, fuel table and ForeFire parameters; the netcdfs are
built by the same code path and arm B overwrites two variables afterwards. All
50 arm-B netcdfs verified tagged and carrying WindNinja winds.

Scored against the IR perimeter
`ngfs/perims/Dry_River_{F5FAD4B1-44AA-4FB5-940A-B16C30619A2B}.geojson`,
`poly_PolygonDateTime` 2026-09-16 07:02:58 UTC, 2078.5 ha (geometry verified
against its own 5136.5-acre property), nearest ForeFire step -3 min:

    source                   area ha  obs x    IoU   elongation
    IR observation            2078.5   1.00  1.000       2.90
    WRF-SFIRE (HRRR)           619.2   0.30  0.235       3.13
    ForeFire <- WRF UF/VF     1537.2   0.74  0.330       1.73
    ForeFire <- WindNinja     1763.5   0.85  0.341       1.73

**The wind source does not fix the shape.** Both arms come out at elongation
**1.73** — identical to two decimals — against an observed 2.90. This repeats
the 09-11 SINLAHEKIN result (modelled 1.38-1.57 vs observed 2.54) on a different
fire in the plains regime. Further wind work cannot close this gap.

**The two models fail in opposite directions.** WRF-SFIRE gets the shape nearly
right (3.13 vs 2.90) but is a third of the observed area. ForeFire gets closer on
area but is far too round. ForeFire's better IoU comes from having more area
roughly in the right place, not from better shape.

**WindNinja costs nothing in skill** — marginally better on both area and IoU.
For the low-latency path that is the answer: dropping the dependency on a
completed WRF-SFIRE run did not degrade the forecast on this fire.

Across all 50 steps the arms diverge gradually: IoU 1.000 at step 2 falling to
0.820 at step 50, with B ending 22% larger.

### Use elongation, not reach, for shape

**Per JH, ignition points in observation files are always suspect.** Reach from
the ignition point inherits that error; elongation (major/minor axis of the
polygon itself) does not. For the record the estimated ignition sits 71 m from
the observed perimeter edge with the centroid 4.3 km away, consistent with a
wind-driven plains fire, but that is not independent confirmation of either.

This also qualifies 09-11 §9 item 15: `reach ~ pSAF^1` transferred across three
fires and is still the physically robust relation, but it is only *measurable*
when the origin is trustworthy. On NGFS-estimated ignitions it is not.

### What this does NOT establish

One fire, one observation time. The 0.011 IoU gap between arms is far smaller
than their common 0.67 shortfall, so this does **not** rank the two wind sources
— it says they are interchangeable at this fire's tolerance. Caveats on record:
arm B winds are hourly against arm A's 30-min; arm A is fire-coupled and arm B
diagnostic, so by 07Z (1537 ha) A has plume indraft and B does not; the
observation time is assumed UTC (02:02 local, consistent with a night IR flight).

---

## 8. HRRR forecast-cycle bug — two uncommitted edits

The HRRR-driven Dry River run failed at GRIB retrieval. Diagnosed by JH: **HRRR
forecasts should start at 0, 6, 12, 18 hours.**

`src/ingest/HRRR.py` had `cycle_hours = 1` with
`grib_forecast_hours_periods = [{'hours':48,'period':1}]` — an inconsistent pair,
because only the 6-hourly cycles reach f48. Confirmed against the bucket for
2026-09-15: t12z and t18z carry 49 files (f00-f48); t15z, t16z, t17z carry 19
(f00-f18). The job needed 26 hours from a 17Z cycle, so f19-f26 had never
existed. `aws` and the bucket were fine throughout; the error listing all 27
files as unavailable was misleading.

A second defect made it worse: the retry loop in `src/ingest/grib_forecast.py`
never stepped back, because `cycle_start` is non-`None` after the first pass and
the `if` branch is taken every time. It retried the identical impossible cycle
three times over 77 seconds. The "supports using previous cycles (up to 2)"
comment was false in practice.

**Edits made on explicit instruction, NOT committed, backups alongside:**

    src/ingest/HRRR.py            cycle_hours = 1 -> 6      HRRR.py.bak-20260918
    src/ingest/grib_forecast.py   remember the given cycle  grib_forecast.py.bak-20260918
                                  and shift back from it

Verified before relaunch: selection lands on 12Z needing f05..f31, retries step
to 06Z (f11..f37) and 00Z (f17..f43), all within f48; all 27 files present. The
relaunch used t12z, fetched 27 GRIBs (12 GB), built 27 of 27 COLMET files and
ran to `Completed` with 53 wrfouts.

Cost note: snapping to 6-hourly cycles trades up to 5 hours of data freshness.
For these jobs it costs nothing — every recent NGFS run is 26-27 hours, beyond
f18, so an hourly cycle could never have satisfied them. **If short-span HRRR
jobs are ever wanted, `cycle_forecast_hours` on `shell-refactoring` is the
finer-grained alternative** (48 for synoptic cycles, 18 otherwise); it exists on
no other branch.

## 9. A live trap in the `new_wrfxpy` install

`/data/jhaley/new_wrfxpy/wrfxpy/src` is a **real directory** with its own
`ingest/HRRR.py` and `ingest/grib_forecast.py`. Only `forecast.py` is a symlink
to `/data/jhaley/wrfxpy/src/forecast.py`. Because Python resolves symlinks when
setting `sys.path[0]`, that one symlink pulls the whole import path into
`/data/jhaley/wrfxpy/src`, so **`new_wrfxpy`'s own copies are dead code**.
Verified by loading the modules the way the launcher does:

    sys.path[0]      : /data/jhaley/wrfxpy/src
    HRRR module file : /data/jhaley/wrfxpy/src/ingest/HRRR.py

Editing `new_wrfxpy/wrfxpy/src/ingest/*.py` has no effect at all, silently. The
two copies already differ (`remote_url` has 3 entries there against 2 in the live
file), which is how the live file was identified from the log.

---

## 10. Open items

1. **The two source edits are uncommitted** (§8). Backups alongside. Not
   committed per repo policy.
2. **Nothing prunes `ingest/HRRRA`** (42 GB) or the new `ingest/HRRR` (12 GB).
3. **`.ff` templates still hardcode `fuelstrans.csv`** — 09-11 §9 item 8, still
   open, unchanged by this session. Both arms used the same table so it cancels
   in §7, but an operational run still does not use the derived table.
4. **Shape, not wind, is where the ForeFire gap lives** (§7). 09-11 §9 item 14
   is partly answered: inputs are ruled out, `pSAF` is a size knob, and now the
   wind field is ruled out for shape. What remains is the missing process —
   spotting — and the anisotropy a uniform speed multiplier cannot reproduce.
5. **The Ranger Road three-cell experiment** (09-11 §9 item 16) was not run.
6. WindNinja hourly fields and all comparison scripts lived in the session
   scratchpad and are gone. `dry_ff_arms.py` (two-arm runner, monkey-patches
   `make_FF_nc`) and `dry_ff_compare.py` (perimeter scoring) are the two worth
   rewriting into the repo if this is repeated.
7. ForeFire outputs are kept at `forefire_A/` and `forefire_B/` beside the HRRR
   Dry River workspace, KMLs suffixed `_A`/`_B` so they open together.

---

## 11. NEXT SESSION — the goal JH set

**Stand up something that can run a set of ForeFire forecasts from only an
ignition point, an ignition time, and GRIB files.** No WRF-SFIRE in the loop.

Why this is now reachable, and what it needs:

- **Winds are solved.** WindNinja ingests raw HRRR GRIBs directly (§1), costs
  ~56 s for a 30 km domain at 50 m on flat terrain (§3), and produces a wind
  field that is interchangeable with coupled `UF`/`VF` at this fire's tolerance
  (§7). That was the hard dependency.
- **What the FF netcdf still takes from the wrfout**: `make_FF_nc` reads
  `NFUEL_CAT` and `ZSF` from `wrfinput_d01` and `UF`/`VF` from the wrfout. Fuels
  can come from the LANDFIRE mosaic and terrain from `LF2020_Elev_220_CONUS`
  directly — both are already staged and the DEM matches `ZSF` to RMS 2.2 m (§1).
  So the netcdf can be built without WRF at all.
- **`windrf` must be applied** to the WindNinja field per fuel cell (§5), which
  requires the fuel map anyway.
- **The `min_timesteps` gate** (`run_forecasts`, 20 steps) exists to wait for
  WRF and is a deliberate temporary measure (09-08 §9). It is the thing to drop
  when this path is built — not before.
- **The timing table** is currently derived from wrfout filenames
  (`make_timing_table`). A GRIB-driven path needs the equivalent from GRIB valid
  times.
- **Set expectations from §7**: such a forecast should land near ForeFire's
  current skill — area within ~0.85x of observed, IoU ~0.34, and elongation too
  round. It will not be better than the coupled run, and on this evidence it will
  not be much worse.
