# ForeFire session handoff — 2026-09-29

Bigfoot end to end, a rework of the apptainer binding, and two findings about fire
breaks — one of which is a deliberate choice rather than a bug.

## 1. TL;DR

- **`wksp_forefire/<grid_code>/` is the new home for ForeFire-only forecasts**, bound
  per job at `/work`. The shared run directory is gone, and with it the collision class
  the cron `flock` was working around (`d814825`).
- **Bigfoot ran end to end** — a VIIRS-only fire with no WRF-SFIRE workspace, 27 steps,
  `--issue-utc` honoured. This closes the `--steps 27` unknown.
- **50 m is the timing sweet spot**: 20 min against 77 at 25 m, for 90% of the area.
- **But 50 m does not restore barrier continuity.** Longest unbroken non-burnable run is
  5.93 km at 50 m against 16.62 km at 25 m — essentially 100 m's 5.24 km. Area degrades
  smoothly with resolution; **connectivity collapses between 25 m and 50 m and then
  barely changes.**
- **Highways are erased on purpose.** `remap_fuel` fills LANDFIRE code 91 with the
  nearest burnable fuel. In Bigfoot's domain that is **13,359 cells against 2,443 cells
  of river kept** — 5.5x more barrier deleted than preserved, at full resolution. This,
  not mesh size, is why models jump highways.
- **Nothing has ever scored the ensemble members.** Every area number this project has
  produced reads the chained step progression; the 28 `final` perimeters are invisible
  to it.

## 2. Per-job binding (`d814825`)

ForeFire ran in one shared directory, and `run_forecasts` deletes every file in it on
entry, so a second invocation wiped a running one's netcdfs — 26 of 47 Union_400257
members came out empty that way. Each job is now bound separately and the shared
directory is gone; the `flock` is belt-and-braces rather than load-bearing, and fires
could now run concurrently.

**The blocker was apptainer's bind syntax.** Bind specs are colon-delimited and cannot
carry a colon in a path, and a wrfxpy `wfc-` workspace name embeds the ESMF timestamp:

    bad bind syntax ".../wfc-Bigfoot_...-2026-09-28_18:00:00-33/forefire/:"

`wksp_forefire/<grid_code>` avoids it because `grid_code` uses underscores. A fallback
binds the nearest colon-free ancestor where a path does contain one, with `--pwd`
supplying the isolation.

**`grid_code` is sanitised for use as a directory name**, because over 28,695 real job
files **7,821 (27%) contain a colon** from the older `18:00:00` convention, plus spaces
(10), `#` (22), parentheses (9), `~` (5) and `&` (1). Mapping `:` to `_` is what the
current convention already does, so old and new jobs converge on one name.

Also replaced `cleanup_ff_run`'s `os.system(f"mv {dir}/* {dest}/.")` with `shutil.move`:
unquoted, the grid_codes containing a space or an `&` do something other than a move.

Verified behaviour-neutral: 128.3 ha against 128.5 for the same configuration under the
old bind, 0.16% — inside the ~3% run-to-run noise.

## 3. Resolution: timing, area and barriers are three different questions

Bigfoot, 27 steps, identical configuration:

| res | cells | netcdf build | s/netcdf | longest barrier | final area |
|---|---|---|---|---|---|
| 25 m | 1,537,600 | **77 min** | 172 | **16.62 km** | 214.7 ha |
| **50 m** | 384,400 | **20 min** | 45 | 5.93 km | **192.9 ha** |
| 100 m | 96,100 | 11 min | 25 | 5.24 km | 128.5 ha |

**Timing.** 50 m is 3.8x faster than 25 m for 90% of the area. Going on to 100 m saves
only 9 more minutes and costs a further 33% of the fire. Cost scales as `cells^0.93`;
area does not.

**Area** degrades smoothly: 215 → 193 → 129 ha.

**Connectivity does not.** It collapses between 25 m and 50 m and then flattens. A river
about 100 m wide is 4 cells at 25 m, 2 at 50 m, 1 at 100 m, and stops being a reliable
barrier below roughly 3 cells. **Total barrier *area* is preserved** (861 ha at 25 m,
842 at 100 m) — the resampling conserves area and destroys topology, so an area-based
check on the fuel map would have called it fine.

JH observed the consequence directly: at 100 m the **fourth** full-duration member jumped
the Trinity River. Measured, that member is **379.6 ha against 169 and 156 for its
neighbours** — 2.2x, with the widest east-west span. Otherwise the river held.

**So mesh choice is not one decision.** If a barrier matters for a fire, only 25 m keeps
it; 50 m does not buy protection. If it does not, 50 m gives 90% of the answer in a
quarter of the time.

**And mesh is coupled to pSAF**, which was calibrated entirely at 25 m. Moving the
default to 50 m implies roughly `0.337 * (214.7/192.9)^(1/1.54) = 0.36`; to 100 m, 0.47.

## 4. The bigger fire-break problem: highways are deleted on purpose

Per JH, models jump highways when they should not. That is not the mesh. `remap_fuel`
fills LANDFIRE codes 0, 91 and 93 with the nearest burnable fuel, with the reason stated
in its own docstring — *"so a fire is not stopped by a road or a field boundary that
LANDFIRE happens to label unburnable"* — mimicking geogrid's nearest fill.

Raw LANDFIRE in Bigfoot's domain:

| code | meaning | cells | fate |
|---|---|---|---|
| **91** | **urban/developed, includes roads** | **13,359** | **filled with fuel** |
| 98 | water (rivers) | 2,443 | kept as barrier |
| 99 | barren | 953 | kept as barrier |
| 93 | agriculture | 1,437 | filled with fuel |

**5.5x more barrier is erased than kept**, at full resolution, before any coarsening.
Rivers survive because water is code 98; highways do not because they are code 91.

The choice is defensible for a subdivision or a dirt road — fire does run through those —
and wrong for an interstate. **LANDFIRE code 91 cannot tell them apart**, so fixing this
properly needs an external roads layer with a class attribute (TIGER/Line, or LANDFIRE's
own roads product) to keep major highways while continuing to fill developed land.

### JH's buffering idea, and a cheaper alternative

Buffering linear barriers once on the LANDFIRE tif would preserve them through
resampling. Two cautions and one alternative worth weighing first:

- **Buffering enough to survive 100 m means a 300 m wide river**, which over-blocks badly.
  The width needed scales with the target mesh, so one buffered raster cannot serve
  several resolutions.
- **The resampling rule is the cheaper lever.** Fuel is warped with `-r near`
  (`forefire_grib.py:515`), so a 30 m river inside a 100 m cell survives only if the cell
  centre happens to land on it — which is exactly how a continuous river becomes a dashed
  one. Warping a **binary barrier mask separately with `-r max`** and OR-ing it back in
  would preserve connectivity at one cell wide, without touching LANDFIRE and without
  choosing a width per mesh. It over-represents barrier area slightly, which is the safe
  direction.
- Whichever route, **do it for roads and rivers together**, since §4 shows roads are the
  larger loss.

## 5. Nothing has ever scored the ensemble members

`ff_growth.track()` skips the `final` geojsons. Those 28 files are the time-lagged
ensemble's full-duration forecasts — the actual forecast product — and **every area
number this project has produced reads the chained step progression instead**: the pSAF
calibration, the 25/50/100 m comparison, the RTMA/HRRR work.

For calibrating against WRF-SFIRE's `TIGN_G` progression the chained track is the right
analogue, so that work stands. But the ensemble **spread** has never been examined, and
outliers like Bigfoot's member 4 — a 2.2x excursion from one barrier breach — are
invisible to everything currently in the repo.


## 5b. The `-r max` barrier mask works at 50 m (`0d1d2f1`)

Tried the cheaper lever from §4 rather than buffering the LANDFIRE tif. A second warp of
the same source with `-r max` marks any target cell containing a barrier cell.

It works on the raw codes because every genuine barrier (92, 94-99) sits numerically
above the two `remap_fuel` deliberately fills -- 91 urban/roads and 93 agriculture -- so
a barrier cannot be outvoted by surrounding fuel. Snow/ice at 92 is the one exception,
outranked by 93; rare in fire season.

| res | near only | **with max** | pieces | non-burnable |
|---|---|---|---|---|
| 25 m | 16.62 km | 16.72 km | 506 -> 381 | 0.90 -> 1.33% |
| **50 m** | **5.93 km** | **14.63 km** | 312 -> 182 | 0.88 -> 1.61% |
| 100 m | 5.24 km | 7.47 km | 171 -> 89 | 0.88 -> 1.94% |

**50 m with the mask recovers 88% of the 25 m barrier continuity at a quarter of the
cost.** 100 m still loses the river even with it, so the mask does not rescue 100 m.
The cost is over-stated barrier width, which is the safe direction.

### The ensemble is where it shows, and the chained track barely notices

| member | 100 m no mask | 50 m no mask | **50 m + mask** |
|---|---|---|---|
| 3 | 169.2 | 286.1 | 207.8 |
| **4** (the river jumper) | **379.6** | 282.4 | **198.5** |
| 5 | 156.3 | 253.4 | 182.4 |

Member 4 goes from 2.2x its neighbours to *below* member 3. Ensemble **max/median falls
3.11 -> 1.66 -> 1.50**. The chained progression moves only 192.9 -> 161.3 ha and never
flags anything -- which is §5's point restated: the diagnostic that catches this is the
one nothing in the repo computes.

**Recommended configuration: 50 m with the barrier mask.**

**It does not touch highways.** The mask restores only barriers `remap_fuel` was already
keeping. Code 91 is still filled by design, 13,359 cells against 2,443 of river in this
domain, so §4 stands unaddressed.

**It moves the calibration again.** pSAF 0.337 was fitted at 25 m with no mask; 50 m +
mask gives 161.3 ha where 25 m gave 214.7, implying roughly
`0.337 * (214.7/161.3)^(1/1.54) = 0.41`. Mesh, mask and pSAF are now a three-way coupled
choice, and the honest way to settle it is against observed perimeters rather than
against each other.

## 6. Open items

1. ~~**Barrier preservation** -- the `-r max` mask~~ **DONE (§5b), for rivers.** Still
   open for **roads**: whether to stop erasing code 91, which needs an external roads
   layer to separate an interstate from a subdivision. That is the larger loss.
2. **Score the ensemble members** (§5). At minimum report the spread and flag outliers;
   a member that is 2x its neighbours is a barrier breach or a bad wind field.
3. **Pick a default mesh and mask, then re-derive pSAF** (§3, §5b). 50 m + mask is the
   candidate, implying pSAF ~0.41. Settle it against observed perimeters: three
   configurations that differ by 25% in area cannot be ranked against each other.
4. Carried from 09-28: `--issue-utc` on every forecast run; perimeter scoring has its own
   calibration (score at containment, IoU 0.25 is a success, check detection dates);
   `ff_domain_connectivity.py` before any per-fire fit; re-run `fmda_fcst_skill.py` once
   the cron has accumulated cycles.

## 7. NEXT SESSION

WRF-SFIRE is running Bigfoot for comparison (per JH). Once it lands, that is the first
genuine ForeFire-vs-WRF-SFIRE test on a fire neither was tuned to — and Bigfoot has a
river barrier, which makes it a useful case for §4.

Commits remain unpushed on `james_ngfs` and `release-fmda-fixes`; JH pushes them manually
(SSH password).
