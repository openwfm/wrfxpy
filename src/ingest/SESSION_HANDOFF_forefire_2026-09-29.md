# ForeFire session handoff — 2026-09-29

Bigfoot end to end, a rework of the apptainer binding, and two findings about fire
breaks — one of which is a deliberate choice rather than a bug.

## 1. TL;DR

- **`wksp_forefire/<grid_code>/` is the new home for ForeFire-only forecasts**, bound
  per job at `/work`. The shared run directory is gone, and with it the collision class
  the cron `flock` was working around (`d814825`).
- **Bigfoot ran end to end** — a VIIRS-only fire with no WRF-SFIRE workspace, 27 steps,
  `--issue-utc` honoured. This closes the `--steps 27` unknown.
- **50 m with rivers masked is the chosen configuration**, per JH: 20 min against 77 at
  25 m, and the `-r max` pass recovers 88% of the 25 m barrier continuity (§5b).
- **pSAF stays at 0.337.** 0.41 was refuted over 20 forward forecasts (§5d), and 0.337
  is independently the best value at the 8-hour horizon: 78% within a factor of two.
- **Skill degrades with lead** -- geo mean 1.51 at 2 h, 1.92 at 8 h, 2.35 at 18 h -- so
  **pSAF has to be fitted at the horizon the product is issued at**, which 0.337 was not.
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

**Decided, per JH: 50 m with rivers masked this way is the resolution to aim for.**

**It does not touch highways.** The mask restores only barriers `remap_fuel` was already
keeping. Code 91 is still filled by design, 13,359 cells against 2,443 of river in this
domain, so §4 stands unaddressed.

**It moves the calibration again.** pSAF 0.337 was fitted at 25 m with no mask; 50 m +
mask gives 161.3 ha where 25 m gave 214.7, implying roughly
`0.337 * (214.7/161.3)^(1/1.54) = 0.41`. Mesh, mask and pSAF are now a three-way coupled
choice, and the honest way to settle it is against observed perimeters rather than
against each other.


## 5c. The overnight cron runs: identical parameters, size-dependent skill (2026-09-30)

JH asked what the ~13 well-matching cron runs had in common. **Nothing in their
parameters** -- all 21 scored runs are identical:

    fire mesh 25 m, domain 30 km, pSAF 0.337, windReductionFactor 0.4,
    perimeterResolution 25, and 25-28 per-step moisture tables each

So the spread from 0.18x to 4.71x is entirely fire-to-fire. What separates them is
**size**:

| WRF-SFIRE size | n | geo mean | median |
|---|---|---|---|
| <100 ha | 6 | **2.37** | 2.60 |
| **100-500 ha** | 10 | **1.10** | 1.27 |
| >500 ha | 5 | 1.43 | 1.34 |

`ratio ~ area^-0.155`, log-log r=-0.356 -- consistent with the -0.198 on the larger
batch. Five of the six worst are under 156 ha, four under 100. Overall n=21, geo mean
1.46, median 1.41, **13 of 21 within a factor of two**, which is JH's 13.

**pSAF 0.337 is well tuned for 100-500 ha and overpredicts small fires by ~2.4x.** That
is the clearest statement yet of why a single constant cannot serve the whole range, and
it is the case for the per-fire pSAF work (09-24 handoff §9 item 2).

**A correction to the plan in §7.** Because every run already shares one parameter set,
testing pSAF 0.41 on these 25 m runs would simply scale all 21 ratios by ~1.21x -- geo
mean 1.46 -> ~1.77, with *fewer* fires inside a factor of two. The 0.41 figure exists to
compensate for the area lost at 50 m with the barrier mask, so it is only meaningful
applied **together with** that configuration, never on its own.


## 5d. The forward-forecast batch: pSAF 0.41 refuted, and calibrate at the horizon you issue (2026-09-30)

20 of the 21 overnight fires re-run as **genuine forecasts** -- `forefire_forecast.py`
from each job file, 50 m with the barrier mask, pSAF 0.41, and `--issue-utc` set an hour
after each ignition so **both winds and fuel moisture** are restricted to what was
published by then. Mill_Creek_RX is a 2025 fire and predates the forecast cache.

**`--issue-utc` had to be extended to the moisture first.** It only filtered the GRIBs;
`read_fmda_hourly` had no issue filter and *prefers analyses over forecasts*, so a
forecast run silently took moisture that assimilated observations made after the forecast
was issued -- exactly what the flag exists to prevent. `issue_utc`/`issue_lag_h` now run
through `read_fmda_hourly` -> `md_series` -> `make_script_set`.

### pSAF 0.41 is refuted

| | n | geo mean | median | within 2x |
|---|---|---|---|---|
| forward, 50 m, pSAF 0.41 | 20 | **1.87** | 2.01 | **45%** |
| cron, 25 m, pSAF 0.337 | 20 | 1.41 | 1.40 | 65% |

**17 of 20 fires rose, median +33%**, against an intended ~0%. The 0.41 came from
Bigfoot alone -- 50 m + mask gave 161.3 ha where 25 m gave 214.7, a 25% loss -- and that
single fire's barrier effect does not generalise. **Reverted to 0.337.**

**The batch also changed three things at once** -- mesh, pSAF and forward-vs-analysis --
so it cannot attribute the shift. The clean follow-up is the same 20 fires at 50 m with
pSAF 0.337 and analysis moisture: one variable from the cron baseline.

### Skill degrades with lead, so calibrate at the horizon you issue

| lead | n | geo mean | median | within 2x |
|---|---|---|---|---|
| 2 h | 15 | **1.51** | 1.38 | **87%** |
| 4 h | 16 | 1.76 | 1.59 | 62% |
| **8 h** | 18 | **1.92** | 1.76 | **67%** |
| 12 h | 18 | 2.27 | 1.99 | 50% |
| 18 h | 18 | 2.35 | 2.20 | 39% |
| 24 h | 20 | 1.90 | 1.95 | 45% |

ForeFire's excess spread **compounds** rather than staying proportional -- Twin_Pine_8
runs 0.91 -> 1.52 -> 2.23, Rx_Angelina 1.79 -> 3.98 -> 4.55. A rate multiplier cannot
correct an error that grows with lead, so **pSAF must be fitted at the horizon the
product is issued at**. The 0.337 in use was fitted at the run endpoint, 24-27 h.

### What that means for the 8-hour product

| pSAF at 8 h | geo mean | median | within 2x | below 1x |
|---|---|---|---|---|
| 0.41 | 1.92 | 1.76 | 67% | 11% |
| **0.337** | **1.42** | **1.30** | **78%** | 17% |
| 0.30 | 1.19 | 1.09 | 72% | 39% |
| 0.26 | 0.95 | 0.87 | 83% | 67% |

**0.337 is the right value for an 8-hour product** on JH's standing posture: best
within-2x at 78%, median 1.30, only 17% underpredicting. 0.26 buys 83% within 2x at the
cost of flipping two thirds of fires to underprediction.

**A reassuring accident:** the 25th-percentile posture that set 0.337 at the 24 h horizon
wants **0.312** at 8 h -- within 8% of the value already in use. Centring the 8 h *median*
would want 0.284, so 0.337 sits deliberately on the overprediction side.

### Two caveats on these numbers

- **Seven of 18 fires had degraded moisture** -- usable FMDA for only 4-8 of 30 steps,
  the rest falling back to the base table at Md=0.1. Four are among the worst ratios
  (Pulaski 5.26, Rx_Angelina 3.98, Comp_265 2.60, Twin_Pine_8 1.52 at 8 h). FMDA forecast
  cycles only start 2026-09-28, so early fires had nothing to draw on. True 8 h skill is
  probably better than shown.
- **The small-fire bias is not a lead-time effect.** Pulaski (25.6 ha) sits at 5.26 and
  Paula_Red_South (12.4 ha) at 4.49 at every horizon. No single pSAF addresses
  `ratio ~ area^-0.2`.


## 5e. Parameter sweeps are ~30x cheaper than full runs (`0278a20`)

Per JH, most of the CPU is WindNinja. Measured over all 21 forward runs it is
**82-98% of the wall clock, median 95%** -- 19 to 46 minutes of netcdf building against
**1 to 6 minutes of ForeFire**. The netcdfs carry winds, fuels and terrain, none of which
depend on the tuning parameters.

`sweep_job()` in `forefire_forecast.py` builds them once and varies only the `.ff`
settings. **Measured on Finley_RX: three pSAF values in 3.1 minutes, against 34 minutes
for its original single run.** A 20-fire by 5-value sweep is therefore about 3 hours, not
two days.

Each variant runs in its own subdirectory with the netcdfs symlinked in, because the
restart template chains through `include[sim_<grid>_<step>.ff]` and two variants sharing
a directory would read each other's perimeters. `run_timing_table` gained `bind_dir` to
make that work under the per-job bind: the links are relative and resolve only if the
**parent** is the mount, so `bind_dir` is the base directory while `--pwd` selects the
variant.

### A caution before building a pSAF ensemble

Finley_RX gives **96.2 ha at pSAF 0.337 and 95.4 at 0.41** -- indistinguishable across a
22% change -- then collapses to **6.8 ha at 0.26**, with only 13 of 28 perimeters.

Its moisture runs 0.103-0.237 against fuel 2's extinction of 0.15, so it is
**moisture-limited, not rate-limited**. pSAF has almost no leverage until it drops far
enough to extinguish the fire outright. An ensemble built by perturbing pSAF would
produce near-zero spread on exactly the fires that are hardest to forecast, and a cliff
rather than a gradient when it does respond.

**So a useful ensemble has to perturb what each fire is actually sensitive to.** Per JH
the candidates are pSAF, ignition location, fuel moisture and other inputs. On this
evidence **moisture is the stronger lever on moisture-limited fires** and pSAF on
rate-limited ones, and which regime a fire is in can be read off in advance by comparing
its Md series against each fuel's moisture of extinction -- the §7 method from the 09-28
handoff, already implemented as `ff_moisture_stall_test.py`. Choosing the perturbation
per fire would give a better-conditioned ensemble than perturbing everything everywhere.

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

**Isolate the mesh.** §5d changed three things at once -- mesh, pSAF and
forward-vs-analysis -- so nothing in it attributes the shift. Re-run the same 20 fires at
**50 m with pSAF 0.337 and analysis moisture**: one variable from the 25 m cron baseline
(geo mean 1.41, median 1.40, 65% within 2x). The fire list is
`fwd_list.txt` in the session scratchpad; regenerate it from `ff_batch_ratio.py 6` if
that is gone.

**Then re-check the 8-hour numbers once FMDA forecast coverage is continuous.** Seven of
18 fires in §5d ran most of their steps on the base fuel table because forecast cycles
only start 2026-09-28, and four of those are among the worst ratios. The cron now runs
four cycles a day, so a week from now that confound is gone.

**Score the ensemble members, not only the chained track** (§5, §5b). Still nothing in
the repo does this, and it is what caught Bigfoot's river crossing.

**Build the perturbed ensemble** (§5e). `sweep_job()` makes it cheap -- about a minute a
variant. Choose the perturbation per fire rather than uniformly: moisture for
moisture-limited fires, pSAF for rate-limited ones, with `ff_moisture_stall_test.py`
deciding which is which. Ignition location is the third candidate and is untested.

**The remaining structural problems, in order of size:**

1. **Highways are erased by design** (§4). `remap_fuel` fills LANDFIRE code 91 with the
   nearest burnable fuel -- 13,359 cells against 2,443 of river kept in Bigfoot's domain.
   The `-r max` mask does not touch it. Needs an external roads layer to separate an
   interstate from a subdivision.
2. **The small-fire bias is not a lead-time effect** (§5d). Pulaski at 25.6 ha sits above
   5x at every horizon. `ratio ~ area^-0.2` needs a per-fire correction, which is the
   standing NGFS detection-field idea (09-24 handoff §9 item 2); run
   `ff_domain_connectivity.py` first and exclude the flagged fires.
3. **Error compounds with lead** (§5d), so a single rate multiplier cannot serve both an
   8-hour product and a 24-hour one. If both are wanted, they need separate values.

**Against observed perimeters, not against WRF-SFIRE.** Everything above is model-vs-model
and WRF-SFIRE is not truth -- on Dome it ran 0.68x an IR perimeter. Score at containment
time, IoU 0.25 is a success, and check detection dates first.

Commits remain unpushed on `james_ngfs` and `release-fmda-fixes`; JH pushes them manually
(SSH password).
