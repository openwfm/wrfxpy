# ForeFire session handoff — 2026-09-09

Branch `james_ngfs`, HEAD `60ce1fa`. Nothing committed.
Read `SESSION_HANDOFF_forefire_2026-09-08.md` first — it carries the module
description, the config, and the traps. This file covers only what is new today.

**All runs described here completed.** 29-set sweep: 1566/1566 steps, zero
errors, 187 min. Resolution study: 8 sets. pSAF x Md grid: 25 sets, 375 steps,
13.3 min. Figures are in `<wksp>/forefire/`.

---

## 1. Observation: `forefire.py` changed outside this session

`src/ingest/forefire.py` has mtime **2026-09-09 05:45:38**, later than the 09-08
handoff (2026-09-08 11:26:21). The file is untracked, so there is no diff to
inspect and no way to say what changed. Verified unchanged in structure: still
1366 lines, 30 functions, and `read_fmda_fmc`, `resolve_md`, `write_fuel_table`,
`moisture_params`, `sweep_ff_params`, `finalize_ff_outputs` all present. Nothing
in this session edited it. **Recorded, not explained.**

---

## 1b. Deleted: `src/ingest/wrfin2forefire.py`

Removed this session, confirmed as scratch work. It held one function,
`convert_nc(nc_input_file)`, the direct ancestor of `make_FF_nc` — same
projection maths, same subgrid slicing, same comment markers. `make_FF_nc` is a
strict superset; the whole delta was:

| | `convert_nc` | `make_FF_nc` |
|---|---|---|
| signature | `(nc_input_file)` | `(nc_path, out_path)` |
| output | hardcoded `forefire_dataset.nc` in cwd | caller-specified |
| **fuel + terrain read from** | **the input file itself** | **`wrfinput_d01` beside it** |
| dead code | two commented-out constant-wind lines | removed |

That third row is why it could no longer work at all. On current output the four
required variables are split across two files:

```
wrfinput_d01   has NFUEL_CAT, ZSF   missing UF, VF
saveout        has UF, VF           missing NFUEL_CAT, ZSF
```

`convert_nc` reads all four from the single file it is handed, so it raises
`KeyError` on either. Its `__main__` pointed at a 2023 TITAN `wrfout_d01`, which
presumably carried all four — likely a different `iofields` configuration then.

Nothing referenced it (checked `src/`, `etc/`, `jobs/`, top-level scripts). It was
untracked, so **the deletion is not recoverable from git** — the delta above is
the record. A session-local copy exists at
`$SP/wrfin2forefire.py.deleted-2026-09-09` but the scratchpad does not survive.

---

## 2. New sweep: SINLAHEKIN fire

Workspace:
`/data/jhaley/wrfxpy/wksp/wfc-SINLAHEKIN_2026-07-26_20_00_00_78D35D3B-F791-4961-AE36-C6D1A4DFF5A0-2026-07-26_18:00:00-30`

Chosen because Ranger Road was an exceptional, very-high-rate-of-spread plains
fire, and this one was observed extensively — observational data is expected
later, so it is the more useful calibration target.

|  | Ranger Road | SINLAHEKIN |
|---|---|---|
| runnable steps | 20 | **54** |
| horizon | 12 h | **27 h** |
| fire grid | 2400², 50 m | 1200², 25 m |
| baseline area (defaults) | 752 ha | **2 606 ha** |
| dominant fuel | 93.6% class 2 | mixed (below) |

A completed baseline run already existed in `<wksp>/forefire/` (54 netcdfs, 108
geojsons, loose in that directory). The sweep writes only into `<tag>/`
subdirectories, so the baseline is untouched.

**Fuel mix is the structural difference.** Ranger Road had one extinction cliff
at `me = 0.15`. SINLAHEKIN is mixed and extinguishes as a staircase:

| fuel | share | `me` |
|---|---|---|
| 5 brush | 32.6% | 0.20 |
| 9 hardwood litter | 20.8% | 0.25 |
| 2 timber grass/understory | 16.9% | 0.15 |
| 10 timber litter + understory | 14.4% | 0.25 |
| 8 closed timber litter | 13.6% | 0.30 |

Nothing burns above Md = 0.30.

**29 sets in three phases**, driver at
`$SP/sinlahekin_sweep.py` (scratch — see §6):

1. `windReductionFactor` — 0.2 … 0.8 (7). Same list as Ranger Road, deliberately,
   so the two fires compare directly.
2. `propagationSpeedAdjustmentFactor` — 0.3, 0.45, 0.6, 0.8, 1.0, 1.2 (6) plus 2
   ROS-matched probes whose pSAF is derived from phase 1's spread rates. Re-tests
   the Ranger Road finding that the two knobs are not interchangeable.
3. Dead fuel moisture `Md` — 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11,
   0.13, 0.15, 0.18, 0.22, 0.28 (14). **Weighted below 0.15**: fire-weather
   conditions are the interest, and the three tail values only need to establish
   whether the fire still moves above 0.15.

FMDA moisture is left **disabled** throughout, so phases 1–2 sit at the same
template 10% the Ranger Road sweeps used and the cross-fire comparison is not
confounded.

Measured cost: 4.5 min for the first set (54/54 steps, no errors); ~2–2.5 h total.

---

## 3. Prior expectation on record

Before phase 3 ran, the expectation stated was that **the fire would hardly move
above Md = 0.15**. The fuel table says it should still move: only the 16.9%
timber-grass class extinguishes at 0.15, while brush, litter and closed timber
carry to 0.20–0.30. But `Etam` damping is severe as `Md` approaches `me`, so a
sharp drop rather than a stop is the competing expectation. The three tail values
(0.18, 0.22, 0.28) decide it.

**Outcome: the stated expectation was closer.** The fire reaches 515 ha at Md 0.18
(15% of baseline) and **exactly zero at 0.22 and 0.28**, despite 65% of the domain
being nominally above its extinction moisture at 0.22. The fuel-table `me` values
overstate how far the fire actually carries — either the surviving fuels are not
contiguous enough to sustain a front, or `Etam` damping drops spread below
`minSpeed` before formal extinction. The fuel-table reasoning was wrong.

---

## 4. Terrain and wind field: the two fires are physically different

Computed from the `windU`/`windV` and `altitude` variables in the FF netcdfs
themselves (8 time steps sampled per fire), not assumed.

| | Ranger Road (plains) | SINLAHEKIN (complex) |
|---|---|---|
| terrain relief | 368 m | **2 123 m** |
| elevation std | 82 m | 449 m |
| mean slope | 1.4° | **16.1°** (p95 33.6°) |
| mean wind speed | 3.2 m/s | 1.1 m/s |
| direction coherence *R* | 0.993 | **0.748** |
| circular std of direction | **6.2°** | **41.5°** |

*R* is the resultant length of the unit wind vectors over the domain: 1.0 is a
single uniform direction. Ranger Road is effectively one wind vector across
120 km. SINLAHEKIN's directions are scattered over ±40°.

---

## 5. Why `windReductionFactor` has far less leverage on SINLAHEKIN

This is the substantive result of the session so far.

Rothermel computes `R = R0 · (1 + phiV + phiP)` (`Rothermel.cpp`), where `phiV`
is the wind term — **the only place `windReductionFactor` acts** — and `phiP` is
the slope term, which wRF cannot touch:

```
phiV = C · (Beta/Betaop)^-E · normal_wind^B      normal_wind already × wRF
phiP = 5.275 · Beta^-0.3 · tan(slope)²
```

Evaluated with the actual fuel-table parameters, each fire's dominant fuel, its
measured mean wind and mean slope, at wRF = 0.4:

| fire | fuel | `phiV` | `phiP` | slope share | `d(ln R)/d(ln wRF)` |
|---|---|---|---|---|---|
| Ranger Road | 2 timber grass | **5.94** | 0.02 | 0% | **1.300** |
| SINLAHEKIN | 5 brush | 1.39 | **2.66** | **53%** | **0.383** |

On Ranger Road the wind term is essentially the whole of the spread enhancement.
On SINLAHEKIN **the slope term is larger than the wind term** — over half the
spread comes from a term wRF has no access to — so the knob has **3.4× less
leverage**.

Two independent supports:

- The predicted Ranger Road elasticity of 1.300 matches the **measured**
  `ROS ~ wRF^1.22` from the 2026-09-08 sweep.
- **The completed sweep measured `ROS ~ wRF^0.40` against the predicted 0.383** —
  an out-of-sample hit, since the prediction came from the wind field and terrain
  alone before any of these runs existed. Area scales as `wRF^0.61` here versus
  `wRF^1.59` on Ranger Road; over the full 0.2–0.8 range the fire grows 2.29x
  versus 8.87x on the plains fire.

**Caveats on the numbers.** `phiV`/`phiP` are evaluated at a single
representative point (domain-mean wind, domain-mean slope, dominant fuel), not
integrated over the domain, so treat 3.4× as an order-of-magnitude statement
rather than a coefficient. It is also an *under*-statement: Rothermel uses the
wind component **normal to the fire front**, and with directions scattered over
±40° the effective normal wind is further reduced, so real leverage is likely
below 0.383.

### Consequence for the calibration plan

**On complex-terrain fires, inverting observed area back to wRF is poorly
conditioned.** The Ranger Road sweep produced a clean monotonic calibration curve
where ±25% area error pinned wRF to ±0.05. Here the same knob barely moves the
answer, so the same observation constrains it far more weakly. Do not assume the
Ranger Road inversion precision transfers.

### A falsifiable prediction, recorded before phase 3 ran

Moisture should **not** lose leverage the same way. `Md` acts on `R0` through
`Etam` and `Qig`, and `R0` multiplies the whole `(1 + phiV + phiP)` bracket — the
slope term dilutes wRF but not moisture. So:

> the Md sensitivity curves should transfer between the two fires far better than
> the wRF curves do, and `cross_fire_sensitivity.png` should show the wRF panels
> diverging sharply while the Md panels stay comparable.

If that fails, the framework above is wrong and should be revisited before any of
it is used for calibration.

---

## 5b. Spotting: the model has it, but not in the way we run it

Investigated because ember transport matters most for exactly the fast
wind-driven plains case. **ForeFire's spotting is a coupled mechanism split into
three stages, and only the first is reachable from a standalone run.**

**Stage 1 — emission. Present.** `src/flux/SpottingFluxBasicModel.cpp` registers
as flux model `"SpottingFluxBasic"` (`:62`). Each burning cell emits a firebrand
flux in kg m-2 s-1, constant for `spottingDuration` seconds after arrival time,
scaled by a per-fuel coefficient `fuel.Spot` (`:78`). Tunables, with defaults:
`spottingDuration` 100 s (`:89`), `nominalSpottingFlux` 1.0 (`:93`),
`lagSpotting` 0.0 (`:96`).

**Stage 2 — transport. Not in ForeFire at all.** The flux is meant to be handed
to a coupled atmospheric model, which lofts and advects the brands and returns a
ground-deposition field named `spotAtGround`.

**Stage 3 — re-ignition. Present, but only in the coupling library.**
`src/CLibForeFire.cpp` consumes `spotAtGround` (`:507`, `:523`). Where deposition
exceeds `reignitionThresholdFromSpotting` (`:526`) — and is below twice it, an
upper gate that skips heavily loaded cells — it draws a uniform random number
against `reignitionMaxProbalilityValueFromSpotting` (`:528`) and on success issues
a literal `startFire[loc=(x,y,0)]` at a random position inside the cell (`:560`).
Spot fires therefore arrive through the same multi-front point-ignition path
described in §5c.

### Two concrete blockers

1. **No transport in the standalone configuration.** Runs here are
   `forefire -i script.ff` inside apptainer against a static netcdf. Nothing
   lofts or advects brands, so stage 2 does not exist and stage 3 never fires.
   The re-ignition code lives in `CLibForeFire`, which is driven by an
   atmospheric model pushing `spotAtGround` in.
2. **No `Spot` column in either fuel table.** `SpottingFluxBasicModel` registers
   `fuel.Spot`, and neither `/home/jhaley/forefire/tests/fuelstrans.csv` nor
   `fuels13.csv` has it. Same failure mode as the `me` incident on 2026-09-08:
   the per-fuel property lookup fails.

The `.ff` templates set no flux model at all (only `propagationModel=Rothermel`),
so none of this is active today.

### Why it matters for comparing against detections

Spotting is the mechanism most likely to matter for fast wind-driven plains runs,
where downwind spot fires outpace the flaming front — and it is **structurally
absent** from these runs. So a detection cluster ahead of the modelled front may
be a real spot fire the model cannot produce, not a parameter error. **Tuning wRF
upward to "catch" spotting would be fitting the wrong knob**, and would corrupt
the wRF calibration for the cases where it is genuinely identifiable.

The planned multi-point initialization sidesteps this empirically: a spot fire
that appears in VIIRS/GOES becomes an ignition point on the next cycle, so
observed spotting enters the forecast without the model having to predict it.
For low-latency work that is arguably the better trade. The physical path would
mean coupling to WRF through `CLibForeFire` instead of the standalone binary —
a much larger change, and it needs a rebuilt container (see 2026-09-08 §1 on why
that is not local work).

---

## 5c. Multi-point ignition is supported — via repeated `lonlat`, not `points`

Checked in `Command::startFire` (`src/Command.cpp:249-671`). The input types
behave differently and the difference matters:

| type | behaviour |
|---|---|
| `points`, `kml`, `polyencoded` | `getPoly()` returns **one** polygon — a single perimeter, not several ignitions |
| `geojson` | parses multiple polygons, but nests polygon *p* inside polygon *p-1* (`:472-488`) — suits a perimeter with holes, not scattered detections |
| `lonlat` / `loc` | one small triangle around a point; **this is the one to repeat** |

Each `startFire[lonlat=...]` finalises the current front, then attaches a **new
front to the containing front** via `addFireFront(t, contfront)` — siblings, not
nested. So N such lines give N independent fronts that evolve separately and
merge naturally on contact. That is the mechanism for detection-driven
initialization; the template needs one line per detection, all at the same `t`.

Two cautions:

- Each ignition is guarded by `striclyWithinDomain(pos) & !isBurnt(pos,t)`, so
  detections outside the domain or on already-burnt ground are **silently
  skipped**. Count accepted versus supplied.
- The triangle is sized at `2 x perimeterResolution` (currently 25 m), so
  detections closer than ~50 m produce overlapping initial fronts.

On the wrfxpy side `read_input` reads only `cfg['ignitions']['1'][0]`, so
multi-point would need that widened too.

---

## 5d. Fuel moisture: the uncertainty dominates, the size-class weighting does not

### The finding that matters: two moisture sources disagree by 0.025

For SINLAHEKIN, the same system produces two different answers for the same fire:

| source | 1h dead moisture |
|---|---|
| FMDA geogrid at ignition (10 km around the point) | **0.0958** |
| WRF `FMC_GC_F`, mean over the burn window | **0.1205** |

A gap of ~0.025, which from the phase 3 curve is worth roughly **1.4x in predicted
area**. For contrast, on Ranger Road the same two sources agreed to **2%**
(0.063 vs 0.065). Agreement is fire-dependent and much worse in the
timber/complex-terrain case, consistent with the known difficulty of estimating
fuel moisture where observations are sparse.

What a moisture error costs, interpolated from the measured curve:

| centre | +/-0.01 | +/-0.02 | +/-0.03 |
|---|---|---|---|
| Md 0.06 | 1.3x | 1.9x | **2.6x** |
| Md 0.08 | 1.3x | 1.6x | 1.9x |
| Md 0.10 | 1.1x | 1.3x | 1.4x |

The penalty is worst exactly where fire-weather conditions sit. **This sets a
floor on achievable forecast accuracy that no amount of parameter tuning can go
below**, and worse, a fitter will silently absorb moisture error into wRF and
return a confidently wrong wind parameter. Worth quantifying this floor per fire
before deciding whether a wRF fit on that fire is meaningful at all.

### Change made: SAV-weighted dead moisture

`read_fmda_fmc` now takes `classes` and `weights` instead of a single
`class_index`; a single class reproduces the old path exactly. `resolve_md`
builds them from new config keys:

```json
"weighting": "sav",              // or "single", which uses class_index
"dead_classes": [0, 1, 2],       // 1h, 10h, 100h
"sav_weights": [2000.0, 109.0, 30.0]
```

Rationale: Rothermel's single dead moisture denotes a surface-area weighted
quantity across size classes, not the 1h class alone.

### But it changes almost nothing, and an earlier claim here was wrong

| fire | 1h only | SAV-weighted | difference |
|---|---|---|---|
| SINLAHEKIN | 0.0958 | 0.0963 | **+0.0005** |
| Ranger Road | 0.0634 | 0.0666 | +0.0032 |

0.5% on the timber fire. **SAV weighting is itself ~94% 1h** — with ratios
2000 : 109 : 30 the fine fuels dominate by construction — so the weighted mean is
essentially the 1h value.

**Correction.** Earlier in this session it was claimed that in timber "the fire
is carried substantially by 10h and 100h fuels", so that a single `Md` misrepresents
it. That is **wrong for spread rate**. Rothermel's spread is fine-fuel dominated
by design; larger size classes govern intensity and residence time, not rate. The
single-`Md` formulation using something close to 1h is well justified even in a
timber domain.

The change is still worth keeping — it is the quantity `Md` actually denotes, it
costs nothing, and it degrades gracefully if 1h is missing or bad. It is **not**
a fix for timber-domain moisture error. That error is observational, per the
0.025 source disagreement above.

Verified: shipped default (moisture disabled) still inert; `class_index` still
honoured under `"single"`; a class outside the tile rejects with a specific
message; the provenance string reports which weighting was used.

---

## 5e. Numerical resolution: NOT converged, and it matters as much as physics

Eight runs on SINLAHEKIN, one `sweep_ff_params` call each so wall clock is
per-set. Driver `$SP/resolution_study.py`.

### perimeterResolution (spatialIncrement held at 5)

| value | runtime | speedup | area | d area | perimeter | d perim |
|---|---|---|---|---|---|---|
| 15 | 10.4 min | 1.00x | 3,954 ha | +0.0% | 48.4 km | +0.0% |
| **25** (template) | 4.5 | 2.30x | **3,515** | -11.1% | 39.7 | -18.0% |
| 40 (ForeFire default) | 2.8 | 3.64x | 3,000 | -24.1% | 30.3 | -37.4% |
| 60 | 2.5 | 4.15x | 2,644 | -33.1% | 25.4 | -47.5% |
| 100 | 2.4 | 4.38x | 2,426 | -38.7% | 21.9 | -54.7% |

### spatialIncrement (perimeterResolution held at 25)

| value | runtime | speedup | area | d area | perimeter |
|---|---|---|---|---|---|
| 2 (ForeFire default) | 7.3 min | 1.00x | 3,191 ha | +0.0% | 36.7 km |
| 5 (template) | 4.5 | 1.62x | 3,515 | +10.2% | 39.7 |
| 10 | 3.6 | 2.01x | 3,733 | +17.0% | 46.3 |
| 20 | 3.0 | 2.43x | 3,960 | +24.1% | 47.0 |

The `sInc = 5` row is the **`pRes_25` run reused** -- template settings are
`pRes 25` + `sInc 5`, so that point is shared between the two families and only
ran once. It therefore does not appear in the study's own `spatialIncrement`
summary table in `resolution.log`; the speedup and delta columns for it are
derived against the `sInc = 2` baseline.

### What this means

**There is no plateau in either parameter.** Area follows `pRes^-0.26` and
perimeter `pRes^-0.42` across the whole range with no sign of asymptoting, so the
solution is **not numerically converged at any setting tested, including finer
than the template**. Perimeter moves about twice as much as area, which is why a
20 m vs 40 m comparison judged on area alone reads as "no effect" -- the two
differ by ~12-15%, easy to mistake for noise without a baseline.

The two parameters push in **opposite directions**: coarsening `pRes` shrinks the
fire, coarsening `sInc` grows it. They partially cancel, which further disguises
the lack of convergence.

**This is as strong a lever as the physics.** `pRes` 25 -> 40 changes area by 15%,
the same as `wRF` 0.4 -> 0.5 on this fire. Consequences:

- A `wRF` calibrated at one `pRes` is not valid at another. **Fix `pRes` before
  fitting anything physical, and record it alongside any fitted parameter.**
- The 2.3x speedup from `pRes` 25 -> 40 is real but is not free: it costs 15% of
  the area. For low-latency previews that may be an acceptable trade, but it must
  be a stated choice, not an implicit one.

Caveat on the measurement: the `initial_area_ha` column is the **first 30-minute
perimeter**, not the ignition triangle, since that is the earliest output
ForeFire writes. It does not isolate the seed-size confound as intended. It does
show the confound is negligible by the first output step (130.8 vs 128.4 ha at
`pRes` 40 vs 60).

---

## 5f. ForeFire vs WRF-SFIRE vs observation -- ForeFire over-predicts ~3x

### The observation

A ground estimate of **453 acres (183 ha) at 2026-07-27 05:50 UTC**, which is
**+8.64 h** after the 21:11:19 ignition. Both model perimeters reportedly fall
inside the GOES detection footprint at that time, which is weak confirmation
given GOES resolution but not nothing.

| source | area at +8.64 h | vs observed |
|---|---|---|
| **observed** | **183 ha** | -- |
| WRF-SFIRE (`TIGN_G`) | 165 ha | **0.90x** |
| ForeFire, template defaults | 599 ha | **3.27x** |

**WRF-SFIRE is within 10% of the observation; ForeFire is 3.3x too large.** An
earlier caution in this session -- that WRF-SFIRE might itself be
under-predicting -- is answered: it is not.

### The full-run discrepancy grows with time

WRF-SFIRE burned area is recovered from `TIGN_G` on the fire subgrid (cells with
an ignition time below the 108,006 s sentinel; 9,051 of 1,440,000 cells at the
end). Verified against `FGRNHFX > 1 kW/m2`.

| hours since ignition | WRF-SFIRE | ForeFire | ratio |
|---|---|---|---|
| 4 | 76 ha | 169 ha | 2.2x |
| 12 | 240 | 915 | 3.8x |
| 20 | 298 | 1,695 | 5.7x |
| 30 | **566 ha** | **3,515 ha** | **6.2x** |

A compounding rate difference, not a fixed offset.

### perimeterResolution cannot be used to tune this away

Tempting, since WRF-SFIRE output is available where observations are not. But
the whole `pRes` range spans only **6.99x down to 4.29x** relative to WRF-SFIRE.
The coarsest setting tested still leaves ForeFire **4.3x too large** -- the knob
moves 39% against a 620% gap. `sInc` makes it worse, not better. **No combination
of the numerical parameters closes it.**

Only **Md ~ 0.18** reproduces WRF-SFIRE's 566 ha, but FMDA says 0.096 and WRF's
own `FMC_GC_F` says 0.121 for this fire -- so matching demands a moisture ~50%
above the wettest defensible estimate. Correcting ForeFire to WRF's actual 0.121
only reaches ~2,700 ha, still 4.8x. **Roughly 20% of the gap is moisture and 80%
is model formulation** (different Rothermel variant, coupled vs static winds,
fire-atmosphere feedback).

The hazard in tuning `pRes` to match: it compensates a physics discrepancy with a
discretisation error, and the two scale differently with fire size, duration and
terrain, so a `pRes` fitted on one fire will not transfer. The defensible version
is to **fix `pRes` by convergence, not by agreement**, and treat the residual as a
model discrepancy to explain.

### And wRF cannot reach the observation at all

At +8.64 h, `wRF` 0.2 -- the bottom of the swept range -- still gives 365 ha,
**2x the observation**. The closest settings are `pSAF` 0.3 (192 ha, 1.05x) and
`Md` 0.15 (220 ha, 1.20x). This is the leverage result of section 5 meeting real
data: on this fire `wRF` cannot span the distance to truth even at its extreme,
so **fitting `wRF` against observed area here is not merely imprecise, it is
infeasible**.

---

## 5g. Joint pSAF x Md grid against the observation -- the calibration works

25 runs, pSAF {0.2 .. 0.6} x Md {0.08 .. 0.18}, truncated to a **10.5 h horizon**
(16 steps) so the +8.64 h observation is bracketed with margin. 375 steps, no
errors, **13.3 min for the whole grid** -- roughly a tenth the cost of the 30 h
runs. Driver `$SP/grid_psaf_md.py`, figure
`<wksp>/forefire/grid_pSAF_Md_vs_observation.png`.

### Area (ha) at +8.64 h; observation = 183 ha

| | Md 0.08 | Md 0.10 | Md 0.12 | Md 0.15 | Md 0.18 |
|---|---|---|---|---|---|
| pSAF 0.2 | 122 | 111 | 91 | 37 | 17 |
| pSAF 0.3 | 223 | **192** | 159 | 64 | 29 |
| pSAF 0.4 | 345 | 299 | 251 | 116 | 41 |
| pSAF 0.5 | 513 | 428 | 345 | *65* | 57 |
| pSAF 0.6 | 679 | 599 | 471 | 219 | 87 |

### The Md = 0.15 column is unusable -- a knife edge, not noise

Fuels **2 and 11 both have `me` = 0.15**, so at Md = 0.15 exactly,
`Mratio = 1.0` and `Etam = 0.000000`. Fuel 2 is 16.9% of the domain sitting
precisely on its extinction point, and floating-point differences flip cells
between burning and not. That column is the only non-monotonic one (pSAF 0.5
gives 65 ha, *less* than pSAF 0.4's 116, which is impossible for a pure speed
multiplier). The other four columns are strictly monotonic.

**Rule: never sweep Md through a value equal to any fuel's `me`.** Offset by at
least 0.005.

### The ridge, from the trustworthy columns

| pSAF | Md that reproduces 183 ha |
|---|---|
| 0.2 | **no match** -- maxes at 122 ha even at Md 0.08 |
| 0.3 | 0.105 |
| 0.4 | 0.139 |
| 0.5 | 0.154 |
| 0.6 (template) | 0.165 |

Degenerate as predicted -- one area at one time cannot separate two isotropic
scalings. But the ridge is **sloped, not flat**, so it rules out pSAF <= 0.2
outright and collapses to a point once moisture is fixed independently.

### With Md pinned at the FMDA estimate (0.096)

| pSAF | area at +8.64 h | vs observed |
|---|---|---|
| 0.2 | 113 ha | 0.62x |
| **0.3** | **198 ha** | **1.08x** |
| 0.4 | 307 | 1.68x |
| 0.6 (template) | 614 | 3.35x |

**pSAF ~ 0.3, half the template's 0.6.** Self-consistent: the ridge says pSAF 0.3
wants Md 0.105 while FMDA independently says 0.096 -- agreeing to within a point
of moisture, inside the FMDA-vs-WRF disagreement of section 5d.

This is the workflow that works: **fix moisture from FMDA, fit pSAF.** Note it is
*not* the wRF fit originally planned -- section 5f showed wRF cannot reach the
observation on this fire at any value in range.

---

## 6. Tooling, and where it lives

Three scripts were written this session, all in the session **scratchpad**, which
**will not survive**:

- `sinlahekin_sweep.py` — the three-phase driver; phase 2's probes derived from
  phase 1's spread rates.
- `sweep_analyze.py` — generic over fire: walks `<wksp>/forefire/<tag>/`, reads
  the chained step perimeters, writes `sweep_curves.csv` and `sweep_sets.csv`
  (area, head radius, mean ROS, and burned-area/head-radius-circle fill, in
  EPSG:5070). Knows both fires.
- `sweep_plots.py` — four figures into `<wksp>/forefire/`:
  `growth_curves_windReductionFactor.png`, `growth_curves_wRF_vs_pSAF.png`,
  `fuel_moisture_sensitivity.png`, and the new `cross_fire_sensitivity.png`
  (both fires, each normalised to its own default run).

**These are worth promoting into the repo** if the sweeps are going to be a
recurring workflow — `sweep_analyze.py` in particular is fire-agnostic and is the
piece that would otherwise be rewritten each time. Not done without instruction.

Colour palettes are the ordinal ramps validated on 2026-09-08 (blue = wRF,
orange = pSAF, aqua = Md), re-validated against the current tooling build this
session. The cross-fire figure encodes fire identity rather than magnitude, so it
uses categorical slots 1 and 2.

---

## 7. Open items carried forward

Everything in §10 of the 2026-09-08 handoff still stands. Added today:

1. The `forefire.py` mtime discrepancy in §1.
2. The analysis scripts are ephemeral (§6).
3. Sweep results not yet analysed or plotted — §5's prediction is untested.
4. `phiV`/`phiP` leverage is a point estimate; a domain-integrated version would
   be more defensible if this is used to justify calibration decisions.

---

## 8. The operational goal, and what it changes

Stated this session, and it should steer future work:

> Low-latency forecasts to about **8 hours**, for potentially **many hotspots
> appearing at once**, to support early-phase response. Long forecasts are not
> the objective; the critical decisions happen early.

This reframes several of today's results.

**The 30 h runs were the wrong unit.** A 10.5 h truncated run costs ~50 s per
parameter set against 4.5 min for 30 h. The whole 25-cell grid took 13.3 min. For
the many-hotspots case that throughput is the number that matters, and it is
already adequate.

**A short horizon is kinder to the numerics.** Resolution sensitivity roughly
halves at 8 h versus 30 h -- pRes spread 1.31x vs 1.63x, exponent -0.154 vs
-0.270 -- because discretisation error compounds with time. Caveat: the flattening
between pRes 60 and 100 at 8 h is **not** convergence; the answer is still moving
at the fine end (15 -> 25 is -7%). The short horizon halves the problem rather
than removing it, so `pRes` must still be fixed and recorded.

**The +8.64 h observation is the primary validation target, not a spot check.**
It sits almost exactly at the operational horizon, which makes the 5g calibration
directly relevant rather than incidental.

**Multi-point ignition (5c) is the natural fit**, not an extra. "Many hotspots at
once" is exactly repeated `startFire[lonlat=...]`, and re-igniting from each new
detection absorbs observed spotting without needing the physical brand model that
5b showed is unreachable standalone.

### Suggested shape for the next sessions

1. Run several more fires through the same short-horizon workflow -- the two done
   so far (plains grass, complex-terrain timber) already behave very differently,
   so the sample needs widening before any parameter choice is trusted.
2. Use the **10.5 h truncated grid** as the standard unit, not 30 h runs.
3. Per fire, record: fuel mix, terrain relief, wind coherence, the `phiV`/`phiP`
   split, and which knob is identifiable. On Ranger Road that was wRF; on
   SINLAHEKIN it is pSAF. **Do not assume it transfers.**
4. Compute the fuel-model ROS ceiling (5f) up front -- it is cheap and tells you
   whether the observed fire is even reachable before any fitting.
