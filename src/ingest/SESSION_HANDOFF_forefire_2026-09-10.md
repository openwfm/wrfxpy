# ForeFire session handoff — 2026-09-10

Branch `james_ngfs`, HEAD `3e5d4a2` (yesterday's commit of `forefire.py`, the
config example and the two prior handoffs). Read
`SESSION_HANDOFF_forefire_2026-09-09.md` first for the sweep results this builds on.

Today is about **shape**, not size: comparing modelled perimeters against a real
IR-interpreted perimeter for SINLAHEKIN.

**Headline:** a ~500 m ignition-point error was costing ~0.26 IoU. Corrected, the
model reaches **IoU 0.70 / Dice 0.82** on a 6-hour forecast. Origin accuracy is
worth more than any parameter in the study.

---

## 1. `forefire.py` edits since the commit, now diffable

The commit gave the module a baseline, and the first diff shows what changed
overnight (4 insertions, 8 deletions):

| change | effect |
|---|---|
| `if len(combined_gdf) > 0:` before `to_file` | guards KML writing on an empty frame |
| **`for gg in g:` -> `for gg in wksp_dirs:`** | **fixes the `run_days` loop bug** |
| `run_days(days2run=40)` -> `7` | a week's window |

The loop fix matters now: with `days2run=7` the loop must actually iterate all
days, where at `days2run=2` running only the last day was tolerable. This also
closes the unexplained mtime noted in 2026-09-09 §1 — the committed file had
`days2run=40`, so that edit was almost certainly the same hand.

---

## 2. The observed perimeter

`/data/jhaley/wrfxpy/ngfs/perims/SINLAHEKIN_{78D35D3B-F791-4961-AE36-C6D1A4DFF5A0}.geojson`

One polygon, `poly_MapMethod = IR Image Interpretation`, `poly_GISAcres` 453.30.

| property | value |
|---|---|
| area | **183.4 ha** (453.3 acres, matches metadata exactly) |
| perimeter | 6.74 km |
| min-rotated-rect | 2617 x 1032 m, **elongation 2.54** |
| long-axis bearing | 175 deg from north |
| centroid | 48.70781, -119.70964 |

Useful context in the attributes: `attr_PrimaryFuelModel` "Timber (Litter and
Understory)", `attr_SecondaryFuelModel` "Timber (Grass and Understory)", and
`attr_FireBehaviorGeneral` **"Extreme / Uphill Runs / Short-range Spotting"** —
independent confirmation of the slope-dominated spread the `phiP` analysis
predicted (09-09 §5) and of spotting, which 09-09 §5b showed ForeFire cannot
produce standalone.

### Timestamps disagree, and I used the wrong one

| stamp | value | +h after ignition | what it is |
|---|---|---|---|
| `poly_PolygonDateTime` | 2026/07/27 03:05:42 | **+5.91** | the IR flight |
| `poly_CreateDate` | 05:47:09 | +8.60 | record created |
| `poly_DateCurrent` | 05:49:07 | +8.63 | record updated |
| `processed_utc` | 14:00:21 | +16.8 | **cache download**, `perim_cache.py` |

`processed_utc` is `datetime.now(timezone.utc)` in
`update_incident_store_geojson` — reliable, but a download time, so only an
**upper bound** on when the perimeter was real. Source dates in these files are
known to be unreliable generally.

**Correction to 09-09 §5f.** That section used +8.64 h, the record-update stamp.
The IR flight is +5.91 h. At the correct time:

| | at +5.91 h (correct) | at +8.64 h (as recorded 09-09) |
|---|---|---|
| observed | 183 ha | — |
| WRF-SFIRE | 112 ha — **0.61x** | 165 ha — 0.90x |
| ForeFire, template | 320 ha — **1.74x** | 599 ha — 3.27x |

So **WRF-SFIRE under-predicts by 39%**; it does not match within 10%. ForeFire
over-predicts by 1.7x, not 3.3x. Both models are wrong in *opposite* directions
with the truth between them — a materially different picture from 09-09 §5f, and
that section should be read with this correction.

The calibration shifts too: with Md pinned at FMDA's 0.096, the best pSAF moves
from ~0.3 to **~0.4** (174 ha, 0.95x), much closer to the template's 0.6.

---

## 3. Shape comparison method

Compared at **matched area**, not matched time — the timestamps disagree by
hours, so "when the model is this big, is it this shape?" is the cleaner
question and removes timing from the shape result entirely.

Metrics in EPSG:5070. **IoU** (Jaccard) is
`area(model ∩ obs) / area(model ∪ obs)`; 1.0 identical, 0.0 disjoint. It
penalises both missed fire and invented fire, so `covered`
(fraction of the observation captured) and `spill` (fraction of the model
outside it) are reported alongside. The Sorensen-Dice index is monotonically
related: `Dice = 2·IoU/(1+IoU)`, always higher.

**Scoring expectation, per the user: for a fully automated system IoU ~ 0.3 is
already good.** An earlier claim in this session that "0.6 is a decent match"
was asserted without support and is not a standard to judge against.

Script: `$SP/shape_compare.py` (scratch).

---

## 4. Shape at the original ignition — the origin dominates

Across 20 runs spanning the full wRF, pSAF and Md ranges, at matched area:

- IoU **0.32–0.46**, tightly clustered — best was `wRF 0.2` at 0.46
- modelled elongation **1.23–1.76** against the observed **2.54**
- bearing within 9–15 deg of observed for the good runs, 1 deg for `wRF 0.8`

**Shape did not break the pSAF/Md degeneracy.** IoU varies by only 0.14 across
every combination; ridge points that fit area equally well score 0.42–0.44 and
are indistinguishable. Both are isotropic scalings, so at matched area they
produce the same shape — as the `R0`-multiplier argument predicted.

Only wRF moves shape (elongation 1.23 at 0.2 to 1.76 at 0.8), being the one
anisotropic knob, but it cannot reach 2.54: the fit is
`elongation ≈ 0.63·wRF + 1.11`, needing **wRF ≈ 2.3, 5.7x the template** and far
past the 1.0 discontinuity where the Andrews cap formula switches. Noisy fit
(0.7 dips below 0.6), so treat as order-of-magnitude.

**12 of 25 grid combinations never reach 183 ha at all** — the whole pSAF 0.2
row and everything at Md 0.18 — so the observation rules out a real part of the
parameter space regardless of timing.

---

## 5. The winds are right; the ignition point is wrong

WRF surface winds over the fire footprint (`UF`/`VF`, what ForeFire reads) blow
**northward throughout**: toward 337–356 deg from ignition to 02Z, then 14–49 deg,
with `v` positive until 08Z. Saveouts are reduced files with no `XLAT` or `U10`,
so `UF`/`VF` is the available and relevant field.

So the model's northward spread is faithful to its wind field. **Wind direction
is not the error.**

Translating the model perimeter to best overlap the observation:

| run | IoU as-is | IoU shifted | shift E | shift N | distance |
|---|---|---|---|---|---|
| wRF 0.4 | 0.41 | **0.63** | -200 m | -500 m | **539 m** |
| wRF 0.8 | 0.32 | 0.61 | -100 | -700 | 707 |
| grid pSAF0.5/Md0.1 | 0.44 | **0.65** | -100 | -500 | 510 |

A ~500 m shift south-southwest lifts IoU from 0.41 to 0.63 — **a bigger single
improvement than any parameter change achieved**. All three runs agree on
direction and magnitude, consistent with the shape being right and only the
origin wrong.

- `input.json` ignition: 48.70860, -119.70411
- implied origin: **48.70375, -119.70506** (539 m away)
- southernmost point of the observed perimeter: 48.69762, -119.70533
- `attr_InitialLatitude/Longitude`: 48.6922, -119.7137 — **1.96 km away** and
  the worst of the candidates

**Reported ignition points in these files are unreliable** and often sit outside
the perimeter, frequently on roads where a responder observed the fire. There is
no standard method behind them. The operational system waits for VIIRS
detections to place the origin better, and this quantifies why: **a 500 m origin
error costs ~0.22 IoU, larger than the spread across every parameter combination
tested.** Origin accuracy buys more than parameter tuning.

### Correction: my "the fire ran south" reading was wrong

From the overlay I claimed the observed fire ran south, opposite the wind. It did
not. The observed perimeter *extends* further south because it *started* further
south. The perimeter shows where the fire is, not which way it went, and the
axis metric I used is direction-blind (it treats 175 and 355 deg as one line).
The fire ran north, as WRF and the model both say.

---

## 6. Grid re-run from the corrected ignition — IoU 0.70

`$SP/grid_ign2.py`, 25 sets, pSAF {0.2..0.6} x Md {0.08, 0.10, 0.12, 0.14, 0.18},
10.5 h horizon, tags `ign2_*`. **Md 0.15 replaced by 0.14** — 0.15 exactly equals
the extinction moisture of fuels 2 and 11, making `Etam` identically zero and
that column numerically unstable (09-09 §5g).

`input.json` is **not** modified; `read_input` is overridden in memory only. The
monitor watches the file's checksum and aborts on any change.

**Circularity caveat, important:** the corrected ignition was *derived from the
observation being fitted*. A good IoU here is **not** independent validation of
the origin. It tests only whether, given a correct origin, the parameters can
match shape as well as size. Real validation needs an origin from an independent
source — early VIIRS detections being the intended one.

### Results — 25 sets, 12.4 min, `input.json` verified unchanged

IoU at matched area:

| | Md 0.08 | Md 0.10 | Md 0.12 | Md 0.14 | Md 0.18 |
|---|---|---|---|---|---|
| pSAF 0.2 | — | — | — | — | — |
| pSAF 0.3 | 0.69 | 0.69 | 0.65 | — | — |
| pSAF 0.4 | **0.70** | **0.70** | 0.68 | 0.62 | — |
| pSAF 0.5 | **0.70** | 0.69 | 0.68 | 0.66 | — |
| pSAF 0.6 | 0.67 | 0.68 | 0.67 | 0.63 | — |

Best six:

| pSAF | Md | h | ha | IoU | Dice | cover | spill | elong | centroid |
|---|---|---|---|---|---|---|---|---|---|
| 0.4 | 0.08 | 6.5 | 191 | 0.70 | 0.82 | 0.84 | 0.19 | 1.61 | 78 m |
| 0.5 | 0.08 | 5.0 | 175 | 0.70 | 0.82 | 0.80 | 0.16 | 1.57 | 84 m |
| 0.4 | 0.10 | 7.0 | 179 | 0.70 | 0.82 | 0.81 | 0.17 | 1.60 | 45 m |
| 0.3 | 0.10 | 10.0 | 188 | 0.69 | 0.82 | 0.83 | 0.19 | 1.42 | 75 m |
| 0.3 | 0.08 | 8.5 | 184 | 0.69 | 0.82 | 0.82 | 0.19 | 1.46 | 48 m |
| 0.5 | 0.10 | 5.5 | 175 | 0.69 | 0.81 | 0.80 | 0.17 | 1.52 | 65 m |

**IoU 0.62–0.70, against 0.32–0.46 from the original ignition.** Better even than
the 0.63–0.65 obtained by rigidly translating the old perimeters, because the
model re-ran from the new origin and grew into the right terrain rather than
being slid across it. Against the "IoU ~0.3 is good for a fully automated
system" bar, **0.70 (Dice 0.82) is a strong 6-hour forecast**: it captures 84% of
the observed area with only 19% of its own footprint outside.

Three things follow.

**Origin dominates everything.** Fixing a 500 m placement error moved IoU by
~0.26. No parameter moves it by more than 0.08. For the low-latency multi-hotspot
workflow this is the single highest-value input, and it is exactly what
VIIRS-derived ignition points would supply.

**IoU still cannot pick parameters.** The grid spans only 0.62–0.70 and the top
six are within 0.01 of each other. Shape does not break the pSAF/Md degeneracy
even with the origin corrected — consistent with both being isotropic scalings.

**The elongation gap is structural, not a placement artifact.** Modelled
elongation is 1.42–1.80 against the observed 2.54, essentially unchanged from
the 1.23–1.76 measured before the origin fix. This answers the open question:
the shape deficit survives correct placement, so it is not an artifact of
comparing misaligned footprints. Short-range spotting and uphill runs — both
reported for this fire, neither available standalone — remain the likely cause.

**Timing offers a weak extra constraint.** Matched area occurs at different
model times: pSAF 0.5 reaches 183 ha at 5.0–5.5 h, pSAF 0.3 not until 8.5–10 h.
If `poly_PolygonDateTime` (+5.91 h) is trusted, that favours **pSAF ~0.5**. If
the timestamp is not trusted — and source dates in these files often are not —
then area and shape together still cannot separate the ridge.

Also: **10 of 25 combinations never reach 183 ha** — the entire pSAF 0.2 row and
all of Md 0.18 — so the observation continues to exclude a real part of the
parameter space independent of any timing assumption.

---

## 7. The timing bug, and the fix

**The intent** (stated by you, and it is the right convention): ForeFire's `t=0`
is the timestamp of the latest wrfout that exists *before* the fire's ignition
time. WRF-SFIRE instead sets its `t=0` to the *first* wrfout of the run, which
can be several hours earlier.

`make_timing_table` was assigning `ig_s = ws`, i.e. seconds from the first
wrfout. For SMOKEHOUSE CREEK (WRF from 18:00Z, ignition 20:41:17Z) that handed
every run a **2 h 41 min head start**; for SINLAHEKIN it was about 3 h.

The fix keys off the step that brackets the ignition and re-references
everything to it:

```python
ws_ignition = 0.0
...
elif (ws < ignition_seconds) and (ws + t_step > ignition_seconds):
    ws_ignition = ws          #this wrfout is ForeFire's t=0
    ig_s = ignition_seconds - ws
else:
    ig_s = ws - ws_ignition   # was: ig_s = ws
```

Verified on SMOKEHOUSE: row 5 is the 20:30Z wrfout and gets
`ign_seconds = 677`, which is exactly 20:41:17 − 20:30:00. Row 6 (21:00Z) gets
1800. Both correct.

**Before this, I twice claimed a wrong magnitude for this bug** — first "2×",
then "corrected" it to 1.25× using a `valid_at` field I had misread. Neither was
right. The real defect only became findable once you stated the intended
convention; I had been inferring it from the code, which is the thing that was
wrong.

### What this invalidates

Sensitivities (elasticities, knee positions, the pRes scaling law) are ratios
between runs that all shared the same offset, so they stand. **Absolute
model-vs-observation numbers do not.** Specifically:

| Result | Status |
|---|---|
| 09-09 §5f, 6.2× vs WRF-SFIRE | invalid, ran with a ~3 h head start |
| 09-10 §2, 1.74×/1.83× oversize | invalid, same cause |
| the pSAF ≈ 0.4 calibration | invalid, rests on the above |
| all wRF / pSAF / Md elasticities | unaffected |
| IoU and the shape/origin work (§3–6) | unaffected — matched-area comparison |

The IoU work is safe precisely *because* it compares at matched area rather than
matched time, which was luck rather than foresight.

## 8. Multi-point ignition

Three new functions in `forefire.py`, plus a `multi_ignition` config block
(`enabled`, `min_separation_m`, `max_per_step`), off by default:

- `read_ignitions(wksp_dir)` — every point in `input.json`, not just
  `ignitions['1'][0]`, sorted by time.
- `thin_ignitions(points, min_separation_m)` — grid-hashed greedy thinning.
  **Not optional**: ForeFire builds each seed as a triangle of
  `2*perimeterResolution`, so co-located seeds make degenerate fronts, which it
  trashes and then segfaults on.
- `ignition_schedule(points, timing_table, max_per_step)` — assigns each point to
  the step whose `[ign_seconds[i], ign_seconds[i+1])` window contains it. This is
  the shape you described: each restart waits on the next wrfout, and detections
  arriving during the wait get added to that step.
- `apply_ignitions(text, entries)` — strips the template's single
  `startFire[lonlat=...]` and inserts the batch before the first `goTo[`.

`sweep_ff_params(..., ignitions=points)` threads it through; `run_forecasts`
picks it up automatically when `multi_ignition.enabled` is true.

**Bug found and fixed during testing:** `thin_ignitions` first scaled each
point's longitude by *its own* `cos(lat)`, which is not a projection — over a
35 km latitude span the same longitude lands tens of km apart and the distance
test silently fails. Measured 492 m separation against an 800 m target. Using one
reference latitude for the whole cloud makes it exact (400/800/1600 m all hit).

### Capacity and front collisions

Measured on SMOKEHOUSE, short runs:

- 500 **distinct** points run fine; 2000 distinct segfault. 500 *duplicated*
  points crash, which is the degenerate-triangle path, not a count limit.
- Fronts merge correctly. Over one step with **no new ignitions**, polygon count
  fell 1552 → 1427 while area grew ~10,000 ha — fronts coalescing, not vanishing.
- Same effect over 8 h at lower seed counts: 42 seeds → 4 polygons.
- 6 steps with 400 ignitions per script completed in 2.6 min.

So collisions are handled; the practical limit is a few hundred *distinct*
concurrent fronts, and thinning is what keeps you under it.

## 9. SMOKEHOUSE CREEK's `input.json` is a hindcast, not a forecast input

The workspace carries **17,645 ignition points**, all at a single time
(2024-02-27 20:41:17Z). You noted these were artificially augmented from
detections toward the fire centre. The numbers confirm that seeding from them
cannot test forecast skill:

| quantity | value |
|---|---|
| convex hull of the 17,645 seeds | 366,396 ha |
| WRF-SFIRE total burn, same workspace | 353,989 ha |
| WRF-SFIRE burn within the **first hour** | 309,493 ha |
| GOES detected footprint, full event | 477,449 ha |

The seed cloud is already the size of the final fire. WRF-SFIRE burns 87% of its
total in hour one because it is *told* where the fire is, not because it spreads
there. Thinning does not rescue this: 800 m thinning drops the seeded area from
43,760 ha to 5,674 ha, discarding 87% of it, because the points are far finer
than the 100 m fire grid.

Note also that the real Smokehouse Creek fire ignited on **2024-02-26**, a day
before this workspace's nominal ignition. The 20:41:17Z timestamp is simply when
NGFS first associated detections with the named incident.

## 10. A real GOES-driven forecast

Source: `ingest/NGFS/NGFS_FIRE_DETECTIONS_GOES-16_ABI_CONUS_2024_02_27_058.csv`
(69,979 rows, 52 columns, 274 scans at the 5-minute ABI CONUS cadence).

### Selection rule

Rows with `known_incident_id = {4A55159B-D06F-4574-A689-CC4CDCCDC097}`: **6,437
detections**, 20:41:17Z–23:56:17Z, 35 scans. Deduplicated at 2 km (one ABI pixel
is 7.4 km², about 2.7 km across, and detections sit on a ~2.4 km grid) → **558
distinct pixels**. Note 800 m thinning is a **no-op** on real GOES data: all 558
survive it. The thinning code matters for augmented point clouds, not for raw
GOES.

### `feature_tracking_id` back-tracing: tried, and rejected

You noted the tracking id is assigned to a location *before* that location is
associated with a known incident, so back-tracing should recover pre-association
detections. It does — but it is contaminated:

- The ids are timestamped strings (`ID-2024-02-27T20:36:30Z_0029`), globally
  unique, **not** recycled integers. 2,156 distinct in the file.
- The 12 ids on Smokehouse rows appear in 15,179 rows spanning the whole day.
- Within 10 km of the origin and before 20:00Z, that pull returns **30 pixels
  present on 218 of the day's 274 scans at steady 200–350 MW FRP**. Those are
  persistent heat sources — gas flares, which the Texas panhandle has many of —
  not fire.
- One id (`…02-26T21:21:30Z_0024`, 2,714 rows) has a mean position 57 km away.

**A persistence filter is required before tracking-id back-tracing is usable in
the pipeline**: reject any pixel detected on more than some fraction of the
preceding scans. This is worth building — it is exactly the failure that would
otherwise seed phantom fires on flares.

### The first scan is three fires, not one

Clustered at 6 km linkage, the 42 detections at 20:41:17Z resolve to:

| cluster | px | centre | span | ΣFRP |
|---|---|---|---|---|
| 0 | 19 | 35.807, −101.207 | 21.4 km | 8,142 MW |
| 1 | 18 | 35.792, −101.466 | 12.4 km | 11,363 MW |
| 2 | 5 | 35.827, −100.921 | 6.9 km | 1,116 MW |

The workspace's single ignition point is cluster 2 — the **smallest** of the
three, at 1,116 MW, with the two larger fires 25 km and 50 km further west. Any
single-point run on this workspace is seeded on the least significant part of the
complex.

### Detection arrival is front-loaded

42 new pixels at 20:41, then 151 more at 21:11, then roughly 10 per scan. The
21:11 jump is NGFS re-associating a batch, not fire growth — it shows up in the
observed area curve as a 35,968 → 162,408 ha step in 30 minutes, so growth rates
must be fitted **after** it or they are meaningless.

### Verification data — use the incident cache, not the daily CSV

`ngfs/incident_data/SMOKEHOUSE_CREEK_{4A55159B-...}.pkl` is much better for
scoring than the daily CSV:

- **12,004 detections running to 2024-02-28 13:12Z**, where the single-day CSV
  stops at 23:56Z after only 3.25 h. It covers the whole 8 h forecast window.
- It carries **true pixel corner coordinates** (`lat_tc_c1..c4`,
  `lon_tc_c1..c4`), so the observed footprint is the union of exact ABI
  quadrilaterals rather than squares of equal area.
- All GOES-16 CONUS; there is **no VIIRS** in this cache, so there is no finer
  observation available for this fire from here.

The observed footprint remains an **upper bound** on burned area: a 7.4 km² pixel
flags if any part of it is hot. Its 477,449 ha against WRF-SFIRE's 353,989 ha
suggests roughly 35% over-coverage, which is the right order for whole-pixel
unions.

### Result: cold start from the first scan

42 seeds at 20:41:17Z, no further information, 8.3 h forward on the corrected
clock, template parameters (wRF 0.4, pSAF 0.6, pRes 100 m, spatialIncrement 5).
Ran in **1.2 minutes** for 17 chained steps — comfortably inside the low-latency
budget.

| h | model ha | observed ha | ratio | IoU |
|---|---|---|---|---|
| 0.3 | 456 | 35,968 | 0.013 | 0.013 |
| 0.8 | 624 | 162,408 | 0.004 | 0.004 |
| 2.3 | 1,224 | 304,089 | 0.004 | 0.004 |
| 4.3 | 2,094 | 447,873 | 0.005 | 0.005 |
| 8.3 | 5,383 | 476,068 | 0.011 | 0.011 |

Over the clean window (0.8–8.3 h, after the re-association step):

- model **634 ha/h**, observed **41,821 ha/h** — a **66× shortfall** in area
  growth rate.
- equivalent-radius rate: model 0.107 m/s, observed 0.547 m/s — **5.1×** in
  linear terms.
- at 8.3 h the model has burned **1.5%** of WRF-SFIRE's total for the same fire.

**Do not read this as a Rothermel result.** §11 shows most of it is an input
defect: 86% of these seeds were placed on non-burnable fuel. The numbers above
are what the workspace produces, not what the model is capable of.

The separate claim that Rothermel cannot reach 10 m/s still stands, but on the
independent evidence of 09-09 §5b — a measured ceiling of **2.66 m/s in fuel 2**,
because `phiV` saturates — not on this run.

## 11. The SMOKEHOUSE workspace's fuel map has the fire scar in it

This is the important finding of the session, and it invalidates §11 as a growth
test.

`wrfinput_d01`'s `NFUEL_CAT` puts **5.2% of the domain in fuel category 14**. In
`fuelstrans.csv`, category 14 has **`e = 0`** — zero fuel bed depth — so
Rothermel returns no spread there. It is effectively "no fuel", and WRF-SFIRE
sets `UF`/`VF` to exactly 0.00 m/s in those cells, which is how I first noticed
it.

Those cells are not scattered roads or fields. They are one contiguous block:

| quantity | value |
|---|---|
| largest fuel-14 connected component | 248,817 ha |
| its span | lon −101.517..−100.083, lat 35.744..36.029 |
| GOES detected footprint, full event | 477,449 ha |
| IoU(block, footprint) | **0.392** |
| fraction of the block inside the fire perimeter | **0.823** |
| fraction of the fire's footprint that is non-burnable | **0.429** |

A 248,817 ha contiguous non-burnable block, 82% contained within the fire
perimeter, is the fire's own scar. Cropland does not align with a fire perimeter.

Consequences:

- **86% of the first-scan GOES seeds (36 of 42) land in it** and cannot spread.
  Of 42, only 6 sit on fuel 2. Across all 558 deduplicated detections, 53% are
  non-burnable.
- **43% of the fire's spread path is blocked**, so even correctly placed seeds
  cannot cross the domain.
- The §11 66× shortfall therefore measures seed placement and a corrupted fuel
  map, not Rothermel's ceiling. §12 quantifies how much of it was placement.

### Cause: wrfxpy's own scars mask (confirmed in the code)

This is not a Landfire artifact — wrfxpy puts it there deliberately, and you
identified it immediately:

- `src/fire_init/tools.py:363` — `NFUEL_CAT[FUEL_MASK] = no_fuel_cat`, with
  `no_fuel_cat = 14` by default (`tools.py:345`).
- `src/fire_init/process_perimeter_masks.py:298-300` supplies that mask from
  `scars_mask.pkl`, built at `:148` as `insidepastperims`.
- `src/forecast.py:343-347` carries `scars_mask.pkl` forward from
  `js.prev_forecast`, so a continuation forecast inherits every earlier scar.

The workspace has `perim1.pkl` and `perim2.pkl` at its top level, so it is a
**perimeter-initialized continuation forecast**: the already-burned area is
masked to no-fuel so WRF-SFIRE does not re-burn it. That is correct behaviour for
what it was built to do.

It also explains the other oddity in §10 — the 17,645 "ignition points" are
perimeter-derived, which is why they paint the observed fire.

**So this workspace cannot test cold-start forecast skill for any model**, and
that includes ForeFire. It was never meant to. Testing the low-latency pipeline
needs an *initial* forecast whose fuel map has no scars mask, ideally one built
at the true ignition (Smokehouse Creek started 2024-02-26, a day before this
workspace's nominal start).

A cheap guard for the driver: on load, report the fraction of the domain in the
no-fuel category and the fraction of the seeds landing in it. Both are one line
to compute and would have flagged this immediately instead of after a full run.

### A pipeline requirement falls out of this

Even on a clean fuel map, a 2.7 km ABI pixel is 27× coarser than the 100 m fuel
grid, so a detection centroid lands on a non-burnable cell often. The pipeline
must **snap each detection to the nearest burnable cell within its own pixel
footprint** (~1.35 km) and drop it if there is none. Without that, seeds are
silently inert — no error, no warning, just a fire that does not grow.

### How I found it, and a correction to my own method

I first mapped lat/lon linearly onto the 2400×2400 grid. That is wrong:
`make_FF_nc` builds the grid in Lambert Conformal metres and only stores the
inverse-projected corners in `BBoxWSEN`, so a linear lat/lon mapping skews the
lookup. Redoing it through the same LCC (`TRUELAT1/2`, `CEN_LAT/LON`,
`a=b=6370000`, spanning ±`nx*DX/2`) changed the first-scan figure from 81% to
86% — the conclusion survived, but the first number was not trustworthy.

I also verified the netcdf row order before trusting any of it: `make_FF_nc`
writes `NFUEL_CAT`/`ZSF`/`UF`/`VF` straight from WRF with no flip, and WRF's
`south_north` index 0 is southernmost, so `j` increases northward. A geographic
transect test was *not* decisive here; reading the code was.

## 12. How much of the shortfall was seed placement? Almost none

`goes_snap` repeats the §10 cold start with one change: each first-scan detection
is moved to the nearest burnable cell within its own pixel footprint (1,350 m),
and dropped if there is none. **33 of 42 were dropped** — the masked scar is
kilometres deep, not a road network, which is independent confirmation of §11.
Nine seeds survived. Same clock, same parameters, 1.1 min to run.

| h | goes_first (42 as detected) | polys | goes_snap (9 on burnable) | polys |
|---|---|---|---|---|
| 0.3 | 456 | 42 | 154 | 9 |
| 2.3 | 1,224 | 42 | 1,043 | 8 |
| 4.3 | 2,094 | 40 | 1,948 | 6 |
| 6.3 | 4,020 | 38 | 3,958 | 4 |
| 8.3 | **5,383** | 4 | **5,424** | 1 |

Nine correctly placed seeds reach the *same total* as forty-two as-detected ones:
**603 ha per seed against 128 ha**, a 4.7× improvement in seed efficiency. So the
33 seeds on masked fuel were genuinely inert — as expected from `e = 0` — and in
`goes_first` the six seeds that happened to land on fuel 2 did nearly all the
work.

**But the growth shortfall is unchanged.** 5,424 ha against an observed footprint
of 476,068 ha. Seed placement was not the binding constraint; the masked domain
is. The fronts spread until they reach the scar wall and stop, and 43% of the
fire's path is wall.

Front merging is confirmed again, and more cleanly than in §8: polygon counts
fall monotonically as fronts coalesce — 9 → 1 for `goes_snap`, 42 → 4 for
`goes_first` — while area rises throughout. No fronts are lost.

### What this leaves

Nothing in §10–§12 measures ForeFire's spread physics, because the domain cannot
support spread. What it does establish is the machinery: the corrected clock, the
detection-driven multi-ignition path, front merging at scale, and a runtime of
about a minute for an 8 h forecast. Those are the pieces the low-latency pipeline
needs, and they work. The physics question — whether Rothermel can approach the
rates you have seen on explosive fires — still rests on the independent 2.66 m/s
ceiling in fuel 2 from 09-09 §5b, and needs an unmasked initial forecast to test
properly.

## 13. ForeFire's geojson packs separate fronts as rings — a reader bug

Found while checking what looked like a catastrophic front loss in the
incremental run: between 7.3 h and 7.8 h its area fell from 94,624 ha to
5,857 ha and its feature count from 510 to 4. **Nothing was lost.** The vertex
count went *up* over the same step, 17,180 → 17,608.

`dumpMode=geojson` writes every separate front as a **ring of a single polygon**
and puts the count in `properties.numberOfPolygons`. At t = 28,800 s one feature
carried **472 rings**. GeoJSON convention says ring 0 is the shell and the rest
are holes, so `gpd.read_file` subtracted 471 fires as holes.

ForeFire signs the rings instead: a positive shoelace area is burning, a negative
one is a genuine unburned island. At that step there was exactly **one** real
island among the 472.

Fixed with `read_ff_geojson(path)` in `forefire.py`, which returns one row per
ring plus a `burning` column, so callers can keep fronts separate for KML or
difference the islands out for areas. `merge_geojson_to_kml` now uses it and
drops the islands — a hole drawn as a KML placemark reads as a fire in Google
Earth.

Verified:

| | rows | area |
|---|---|---|
| `gpd.read_file` on the 472-ring step | 4 | 5,857 ha |
| `read_ff_geojson` | 475 | **106,987 ha** |

Single-front output is unchanged — SINLAHEKIN `windReductionFactor_0.4` gives
136.9 / 169.1 / 205.3 ha either way — so nothing in §1–§6 is affected. **This bug
only bites multi-front runs, which means it arrived with §8 and would have
corrupted every KML the multi-ignition path produced.**

### Corrected results for §10

| variant | area at 8.3 h | peak IoU | fronts |
|---|---|---|---|
| `goes_first` — 42 seeds as detected | 5,726 ha | 0.012 | 35 |
| `goes_snap` — 9 seeds on burnable | 5,456 ha | 0.011 | 2 |
| `goes_incr` — fed every scan, 558 seeds | 90,495 ha (peak 106,869 at 7.8 h) | **0.218** | 429 |

So **the ongoing detection feed is worth about 19× in area and 18× in IoU** over
a cold start. That is the single most useful number for the pipeline: for a fire
like this, continuing to ingest detections matters far more than any parameter
choice measured in earlier sessions. And front merging is orderly throughout —
counts fall monotonically (558 → 429) while area rises.

Cost: `goes_incr` took 27.3 min against 1.2 min for `goes_first`, because each
restart serialises the whole front state and the `.ff` scripts reach ~1 MB. That
is the scaling limit to watch, not the front count itself.

## 14. The unmasked workspace — the mask was costing 5.5×

`wfc-SMOKEHOUSE_CREEK_..._points-2024-02-27_18:00:00-30` (no `_behave`) is the
initial forecast: **no `perim1.pkl`/`perim2.pkl`, and category 14 is 0.92% of the
domain against 5.2%.** It has 61 wrfouts plus saveouts through 2024-02-29, and
its `input.json` carries the same 17,645 points but across **21 distinct times**
(20:41:17–22:41:17) rather than all at one.

The seeds confirm the diagnosis of §11 outright:

| | masked `_points_behave` | unmasked `_points` |
|---|---|---|
| first-scan seeds dropped (no burnable fuel within 1,350 m) | **33 of 42** | **0 of 42** |
| all 558 seeds: snapped / dropped | 303 non-burnable | **2 snapped, 0 dropped** |

The detections were always on good fuel. The masked workspace was rejecting them.

Cold start, same 42 seeds, same clock, same parameters:

| | masked | unmasked |
|---|---|---|
| area at 8.3 h | 5,726 ha | **31,284 ha** |
| IoU | 0.012 | **0.063** |
| growth, 0.8–8.3 h | 634 ha/h | **3,960 ha/h** |
| shortfall vs observed 41,821 ha/h | 66× | **10.6×** |
| runtime | 1.2 min | 3.4 min |

So the scars mask alone accounted for a factor of **5.5** in area. The remaining
~10× gap is the real question, and it is now being asked on a domain that can
actually burn.

## 15. Estimating rate of spread from the detection feed

Your proposal: the distance from the first detection to the furthest later
detection, over the elapsed time, estimates the maximum spread rate, and that
gives a target for `windReductionFactor` or a `UF`/`VF` scaling. The method
works, but it needs two guards, and finding them was most of the work.

### Guard 1: cluster the first scan

Applied naively, the method gives **1.09 m/s** here — from the first-scan
centroid to the furthest last detection, 64.8 km in 16.5 h. That number is
meaningless, because the first scan is already **three fires spanning 96 km**
(§10). Distance from a single origin measures cluster separation, not spread.

Clustering scan 1 at 6 km linkage and assigning each later detection to its
nearest origin gives per-cluster rates instead.

### Guard 2: require connectivity, or the answer is nonsense

Per-cluster is still not enough. Cluster 2 (the small eastern one, 5 px) reported
**2.65 m/s** over 8.3 h and **9.96 m/s** with a p98 edge — physically absurd. The
test that exposes it: at each scan, check what fraction of a cluster's assigned
detections connect back to its origin through a chain of ≤4 km hops.

| cluster | px | connectivity (8.3 h) | reach | rate |
|---|---|---|---|---|
| 0 | 19 | 0.30 | 19.5 km | 0.750 m/s |
| 1 | 18 | **0.73** | 12.1 km | **0.197 m/s** |
| 2 | 5 | **0.05** | 90.3 km | 9.96 m/s ← junk |

Cluster 2's leading edge oscillates — 70 → 90 → 38 → 47 → 34 → 74 → 36 km — with
gaps of 24–52 km behind it. It is not a front. The nearest-origin rule was
handing it the *main* fire's eastern head, because cluster 2 happens to be the
easternmost origin. A front that genuinely spreads stays connected; a newly
associated patch appears at a distance with a gap behind it.

**Connectivity is the discriminator, and any pipeline use of this method needs
it.** Without it the method silently returns the diameter of the whole fire
complex divided by elapsed time.

### The window matters too

Over the shorter 3.3 h window connectivity is much better, because the complex
has not yet fragmented and NGFS has not yet swept in distant parts:

| cluster | connectivity (3.3 h) | rate |
|---|---|---|
| 0 | **0.80** | **0.827 m/s** |
| 1 | **0.95** | **0.500 m/s** |
| 2 | 0.20 | 9.96 m/s ← still junk |

So the defensible observed target is **0.5–0.83 m/s**, from the two coherent
clusters over 3.3 h. Both sit **below** ForeFire's 2.66 m/s ceiling in fuel 2
(09-09 §5b), so this is reachable by wind — which was not obvious beforehand.

Note this is an *average* rate along the fastest axis. It is a lower bound on the
instantaneous maximum, because a 2.7 km pixel locates the front only to within
half a pixel and the fire can burn between scans without lighting a new pixel.
Peak rates cannot be extracted from GOES this way — scan-to-scan differencing
gave 23–55 m/s, which is entirely association noise.

### Caveat on this fire specifically

Smokehouse Creek on Feb 27 is a mature complex, not a set of fresh ignitions, so
even the "coherent" clusters are parts of an existing fire. The method will work
far better on a genuinely new detection with a single origin, which is also the
case the low-latency pipeline actually serves.

## 16. `wind_scale`: scaling UF/VF in the netcdf

Added to `make_FF_nc(nc_path, out_path, wind_scale=1.0)` and driven from
`cfg['wind_scale']`, default 1.0.

**Why UF/VF and not `windReductionFactor`.** They multiply the same midflame
wind, so for the spread term they are interchangeable — but `wRF` also gates the
Andrews/Cruz/Rothermel wind cap at `Rothermel.cpp:217`, which is active only for
`wRF < 1.0`. Reaching a large multiple by raising `wRF` from 0.4 would cross 1.0
and switch the cap off, changing two things at once. Scaling `UF`/`VF` changes
only the wind, so the experiment stays clean.

Scaled files are cached under `nc_<grid_code>_w<scale>/` — the scale goes in the
**directory** name so the filenames the `.ff` scripts reference stay unchanged —
and the "reuse the workspace copy" shortcut is skipped when the scale is not 1,
since the workspace copy is unscaled.

### Two bugs found while validating it, both mine

1. I put the scaling in **`ff_ideal_nc`** instead of `make_FF_nc`. Both functions
   contain a near-identical `u_matrix = np.array(wrf_in.variables['UF']...)`
   block, and `ff_ideal_nc` is defined first. Worse, `ff_ideal_nc` has no
   `wind_scale` parameter, so the edit left a latent `NameError` on its fallback
   path. **If you are editing either of these, note that the UF/VF read appears
   twice in the file.**
2. The probe reported identical rates at every scale and I nearly read that as
   the Andrews cap saturating. It was not — the netcdfs were byte-identical
   (same md5 across the scale-1, `_w2` and `_w4` caches), which is the check that
   settled it. Compare the *inputs* before theorising about the physics.

## 17. The wind probe: what scaling actually buys

`wind_scale` swept over 1, 2, 3, 4, 6 on the unmasked workspace, 42 first-scan
seeds snapped to burnable fuel, 3.3 h horizon, template parameters otherwise.

| scale | mean \|w\| | area at 3.3 h | vs ×1 | cluster 1 edge | runtime |
|---|---|---|---|---|---|
| 1 | 3.70 m/s | 8,430 ha | 1.00 | 0.128 m/s | 0.7 min |
| 2 | 7.40 | 20,239 | 2.40 | 0.392 | 1.3 min |
| 3 | 11.10 | 34,364 | 4.08 | 0.401 | 1.9 min |
| 4 | 14.80 | 45,987 | 5.45 | **0.487** | 2.3 min |
| 6 | 22.19 | 69,173 | 8.21 | 0.557 | 3.6 min |

Observed target, cluster 1 over the same window: **0.500 m/s** at 95%
connectivity (§15).

### Findings

- **Area goes as `scale^1.17`** and is cleanly monotone. There is **no Andrews
  cap saturation** even at a 22 m/s mean wind, despite the cap being active at
  `wRF` 0.4. I claimed a plateau after seeing scales 1–3, where the 2→3 step
  happens to be flat; scales 4 and 6 disproved it. Three points were not enough.
- **`wind_scale ≈ 4.3` reproduces the observed leading-edge rate**, interpolating
  between 0.487 at scale 4 and 0.557 at scale 6.
- Cost scales gently: 3.6 min at scale 6 against 0.7 min at scale 1, so this is
  affordable inside a low-latency budget.

### But 4.3× is not a wind correction

A mean *fire-level* wind of ~16 m/s is not physically credible. `UF`/`VF` are
coupled winds near flame height, which are normally **lower** than the 10 m wind,
and 16 m/s at flame height would imply a 10 m wind well above what this event is
reported to have had. So the honest reading is:

**~4.3× is the total multiplier the spread formulation needs, expressed in the
only units this experiment varied.** It is not evidence that WRF's winds are 4.3×
too low. Some plausible split: WRF fire-level winds genuinely low in complex
coupling, `pSAF` at 0.6 being conservative, moisture too high, and — the piece
Rothermel structurally cannot supply — spotting, which this fire certainly did
(09-09 §5b).

Treating the whole 4.3× as wind would be calibrating one error with another. It
is a usable *effective* tuning for low-latency work, and should be labelled that
way rather than as a wind bias correction.

### What to do next with it

1. Split the multiplier: sweep `pSAF` and `Md` at `wind_scale` 1 on this unmasked
   workspace and see how much of the 4.3× they can absorb at physically sensible
   values. Both scale `R` more directly than wind does, and 09-08/09-09 showed
   they transfer between fires whereas `wRF` does not.
2. Re-derive the target on a fire with a **single** origin. Every rate here comes
   from one cluster of a mature complex; the method in §15 will be much better
   constrained on a genuinely new detection, which is also the pipeline's real
   case.
3. Check the equivalence numerically: `wind_scale` 4.3 at `wRF` 0.4 implies an
   effective midflame factor of 1.72, which is unreachable via `wRF` alone
   without crossing the cap gate at 1.0. That is the concrete reason the knob was
   added to the netcdf rather than the script.

## 18. The NWS summary settles the wind: a double reduction, not a WRF error

You supplied the NWS Amarillo event summary. The relevant facts for Feb 27:

- morning: west to southwest winds **25–35 mph**, gusts exceeding 58 mph by late
  morning
- afternoon: frequent gusts **65–70 mph**, some areas sustained **40–45 mph** for
  a few hours
- relative humidity **12–15%**, afternoon temperatures in the 70s F
- a cold front pushed in from the north; **"fires that once were moving from west
  to east soon started to move from north to south behind the front"**
- roughly **9 named fires**, several merging; Smokehouse Creek finished at
  **1,058,482 acres = 428,350 ha**
- a forward flank **~100 miles long**

### WRF's winds are right in both direction and magnitude

**Direction and frontal timing:** averaged over the fire footprint, WRF holds
westerly (252–258°) through 23:30Z, swings NW at 00:30Z, and is northerly
(351–12°) from 01:00Z onward. That is the described frontal shift, in the right
place in the forecast window. WRF's synoptic evolution is correct.

**Magnitude:** the namelist settles it.

```
fire_wind_log_interp = 1
fire_use_windrf      = 2
windrf = 0.36, 0.36, 0.44, 0.55, 0.42, ...
```

`fire_use_windrf = 2` means **WRF-SFIRE has already applied the fuel-dependent
wind reduction factor**, 0.36 for category 2. So `UF`/`VF` are already midflame
winds. The implied 10 m wind is

    4.15 / 0.36 = 11.5 m/s = 25.8 mph

which lands exactly in the NWS's "sustained 25 to 35 mph". **WRF's wind magnitude
is not the problem.**

### The problem is that ForeFire reduces the wind a second time

ForeFire's `windReductionFactor = 0.4` multiplies `UF`/`VF` again, so the total
reduction from the 10 m wind is

    0.36 (WRF-SFIRE) x 0.40 (ForeFire) = 0.144

where it should be 0.36. The effective midflame wind is **1.48 m/s** on a day
that warranted 4–6 m/s. That is a factor of **2.5**, i.e. `1/0.4`, and it is most
of the 4.3× the wind probe needed (§17).

**Correction, and it supersedes §17's framing.** I wrote there that 4.3× "is not
a wind correction" because 16 m/s at flame height is not credible. That reasoning
was wrong: at `wind_scale` 4.3 the netcdf carries 15.9 m/s, but ForeFire then
applies `wRF` 0.4, so the midflame wind is 6.4 m/s — entirely credible for a
25–45 mph day. I compared the pre-`wRF` field against a flame-height expectation.

**The right fix is `wRF = 1.0`, not a wind scaling** — because the reduction is
already in the data. But `wRF ≥ 1.0` switches off the Andrews cap gate
(`Rothermel.cpp:217`), so the clean way to express it without changing two things
is `wind_scale = 2.5` with `wRF` left at 0.4. Those are numerically equivalent in
the spread term and differ only in the cap.

### The residual after the double reduction is ~1.7x

4.3 / 2.5 = **1.7**, and that is the part the fuels have to carry — which is
where you pointed.

## 19. The fuels: your hypothesis measured

Both fuel corrections you proposed are supported by the NWS text independently of
any tuning argument.

**Loading.** *"A good growing season through spring and summer of 2023 led to
decent grass fuel loading compared to the previous couple of years."* The domain
is **86% NFUEL_CAT 2**, "timber (grass and understory)" — a model presuming a
timber overstory the Texas panhandle largely lacks, and built from a static map
that cannot know 2023 was a good growing year. Anderson 3, tall grass, has 2.7×
the bed depth, a lower packing ratio (so a larger `phiV` from the same wind) and
`me` 0.25 against 0.15.

**Moisture.** *"Two weeks without precipitation"*, dormant vegetation, RH 12–15%
at 70s F. Simard EMC for those conditions is **0.030–0.041**, against the fuel
table's default **Md = 0.10**. But moisture is a weaker lever than it looks here:
`Etam` only rises 0.501 → 0.658 from Md 0.10 → 0.03 in fuel 2, i.e. **1.31×**,
because 0.10 is already well below `me` 0.15. The 09-08 knee is above this range,
not below it.

`fuel_remap` was added to `make_FF_nc` (and `cfg['fuel_remap']`) to test the type
change; scaled/remapped netcdfs cache under `nc_<grid_code>[_w<scale>][_f2to3]/`.

### Measured, all at wind_scale 1 unless noted, 3.3 h, same seeds and clock

| case | fuel | Md | wind | area | vs base |
|---|---|---|---|---|---|
| `fuel_base` | 2 | 0.10 | 1 | 8,430 ha | 1.00× |
| `fuel_tall` | **2→3** | 0.10 | 1 | **19,625 ha** | **2.33×** |

Tall grass alone gives **2.33×**. With the 2.5× double-reduction correction that
is 5.8×, against the 4.3× needed — so between them the levers more than cover the
gap, and the question becomes how to split it rather than where to find it.

### Where the ForeFire fuel table comes from, and what is actually wrong

I first called the ForeFire/WRF-SFIRE table mismatch a defect and suggested
raising category 2's load as "correcting a transcription". That was too broad.
You identified the cause: **`fuelstrans.csv` is mapped from the 40-category
Scott & Burgan models into the closest of the 13 Anderson categories**, so it is
a translation, and its parameters are not expected to equal Anderson's.

The evidence supports that, and narrows the problem to one cell.

**`fuelstrans.csv` shares a source with WRF's `namelist.fire`:**

| column | agreement with `namelist.fire` |
|---|---|
| `me` | **identical** for all 11 comparable categories |
| `e` | `fueldepthm` rounded, all 9 (0.305→0.3, 0.762→0.8, 1.829→1.8, 0.061→0.1) |
| `Sigmad` | `fgi` within ±20% for **7 of 9** |
| `sd` | matches neither — binned to six levels (7800/6500/5500/4500/4000/3500) |

The `sd` binning is the translation fingerprint: the values are SI (1/m) but do
not equal Anderson's characteristic SAV in either 1/ft or 1/m, and they collapse
13 distinct values into six. That is a deliberate lossy mapping, not an error.

**There is also a second table, `fuels13.csv`**, which *is* a direct
transcription of `namelist.fire` — `Sigmad` = `fgi` (0.166, 0.896, 0.674, …),
`e` = `fueldepthm` exactly, `sd` = `savr` exactly. But it is **not usable as-is**:
its `sd` is in **1/ft** while every other column is SI, so SAV would be 3.28×
too small, and row 13's load is `1.00E-07`. `fuelstrans.csv` is the right file to
be using.

**The one thing that does look wrong is category 2's load:**

| | value |
|---|---|
| `fuelstrans.csv` | **0.400** |
| `fuels13.csv` | 0.896 |
| WRF `fgi` | 0.897 |
| ratio to `fgi` | **0.45×**, against 0.71–1.41× for every other category |

Category 9 is the only other outlier (1.41×). Since the SMOKEHOUSE domain is
**86% category 2**, this single cell carries the domain. Whether correcting it
helps is not obvious from the formulation — raising the load raises reaction
intensity but also the bulk density in `R0`'s denominator — so it was measured
rather than argued (§20).

**Caveat on the whole comparison:** because `fuelstrans.csv` is a 40→13
translation and WRF-SFIRE uses Anderson 13 directly, the two models are not
simulating identical fuels even given the same `NFUEL_CAT` map. That is a
structural limit on ForeFire-vs-WRF-SFIRE comparison, including the pRes tuning
idea from 09-09, and it is not something a single edit fixes.

## 20. Corrections: the ring reader, and the Andrews cap

### `read_ff_geojson` was wrong twice, and the second version is committed broken

`daecd00` decided burning-vs-island by **ring winding** (positive shoelace area =
burning). That is not reliable: ForeFire does not wind merged fronts consistently,
so a **31,079 ha burning front came through negative**, was taken for an unburned
island and differenced out. The symptom was area oscillating between steps in the
`wfix_tall` run:

    3,473 -> 11,837 -> 22,234 -> 33,257 -> 9,003 -> 45,028 -> 12,632 ha

Every collapse coincided with exactly one "island" appearing, and at those steps
the island was 4.7–5.1× the size of the largest "shell" — which is the tell.

**The correct rule is nesting depth under the even-odd convention**: a ring
enclosed by an even number of other rings is burning, odd means an unburned
island inside a burn. That holds whatever the winding. `read_ff_geojson` now
computes containment depth via representative points with a bounding-box
prefilter. With it, `wfix_tall` grows monotonically to 47,488 ha.

**This is uncommitted and supersedes the committed version.** Anything scored
between `daecd00` and this fix that had merged fronts is suspect.

### Numbers that change

| result | was | now |
|---|---|---|
| `wfix_tall` (wind fix + tall grass) | 12,632 ha | **47,488 ha** |
| `goes_incr` final, masked, 8.3 h | 90,495 ha | **120,763 ha** |
| `goes2_first` final, unmasked | 31,284 ha | 31,305 ha |
| everything single- or few-front | — | unchanged |

Sections 1–6 are unaffected: SINLAHEKIN never had inner fronts, and its areas are
identical either way (136.9 / 169.1 / 205.3 ha).

### The Andrews cap is NOT binding — measured at last

I invoked the cap three times this session as a reason to prefer `wind_scale`
over `wRF = 1.0`, and once claimed it explained SINLAHEKIN's flattening (it did
not). It is now measured directly. Both configurations give the same midflame
wind and differ *only* in the cap gate at `Rothermel.cpp:217`:

| | midflame | cap | area |
|---|---|---|---|
| `wind_scale` 2.5, `wRF` 0.4 | ×1.0 | **active** | 27,691 ha |
| `wind_scale` 1.0, `wRF` 1.0 | ×1.0 | **off** | 27,868 ha |

**0.6% apart.** The cap does not bind on this fire, the two knobs are
interchangeable, and **`wRF = 1.0` is the simpler and more honest fix** than a
wind scaling. `wind_scale` remains useful for sensitivity work, but it is not
needed to express the double-reduction correction.

## 21. Where the correction should come from: the balance sheet

Target: **ROS 2.34×** (area 5.46×, ~46,000 ha at 3.3 h), from the
connectivity-guarded leading edge of the best-connected first-scan cluster (§15).

| correction | area | ROS | basis |
|---|---|---|---|
| category 2 load 0.400 → 0.896 | 1.21× | 1.10× | correctness — outlier cell |
| Md 0.10 → 0.035 | 1.83× | 1.35× | **observed** — NWS RH 12–15% |
| `wRF` 0.4 → 1.0 | 3.28× | 1.81× | **correctness** — `namelist.fire` |
| `wRF` fix + Md 0.035 | 4.77× | **2.18×** | both justified from inputs |
| `wRF` fix + Anderson 2→3 | 5.63× | 2.37× | needs a reclassification |
| Anderson 2→3 + Md 0.035, no wind fix | 4.05× | 2.01× | fuels only |

**The two corrections that need no hypothesis — the `wRF` double reduction and
moisture from observed RH — reach 93% of the target on their own.** The residual
7% is well inside the target's uncertainty, which derives from one cluster of a
complex the NWS confirms was merging ("roughly 9 named wildfires with some of
these merging"). So a fuel reclassification is **not required** to explain this
fire, which matters given it is the least defensible of the levers.

Levers are **sub-multiplicative** when combined: tall+dry+wind×2 measures 6.92×
where multiplying the singles gives 9.72×. Combinations must be run, not inferred.

## 22. Recommendation: correct the inputs, leave the fuel categories alone

Your position, and it is the right one: *"It's good to avoid any unjustified fuel
category changes. That's black magic that might work for one fire but fail in
every other circumstance."* The measurements support it — the fuels are not
needed for this fire.

### What to change

| # | change | basis | ROS |
|---|---|---|---|
| 1 | `windReductionFactor` 0.4 → **1.0** | `namelist.fire` has `fire_use_windrf = 2` with `windrf = 0.36`, so WRF-SFIRE already reduced the wind. ForeFire reducing it again is a double count. | 1.81× |
| 2 | `Md` 0.10 → **~0.035** | Simard EMC from the NWS's RH 12–15% at 70s F, vegetation dormant, two weeks dry. | 1.35× |
| | **together** | | **2.18× of the 2.34× target (93%)** |

Neither is a tuning choice. (1) is a code correctness fix demonstrable from the
namelist without reference to any fire; (2) is an observed input replacing a
placeholder. The residual 7% sits inside the target's own uncertainty, which came
from one cluster of a complex the NWS confirms was merging.

Use `wRF = 1.0` rather than `wind_scale` 2.5: they are equivalent (§20), and the
Andrews cap does not bind, so the simpler expression is the honest one.

### What NOT to change, and why

- **Fuel category remapping (Anderson 2→3).** Gives 1.53× and is not needed. The
  physical argument is real — the NWS notes a good 2023 growing season, and
  category 2 presumes a timber overstory the panhandle lacks — but it is a
  judgement that would not transfer to another fire. `fuel_remap` was implemented
  in `make_FF_nc` and is **deliberately left uncommitted**; the experiments are in
  the handoff, the facility is not in the repo.
- **Category 2's fuel load.** I called 0.400 an outlier cell and recommended
  raising it to 0.896 as "correcting a transcription". That was overstated twice
  over. It is ambiguous:

  | | value |
  |---|---|
  | Anderson FM2 **1-hour** load | 0.448 |
  | Anderson FM2 **total** load (= WRF `fgi`) | 0.897 |
  | `fuelstrans.csv` | **0.400** |

  8 of 9 categories in `fuelstrans.csv` follow the **total** convention and
  category 2 alone follows **1-hour** — so it is inconsistent with its
  neighbours, but 0.400 is exactly what a 1h-only reading gives, so it may be
  intentional. **Measured effect: 1.10× in ROS.** Not worth guessing over.

### Where `fuelstrans.csv` comes from (confirmed)

You identified it: mapped from the 40-category Scott & Burgan models into the
nearest of the 13 Anderson categories. `wildfire_ROS_models/fuels_database.py`
(`CB2005_t7_csv`) carries the S&B table and makes the fingerprint clear —
category 2's `sd` 6500, `e` 0.30, `me` 0.15 is essentially **GR2 fully cured**
(6562, 0.305, 0.15), and the binned SAV levels 7800/6500/5500/4500 map onto S&B's
2200/2000/1800/1500 1/ft.

One observation worth keeping, not acting on: GR2 cured carries a load of
**0.247 kg/m²** at that depth, against `fuelstrans`' 0.400 — a bulk density of
0.809 vs 1.333, so category 2 as shipped is denser than any cured S&B grass
model, and bulk density sits in `R0`'s denominator. That is a lead for whoever
reconciles the table properly, with a source, rather than a change to make now.

**Structural caveat that no edit fixes:** because `fuelstrans.csv` is a 40→13
translation and WRF-SFIRE uses Anderson 13 directly, the two models are not
simulating identical fuels even given the same `NFUEL_CAT` map. Any
ForeFire-vs-WRF-SFIRE comparison inherits that, including the pRes tuning idea
from 09-09.

Grass-fire behaviour is an active research area and the loading/curing question
here sits squarely in it; that is the reading to do before touching the table.

## 23. Open items

Carried from 09-09 §7, plus:

1. 09-09 §5f is wrong on timing and should be read with §2 above.
2. **Answered by §6:** the elongation gap survives the origin correction
   (1.42-1.80 vs observed 2.54, against 1.23-1.76 before). It is structural, not
   a misalignment artifact.
3. `attr_FireBehaviorGeneral` reports short-range spotting for this fire, which
   ForeFire structurally cannot produce standalone (09-09 §5b). Some residual
   shape error is expected to be unreachable by any parameter.
4. Scratch scripts remain ephemeral: `shape_compare.py`, `overlay_plot.py`,
   `grid_ign2.py`, `grid_plot.py`. `sweep_analyze.py` is still the one most worth
   promoting into the repo.

---

### Added this session

5. **Recompute the absolute numbers on the corrected clock.** 09-09 §5f
   (6.2× vs WRF-SFIRE), 09-10 §2 (1.74×/1.83×) and the pSAF ≈ 0.4 calibration all
   ran with a ~3 h head start (§7). Sensitivities and the matched-area IoU work
   are unaffected; those three need redoing before they are quoted.
6. **Use the unmasked workspace for anything Smokehouse.** `_points` is the
   initial forecast; `_points_behave` is a continuation with the scar masked to
   no-fuel and is unusable for forecast skill (§11, §14). Check any other
   retrospective case the same way: `perim1.pkl`/`perim2.pkl` present, or
   category 14 above ~1%, means masked.
7. **Done, and worth keeping:** the no-fuel diagnostic. `goes_run2.py` reports
   the category-14 fraction of the domain and the snapped/dropped seed counts at
   load. Two lines; promote them into `forefire.py` rather than leaving them in
   a scratch script.
8. **Done:** seeds are snapped to the nearest burnable cell within their own ABI
   pixel footprint (1,350 m). Keep it even on clean fuel maps — a 2.7 km pixel is
   27× coarser than the 100 m fuel grid. On the unmasked domain it moves only 2
   of 558 seeds, so it costs nothing when it is not needed.
9. **Build a persistence filter before using `feature_tracking_id` back-tracing**
   (§10). Without one it seeds phantom fires on gas flares — 30 such pixels sit
   within 10 km of this fire, detected on 218 of the day's 274 scans at steady
   200–350 MW.
10. **The detection feed is the biggest lever measured so far** — 19× in area,
    18× in IoU over a cold start (§13). Ingesting detections continuously
    outranks every parameter choice from 09-08/09-09. Design the pipeline around
    that.
11. **Watch restart serialisation, not front count.** 558 fronts are fine
    numerically, but each restart writes the whole front state and the `.ff`
    scripts reach ~1 MB, taking the incremental run to 27.3 min against 1.2 min
    for the cold start. Front *count* was never the limit; state I/O is.
12. **Partly answered (§16, §17): the residual gap needs an effective 4.3×.**
    `wind_scale` is implemented and swept; area goes as `scale^1.17` with no
    Andrews saturation to 22 m/s, and `wind_scale ≈ 4.3` reproduces the observed
    leading-edge rate. But 16 m/s at flame height is not credible as a wind
    correction, so 4.3× is the multiplier the *formulation* needs, not a wind
    bias. **Next: sweep `pSAF` and `Md` at `wind_scale` 1 on the unmasked
    workspace to see how much of the 4.3× they absorb at sensible values**, then
    re-derive the target on a single-origin fire.
13. **All committed:** `d6b03a8` timing fix + multi-point ignition, `daecd00`
    the geojson ring reader, `282fa86` `wind_scale`. Handoff sections 7–17 are in
    those commits.
14. New scratch scripts, still ephemeral: `goes_run.py`, `goes_run2.py`,
    `goes_snap.py`, `goes_verify3.py`, `goes_verify4.py`, `goes_plot.py`,
    `goes_fuel.py`, `ros_estimate2.py`, `wind_probe.py`. The two worth promoting
    are `ros_estimate2.py` (the connectivity-guarded ROS estimator from §15,
    which the pipeline needs) and `goes_verify3.py` (scoring against GOES pixel
    polygons).
15. **`goes2_incr` was still running** when the session ended — 558 seeds on the
    unmasked workspace, considerably slower than the masked run because the fire
    actually spreads. Score it with `goes_verify4.py`.
