# ForeFire session handoff — 2026-09-22

Branch `james_ngfs`, **nine commits this session, none pushed**:

```
af212f4  Add the WindNinja against WRF-SFIRE wind comparison
1539a22  Add a driver to run both ensembles to a common end time
e5e5cb3  Record the 2026-09-22 session: the first complex-terrain case
1172295  Add Dome to the case table; roughness falsified, B/A is not
531de4a  Defer Cotton 2 until the cache covers its wind reversal
9302b05  Add perimeter scoring against an observed IR perimeter
77171c4  Score at a target time, not only the observation timestamp
eb5b8ff  Record Tabor: the first run with no WRF-SFIRE in the chain
```

Follows `SESSION_HANDOFF_forefire_2026-09-21.md`, whose §10 set the goal: **test the
GRIB-driven path in more complicated terrain.** Done, on two fires — and a third that
went further, running with no WRF at all.

**Headline 1, complex terrain:** **ForeFire tolerates WindNinja's direction error and
is sensitive to its speed bias**, and the speed bias does not go one way. Silver
(Selkirk Mtns) ran 24% slow and produced a fire 28% *smaller* than the coupled run;
Dome (Yosemite), rougher still, ran 18% fast and produced one 16% *larger*. Both times
the area ratio tracked the speed ratio, and both times a 36-49 deg direction RMS moved
the ensemble centroid under 520 m on a multi-kilometre fire. **Terrain roughness
predicts none of it** (§1).

**Headline 2, the goal itself:** Tabor ran **end to end from an ignition point, a time
and GRIB files, with no WRF-SFIRE anywhere**, and was scored against a real IR
perimeter (§6). That is 09-18 §11 delivered. Both measurement tools 09-18 asked for now
exist (§2), so its open item 6 is fully closed.

---

## 1. Complex topography — the case table

**This is the table to grow.** One complex-terrain case is not a result, it is a
hypothesis. Terrain roughness is `ZSF` std over the fire grid, which is the number
09-18 §4 used to separate the converged from the unconverged regime.

    fire        ZSF std   wind B/A  dir RMS   anomR   area HRRR/WRF  area/wind  centroid
    Dry River     26.7 m    1.04       3.7     0.129   (not run)         -          -
    Union         ~20 m*    -          -       -        0.99 +           -          -
    Red Bank      31.5 m    -          -       -        1.035            -          -
    Silver       364.7 m    0.765     49.1     0.109    0.717          0.937      514 m
    Dome         504.1 m    1.175     35.9     0.396    1.164          0.991      455 m
    Cotton 2     184.4 m    1.119     54.3     0.20*    1.082          0.967        -

    * Union: relief 21-93 m; std not computed.  Cotton 2 anomR is a rough mean over
      hours with R >= 0.3; its wind reversal makes several hours unusable.
    + Union is an ensemble mean against a single WRF chained-track perimeter, not
      ensemble against ensemble.  Red Bank, Silver and Dome are ensemble against
      ensemble, each at a single common valid time.

**Roughness is falsified as the predictor; `B/A` is not.** The first version of this
table guessed that rougher terrain would mean a lower area ratio. **Dome breaks
that** — it is the roughest case and errs the *other* way, 1.16 against Silver's
0.72, with *better* direction agreement and by far the best anomaly correlation of
any case including flat ground. Nothing in the wind columns is monotonic in `ZSF`
std.

What does hold, on two complex cases sitting either side of 1.0:

- **the area ratio tracks the speed ratio**, area/wind 0.937 and 0.991;
- **the direction error does not matter** — 49 deg and 36 deg RMS both moved the
  ensemble centroid only 455-514 m on multi-kilometre fires.

**How strong that is: n=2.** Dome's area ratio was predicted at ~1.10 before the
ensembles ran, from Silver's own area/wind ratio, and came out 1.164 — right sign,
right magnitude. But two points always fit a line, and the prediction borrowed its
constant from one of them. The *sign* following `B/A` is established; the
*coefficient* is not.

**Why it would matter if it survives.** `B/A` is measurable from **one netcdf and
one wrfout** with `wn_vs_wrf.py` (§2) — minutes, against hours for an ensemble. If
it predicts the area ratio, a GRIB-driven forecast has a knowable, bounded error
that can be estimated before it is run, and possibly corrected. That is the single
most valuable thing left to test.

Gaps to fill, cheapest first: `B/A` and direction statistics for **Union and Red
Bank**, whose runs already exist and whose comparisons predate the tool. Two flat
points with wind statistics would show whether area/wind stays near 1 at `B/A`
near 1, which is the weakest part of the claim.

## 2. The measurement tool

`src/ingest/wn_vs_wrf.py`, commit `af212f4`. Compares a WindNinja-derived ForeFire
netcdf against the wrfout's `UF`/`VF`. Both carry `windrf` — the wrfout's by
construction, the netcdf's because `apply_windrf` put it there — so nothing is
divided out, which removes the category-14 masking problem 09-18 §5 had going the
other way.

WRF's fire grid is Lambert and the netcdf's is UTM, so `FXLONG`/`FXLAT` are
projected into the netcdf's UTM and the netcdf sampled bilinearly there; statistics
are computed on WRF's own grid. Per hour it prints the speed ratio, direction bias
and RMS, the fraction of cells more than 90 deg out, the speed-anomaly correlation,
the peak ratio, and **the circular resultant length R next to them** — because
below about 0.3 the direction columns mean nothing. That is the trap 09-18 §2
recorded, and the number that detects it is now printed beside the numbers it
invalidates instead of being left to memory.

This closes 09-18 §10 item 6 for the wind half. `src/ingest/ff_score_perim.py`
(`9302b05`, `77171c4`) closes the perimeter-scoring half, so that item is now fully
retired. It reports area, obs ratio, IoU and elongation, and validates against 09-18's
own figures: on the Dry River perimeter it returns 2078.5 ha and elongation 2.90
against that handoff's recorded 2078.5 and 2.90.

## 3. Silver — what complex terrain did to the winds

`Silver_2026-09-20_18_00_00_68A6E77D-AC4A-4E0B-A7FB-2CB72D7ED142`, Selkirk Mountains,
northern Idaho. Ignition 2026-09-20 19:51:19Z at 48.636665, -116.661667. 30 km at
25 m. **`ZSF` std 364.7 m, relief 637-2300 m** — the Wind River regime (472 m,
measured as unconverged), not Dry River (26.7 m).

Solved at **25 m, the fire resolution**, with no `--ascii_out_resolution`
interpolation: ~165 s per WindNinja solve, ~3 min per step including warps and the
netcdf write, 26 steps. 09-21 §10 said not to carry `--wn-mesh 50` over and that
was right — at 30 km the budget allows 25 m (1.44 M cells) but not 12.5 m (5.76 M,
fails).

All 26 hours, every one with R >= 0.59 so none is the light-wind artefact:

    B/A            mean 0.765   range 0.614 to 0.956
    dir RMS        mean 49.1    range 28.1 to 82.2 deg
    anomaly corr   mean 0.109   range -0.130 to 0.406
    peak B/A       mean 1.50    range 0.92 to 2.00

**Three things, in decreasing confidence:**

1. **WindNinja is systematically slow and never once exceeds WRF.** Flat ground gave
   1.04; here 0.765.
2. **Direction RMS is 13x the flat-terrain value** (49.1 against 3.7 deg).
3. **Anomaly correlation is 0.109 — no better than flat ground's 0.129.** Terrain
   does *not* give the two models shared spatial structure. An early single-hour
   reading of 0.406 suggested otherwise and was an outlier; do not repeat that
   inference from one hour.

**A clean diurnal signature.** Direction bias sweeps a full arc — near zero in the
evening, +33 deg by 05Z, positive through the night and morning, -25 deg by 16Z,
back toward zero by 20Z. Peak ratio tracks it: 2.0 in the evening, ~1.0 overnight,
2.0 again the next evening. WindNinja is at its best in the well-mixed evening and
degrades as the boundary layer decouples, which is what a mass-conserving solve with
no stability or thermal physics would do.

**Hypothesis tested and rejected: `--diurnal_winds` does not fix it.** It defaults
off and `forefire_grib.py` never passes it, which looked like the obvious
explanation. Enabled at 04Z, the worst hour, it helps speed slightly and **makes
direction worse**:

    04:00Z          B/A     bias   dirRMS   >90     anomR   peak B/A
    non-diurnal    0.681   +16.7    47.5    6.4%   -0.104     1.01
    diurnal        0.731   +18.3    55.3   10.7%   +0.017     1.21

Do not turn it on for this. The overnight divergence remains unexplained.

## 4. Silver — what it did to the fire

Both ensembles at a common valid time of 2026-09-21T21:00Z. The orchestrator
computed that end as `min(last wrfout, last GRIB)` and trimmed both tables to it;
here the two data ends coincided naturally.

    source                    members   mean      median     sd     CV    spread
    WRF-SFIRE (30-min)           51     2420.8    2387.2   238.7   9.9%   1.62
    HRRR+WindNinja (hourly)      26     1734.5    1699.4   134.3   7.7%   1.43

    HRRR/WRF  mean 0.717   median 0.712
    mean centroid offset  514 m

**The two wind errors behave completely differently, and this is the useful part:**

- **The 49 deg direction scatter did essentially nothing.** Centroids 514 m apart on
  a fire of ~2.5 km radius — the same fire in the same place. ForeFire integrates
  over direction noise, as it did on flat ground against an uncorrelated field.
- **The 24% speed deficit propagated almost one-for-one into area.** `B/A` 0.765 in,
  area ratio 0.717 out. That is the speed deficit working through Rothermel, not a
  steering or shape effect.

**Why that matters beyond this fire:** a systematic speed bias is the kind of error
that can be characterised across cases and, in principle, corrected. Uncorrelated
direction scatter could not be. So the table in §1 is worth filling in: if `B/A`
predicts the area ratio across fires, the GRIB-driven path has a known, bounded
error rather than an unreliable one.

Also note the HRRR ensemble is **tighter** than WRF's here (CV 7.7% against 9.9%),
the opposite of Red Bank, where 09-21 §6 found it more dispersed near the ignition.
That near-ignition dispersion therefore does **not** generalise.

## 5. Dome — the second complex case, erring the other way

`DOME_2026-09-15_17_00_00_AACFD673-4C1B-4CDF-B9DD-128E34BD7272`, Yosemite. Ignition
2026-09-15 18:11:53Z at 37.568215, -119.615585. 30 km at 25 m, solved at 25 m.
**`ZSF` std 504.1 m, relief 571-3006 m** — the roughest case measured, above the
Wind River reference (472 m) that 09-18 §4 found unconverged at any mesh.

California, which per JH matters operationally: other WRF-SFIRE developers have
CalFire relationships and performance across California is always of interest. So
the hardest terrain and the highest operational scrutiny coincide.

Winds, all 28 hours, every one with R >= 0.76:

    B/A            mean 1.175   range 0.795 to 1.494
    dir RMS        mean 35.9    range 21.2 to 76.1 deg
    anomaly corr   mean 0.396   range 0.250 to 0.562
    peak B/A       mean 2.01    range 1.15 to 3.21

**Opposite to Silver on speed, better on everything else.** WindNinja runs *fast*
here (1.175 against Silver's 0.765), its direction RMS is lower (35.9 against 49.1),
and its anomaly correlation of 0.396 is by a wide margin the best of any case
measured — better than flat ground's 0.129. The plausible reading is that where
terrain forcing genuinely dominates, both models respond to the same topography, and
Silver's mid-range relief is the awkward regime where WRF's boundary-layer structures
matter but WindNinja's terrain response does not yet dominate. **That is a hypothesis
and needs the cases in between.**

The diurnal signature repeats from Silver, offset: daytime hours run 1.32-1.49,
overnight hours fall to 0.795-1.012. Including the night is what pulls the full-window
`B/A` down from a daytime 1.438 to 1.175 — **so always use the full window against
the area ratio, never the daytime figure.**

Ensembles, both at a common valid time of 2026-09-16T21:00Z (HRRR data ran an hour
longer and was trimmed):

    source                    members   mean      median     sd     CV    spread
    WRF-SFIRE (30-min)           54     3370.9    3366.9   202.5   6.0%   1.31
    HRRR+WindNinja (hourly)      27     3922.8    4020.1   387.2   9.9%   1.50

    HRRR/WRF  mean 1.164   median 1.194
    mean centroid offset  455 m

**The area ratio was predicted at ~1.10 before the run** — Silver's area/wind of
0.937 applied to Dome's 1.175 — and came out **1.164**. Right sign, right magnitude.
See §1 for how much weight that carries: the prediction borrowed its constant from
the only other point.

Note the ensemble spread inverts again: the HRRR ensemble is *more* dispersed here
(CV 9.9% against 6.0%) where on Silver it was *tighter* (7.7% against 9.9%), and on
Red Bank more dispersed. **Ensemble spread behaviour does not generalise across
cases** — three fires, three different orderings.

## 6. Tabor — the first run with no WRF-SFIRE at all

**This is the 09-18 §11 goal reached end to end:** a fire forecast from an ignition
point, an ignition time and GRIB files, with no WRF run anywhere in the chain, scored
against a real IR perimeter.

`Tabor_2026-09-20_06_00_00_F1623A08-608E-4150-816C-C14A302C5867`, Ozarks, Missouri.
Ignition 2026-09-20 07:40:30Z at 36.780685, -92.112183. **DEM std 36.6 m** — a
Red Bank-class flat case. Not seen by GOES; VIIRS only.

**The workspace was created by hand**, in
`/data/jhaley/new_wrfxpy/wrfxpy/wksp/wfc-Tabor_..._-33/`. The job file in
`new_wrfxpy/wrfxpy/jobs/` already carried everything `read_input` needs — `start_utc`,
`ignitions.1[0].time_utc` and `.latlon`, `grid_code` — so it was copied to
`input.json` with `end_utc` extended from the job's 2026-09-20 17:00 to 2026-09-21
17:00, to reach the perimeter. **Stopping at the job's own window would have made the
fire unscoreable.** Domain 31 km at 25 m, from the job's `domain_size` 31 and
`subgrid_ratio` 40. 37 gap-free GRIBs, 35 netcdfs, 35 members.

### Scoring a suppressed fire

Tabor was fought: `attr_FireStrategyFullSuppPrcnt` 100, `attr_PercentContained` 99,
behaviour "Minimal / Flanking / Smoldering", **contained 2026-09-20 23:59** — 17 hours
*before* the perimeter's own timestamp. Per JH an estimated **98% of fires have
suppression applied**, which no model here represents.

    observed                      158.9 ha   elongation 1.93   (99.9% of stated acreage)
    at containment  09-21 00:00    70.1 ha   obs_x 0.44   IoU 0.225   elong 1.16
    at observation  09-21 18:00   255-302    obs_x 1.60-1.90  IoU 0.29-0.36  elong 1.01-1.26
    crosses 158.9 ha at 09-21 09:00, about 9 h after containment

**Neither endpoint is a clean skill measure and both are biased in known directions.**
The containment score is low largely because the model ignites from a *point* at
07:40 while the real fire was already **50 acres at discovery 12.6 h earlier** — that
is a missing head start, not slow spread. The observation-time score is high because
the model grew for 17 h after the real fire was held.

**Per JH, IoU around 0.25 or better counts as a success.** That is the calibration
anchor this work did not previously have, and it reframes earlier numbers: Tabor's
0.29-0.36 and Dry River's 0.330/0.341 (09-18 §7) are successes, not the mediocre
results they were written up as.

### The result that is not confounded

    fire                observed    modelled
    SINLAHEKIN (09-11)    2.54      1.38-1.57
    Dry River  (09-18)    2.90      1.73
    Tabor      (here)     1.93      1.01-1.26

**The shape deficit holds, and this is the first time it has been shown in the
standalone GRIB-driven path.** ForeFire produces a fire that is too round whether the
winds come from coupled WRF or from HRRR through WindNinja.

Elongation is also the metric least corrupted by suppression here. Suppression caps
how far a fire gets; there is no obvious reason it should make the real fire
*rounder* than the model, and flank-first attack — which the behaviour fields record
— would tend to make it narrower, widening the gap rather than explaining it.

### Why the ignition time was so far off, and what to do about it

Per JH, **VIIRS spatial resolution is decent but its temporal resolution is not**, so
a fire seen only by VIIRS gets an ignition time pinned to an overpass rather than to
ignition. Tabor's 12.6 h offset is that, not an error in the job.

His suggested direction: for VIIRS-only fires, **estimate the state of the fire at
the time of overpass and run the model forward from that estimate**, rather than
igniting a point at the overpass time. That would remove the largest confound in the
containment-time score above, and it is a different initialisation problem from
anything `forefire_grib.py` currently does.

## 7. Cotton 2 — winds, mesh and moisture

`COTTON_2_2026-09-21_20_00_00_283A99DF-AFD8-4F6A-8CA4-F827DA0AFBCD`, Napa County,
California. Ignition 2026-09-21 21:26:54Z at 38.63331, -122.06932. 30 km at 25 m.
**`ZSF` std 184.4 m** — the mid-roughness point between Red Bank (31.5) and Silver
(365), and a second California case. Deferred a day on 09-22 because its ~164 deg
wind reversal began on the last hour the f03 cache reached; by 09-23 the cache
covered the whole window with 31 gap-free GRIBs.

Five runs, each differing from a neighbour in exactly one thing:

    forefire_md_const     WRF winds,  constant Md
    forefire              WRF winds,  per-step Md
    forefire_hrrr_m25     HRRR winds, per-step Md, 25 m solve
    forefire_hrrr_m50     HRRR winds, per-step Md, 50 m solve
    forefire_hrrr_m100    HRRR winds, per-step Md, 100 m solve

**Reference:** WRF-SFIRE's own fire, from wrfout `FIRE_AREA` (verified cumulative —
monotonic, max 1.0, nonzero-cell count ~ sum, while `FUEL_FRAC_BURNT` is the
per-timestep rate): **817.9 ha**, stalling 09:00-18:00Z at 0.7-4.9 ha/h as its
`FMC_GC_F` climbs past 0.16.

### Mesh: it barely matters, and that is the operationally useful result

    solve mesh   ensemble mean   ratio to 25 m   build time
      25 m          4719.9 ha         -            ~65 min
      50 m          4656.0 ha       0.986          ~18 min
     100 m          4636.6 ha       0.982           ~5 min

**Under 2% across a 4x mesh range, for 13x the cost.** The prediction in §1 held:
fire area follows the bulk wind and ignores the field's detail.

Direction statistics are **even less sensitive** — mesh changes them by under 0.5%:

                     dir RMS mean   at reversal 10:00Z   >90 deg   peak B/A
      25 m               54.3            118.3            57.5%      2.97
      50 m               54.2            118.1            57.4%      2.59
     100 m               54.1            117.8            57.4%      2.26

**The only thing the mesh changes is the extremes**: peak B/A falls ~3.0 -> ~2.3 from
25 m to 100 m, which is 09-18 §4's "peak speed keeps climbing as the mesh refines"
seen from the other end. Peaks move 30%, area moves 2%. **Mesh buys extremes and the
fire ignores extremes.**

### The wind reversal is not resolved, at any scale

At 10:00Z, the onset, still coherent at R = 0.53 so the numbers mean something:
**direction RMS 118 deg with 57% of cells more than 90 deg wrong**, and the bias
swings +69 deg then -63 deg in consecutive hours. Through 11:00-14:00 R falls to
0.18-0.28, below the 0.3 floor, so those hours cannot be judged at all.

**Refining the mesh does not help.** If capturing a synoptic wind shift matters, this
is a limit of the method, not of the resolution it is run at.

### Winds

`HRRR/WRF = 1.082` at 25 m, both moisture-on. Cotton 2 sits between the flat cases
(~1.0) and the complex ones (Silver 0.72, Dome 1.16), consistent with §1.

### Moisture: real, and smaller than it first appeared

Same WRF winds, constant against per-step Md:

    chained track      MdConst 5325.9 ha, 0 of 53 hours stalled
                       MdVary  1936.2 ha, 8 of 53 hours stalled
    ensemble mean      MdConst 5587.5 ha  ->  MdVary 4362.4 ha   (0.781)
    ensemble spread    1.40 -> 6.15

**The stall is reproduced and it is partial, as predicted from per-fuel `me`.**
MdVary decelerates from 09:00, bottoms at ~2.2 ha/h around 15:00 and recovers by
16:30 — it never stops. Cotton 2 is 56% cat 2 (`me` 0.15) and 25% cat 5 (`me` 0.20),
so the majority arrests while a quarter creeps. A single global `me` would have
predicted a hard stop; WRF-SFIRE's own residual 0.7-1.7 ha/h shows the same shape.

**The track and the ensemble mean differ by 2.3x on the same run**, because
later-launched members start after the moisture peak and never stall — which is also
why spread explodes from 1.40 to 6.15. **For an 8 h product launched into a rising
moisture window the member that matters is the one starting at forecast time, not the
mean.**

**A correction.** An earlier reading of this said moisture "moves area by 550%" and
dwarfed the wind effects. **That was wrong** — it confused the size of the gap with
the size of the fix. Per-step moisture closes the gap to WRF-SFIRE from 6.8x to 5.3x
on the ensemble mean. It is clearly right, and it reproduces a mechanism ForeFire
structurally lacked, but it does not explain the overprediction.

### The overprediction is the open problem

    WRF-SFIRE        817.9 ha
    ForeFire         4362 ha ensemble mean / 1936 ha track, with moisture on

**A 5.3x overprediction remains after winds (8%), mesh (2%) and moisture (22%) are
accounted for.** Per JH this tendency for ForeFire to overpredict relative to
WRF-SFIRE has been noticed repeatedly, so it is a pattern rather than a Cotton 2
quirk, and it is now the largest known discrepancy in this work.

**JH's suggested direction:** determine which observables make the use of **wind
adjustment factors** sensible. That is a different lever from anything tried so far —
`windrf` and ForeFire's own `windReductionFactor` are currently taken from the
namelist and the template without being tied to anything measured. 09-18 §5 warns
that `UF`/`VF` already carry `windrf`, and 09-18 §7 that a `wRF` calibrated against
WindNinja means something different from one calibrated against coupled winds, so the
two must not be mixed.

### A hazard that cost three mistakes this session

**`cleanup_ff_run` empties the run directory it is given** (`mv {run_dir}/* {dest}/.`).
That consumed:

- the first 25 m ensemble, when the 50 m variant wrote into the same `forefire_hrrr`
  that still held it — identical filenames, silently overwritten;
- `timing_table_grib.csv`, which a later WRF rerun needed and could not find.

Pass a directory you are willing to lose, never write two runs into one destination,
and rename results immediately rather than leaving a directory occupied for the next
stage to walk into.

**Also:** moisture was enabled *mid-sweep*, so the original 25 m run was constant-Md
while 50 m and 100 m were not. The 25 m was rebuilt with moisture on rather than
comparing across a changed configuration. Do not change configuration inside an
experiment.

## 8. Silver-specific caveats that do not affect the ratio

- **`fire_init` was not invoked for this fire**, confirmed by measurement: cells
  inside the 2026-08-27 perimeter are 9.6% `NFUEL_CAT == 14` against 13.3% outside.
  A burn mask would have made inside ~100%. Per JH the module can acquire IR
  perimeters and mask consumed fuel to unburnable, but it did not run here, so both
  models burn through the August scar. The domain is sub-alpine with considerable
  talus near treeline, which is why there is a natural cat-14 background at all.
- **Absolute areas are therefore not realistic for this fire.** The ratio is,
  because both models use the identical fuel map.
- **The only perimeter available predates the forecast by 24 days** (§9), so nothing
  here is scored against observation.

## 9. A perimeter can be much older than the forecast

`ngfs/perims/Silver_{68A6E77D-...}.geojson` has `poly_PolygonDateTime`
**2026/08/27 18:56** against a forecast window of 2026-09-20 18:00Z to
2026-09-21 21:00Z. The fire was discovered 2026-08-23.

Per JH this is common rather than exceptional: **a slow-burning fire with a low heat
signature is invisible to GOES, so an NGFS detection fires when it flares up, not
when it started, and forecasts are regularly made for fires already weeks old.**
Silver was picked up by VIIRS; `ngfs_start2.py` in the `new_wrfxpy` testing install
created a job for it on 2026-08-23 (`ngfs/cron_ngfs_2026_08_23T143304.log`).

**Two consequences.** For scoring: check `attr_FireDiscoveryDateTime` and
`poly_PolygonDateTime` against the forecast window *before* choosing a validation
fire, not after. For fuels: a re-detected fire is modelled spreading through cells
already consumed unless `fire_init` runs, which inflates absolute area for both
models.

What the perimeter does establish is location, and it checks out: geometry is
self-consistent (709.1 ha computed against 714.7 ha stated, 99.2%) and the forecast
ignition sits inside it, 108 m from the nearest vertex.

## 10. The two fuel paths differ cell by cell

Measured on Silver, and it applies to **every** comparison including Union and Red
Bank, where it was never checked:

    NFUEL_CAT, WRF geogrid against forefire_grib's gdalwarp path
      exact category match     64.5%
      burnable / not-burnable  95.4%
      burnable fraction        WRF 86.8%   mine 86.8%
      category histograms      agree to <0.2% on every major category

Same distribution, a third of the cells shuffled. geogrid uses
`nearest_neighbor+average_16pt+search` with `dominant_only`; `forefire_grib` uses
gdalwarp `near`. They make different sub-pixel choices mapping a 30 m source onto a
25 m grid. Union and Red Bank agreed to 1-3.5% on area while carrying the same
difference, which is good evidence that spread responds to the aggregate
distribution rather than to cell-level assignment — but it is a confound that should
have been measured on the first fire, not the third.

Terrain, by contrast, agrees closely: the LANDFIRE DEM cut reproduced WRF's `ZSF`
range to 637-2299 m against 637-2300 m. That is the first check of the `altitude`
field in real relief, which 09-21 §10 flagged as untested on flat fires.

## 11. Open items

1. **Nothing is pushed**, on either branch, now across two sessions.
2. **§1 is two complex cases**, sitting either side of 1.0. The *sign* of the area
   error follows `B/A`; the coefficient is fitted to one point. Roughness is
   falsified as a predictor.
3. **`B/A` and direction statistics are missing for Union and Red Bank** — cheap to
   add now that the tool exists, and they would make the table meaningful.
4. **The overnight divergence is unexplained** (§3, §5) and `--diurnal_winds` is
   ruled out as the cause. It now appears on both complex fires, offset but the
   same shape, so it is a property of the method rather than of one case.
5. **Whether `B/A` predicts the area ratio quantitatively** — the sign is
   established on two cases, the coefficient is not. If it holds, a speed
   correction becomes possible; note 09-18's warning that a `wRF` calibrated
   against WindNinja means something different from one calibrated against coupled
   winds, so do not mix them.
6. Archived runs still carry wrong `valid_at` and mostly-empty ensembles from the
   09-21 §3 and §4 bugs. Union, Red Bank, Silver and Dome are redone; the rest are
   not, and how far back to go is still undecided.
7. ~~`dry_ff_compare.py`~~ — **done**, as `src/ingest/ff_score_perim.py` (§2, §6).
   09-18 §10 item 6 is fully closed.
8. **Ensemble spread does not generalise** (§5). Three fires give three orderings:
   HRRR more dispersed on Red Bank and Dome, tighter on Silver. Do not read a
   spread difference as meaningful without more cases.
9. **ForeFire overpredicts against WRF-SFIRE by ~5.3x on Cotton 2**, after winds,
   mesh and moisture are accounted for (§7). Per JH this tendency has been noticed
   repeatedly. **This is now the largest known discrepancy in this work.** His
   suggested direction: work out which observables make **wind adjustment factors**
   sensible, which is a lever nothing so far has touched.
10. **`cleanup_ff_run` empties the run directory it is given** (§7) — it destroyed an
    ensemble and a timing table this session. Never point two runs at one
    destination, and rename results immediately.
11. **Per-step moisture is on by default now** (`etc/forefire.json`, untracked).
    CONUS only, deliberately; off-grid fires keep the default table. The cron moved
    from every 20 minutes to `10 0 * * *` after interfering twice.
12. **Md is a domain-mean scalar per step**, not a field. The real moisture varies
    across a 30 km domain and ForeFire is being handed one number for all of it.

## 12. NEXT SESSION

1. **Finish Cotton 2** (§11 item 9). The cache should now cover its wind reversal,
   14 netcdfs of the pre-shift window are already built, and it is both a
   mid-roughness point for §1 and the only case so far with a large wind shift in it.
2. **Fill in §1 cheaply.** Run `wn_vs_wrf.py` on more fires spanning roughness — one
   netcdf and one wrfout per fire, minutes rather than the hours an ensemble costs.
   **Union and Red Bank first**, since their runs already exist and two flat points
   with wind statistics would test the weakest part of the `B/A` claim: whether
   area/wind stays near 1 when `B/A` is near 1.
3. **Then run ensembles only where the winds say it is interesting**, i.e. where
   `B/A` departs from 1. Both complex cases so far were worth it; a case with
   `B/A` ~ 1 would mostly confirm what the winds already said.
4. **Get a third point off the fitted line.** The `B/A` coefficient is currently
   fitted to Silver and checked on Dome. A case whose area ratio is predicted from
   *both* of them, in advance, is what turns this from consistent into established.
5. **Pick at least one complex-terrain fire where `fire_init` did run**, so absolute
   areas mean something and a perimeter score is possible — and check
   `poly_PolygonDateTime` against the forecast window first (§9).
