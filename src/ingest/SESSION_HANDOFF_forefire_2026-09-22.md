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

    fire        ZSF std*  wind B/A  dir RMS   anomR   area HRRR/WRF  area/wind  centroid
    -- scored against observation (§9) --
    Dome        obs_x 6.00  IoU 0.165   elong 1.16 vs observed 2.77
    Red Bank    obs_x 4.89  IoU 0.204   elong 1.76 vs observed 1.83   (timestamp +33 h, wide bar)
    Tabor       obs_x 1.6-1.9 IoU 0.29-0.36 elong 1.01-1.26 vs 1.93   (suppressed, ignited 12.6 h late)
    Silver / Cotton 2 / Union: no usable observation
    Dry River     26.7 m    1.04       3.7     0.129   (not run)         -          -
    Union         ~20 m*    -          -       -        0.99 +           -          -
    Red Bank      31.5 m    -          -       -        1.035            -          -
    Silver       364.7 m    0.765     49.1     0.109    0.717          0.937      514 m
    Dome         504.1 m    1.175     35.9     0.396    1.164          0.991      455 m
    Cotton 2     184.4 m    1.119     54.3     0.20*    1.082          0.967        -

    * DOMAIN-WIDE std, which overstates what the fire experiences by 14-51% and not
      uniformly -- see §8.  Near-fire (<=5 km) values: Red Bank 21.5, Cotton 2 162.3,
      Silver 241.0, Dome 377.9.  Ordering is preserved; use near-fire in new rows.
      Union: relief 21-93 m; std not computed.  Cotton 2 anomR is a rough mean over
      hours with R >= 0.3; its wind reversal makes several hours unusable.
      Silver, Dome and Red Bank area ratios were measured PRE-MOISTURE.  Dome's
      moisture-on 25 m ensemble is 5139.5 ha against 3922.8 pre-moisture (§8), so its
      area ratio would move; it has not been recomputed against a moisture-on WRF arm.
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

## 8. Dome revisited — mesh and moisture, and a correction to the roughness metric

Dome re-run 2026-09-23 with **per-step moisture** at three WindNinja solve meshes
(25/50/100 m, fire grid fixed at 25 m). The existing 25 m netcdfs were reused —
**moisture lives in the fuel table, not the netcdf** — so only 50 m and 100 m needed
building. The pre-moisture 25 m run is preserved as `forefire_hrrr_premoisture`.

### Mesh: insensitivity holds, and tightens

    25 m  5139.5 ha     50 m  5175.6 ha     100 m  5177.3 ha
    50/25 = 1.007       100/25 = 1.007      (Cotton 2: 0.986, 0.982)

**Under 1% across a 4x mesh range**, tighter than Cotton 2's 2%. Two fires now agree
that the fire ignores the solve mesh. A 100 m solve at ~5 minutes remains the right
operational choice against ~65 minutes at 25 m.

### Moisture *increased* area by 31%, and the reason matters

    25 m, constant Md   3922.8 ha
    25 m, per-step Md   5139.5 ha     ratio 1.310
    (Cotton 2: 0.781, with a 13-hour stall)

**Opposite in sign to Cotton 2, and a prediction of mine was wrong.** I expected
moisture to damp Dome modestly. It did the reverse, because:

    Dome      Md 0.055-0.106  mean 0.079   80% of steps BELOW the table default 0.1
    Cotton 2  Md 0.064-0.202  mean 0.118   57% of steps ABOVE it

**The constant `Md = 0.1` was never neutral.** It was wetter than reality at Dome and
drier than reality at Cotton 2's peak, so switching to real moisture speeds one fire
up and slows the other down. **Per-step moisture does not systematically shrink
fires; it removes an arbitrary constant.** That is a better argument for the feature
than "it makes fires smaller", and it means ForeFire's overprediction (§7) cannot be
blamed on moisture in either direction.

The no-stall prediction *did* hold: Dome is 65% timber litter (cats 8-10, `me`
0.25-0.30) and Md peaks at 0.106, nowhere near extinction.

**Per JH: below about 5-6% dead fuel moisture a fire is a good bet for big growth,
and there is a tipping point below which fires become largely unpredictable.** Dome's
minimum is **5.53%**, right at that threshold — so the 31% increase is not a modelling
curiosity, it is the model being handed conditions that really do drive large growth.

**Caveat carried from before the run:** the SAV weighting makes Md 93.5% a 1-hour
quantity, which suits Cotton 2's grass and fits Dome's timber litter poorly, where 10h
and 100h classes carry more of the behaviour. Part of the 31% is that configuration
choice and this run cannot separate it from physics.

### Correction: the case table's roughness column is misleading

**Per JH the Dome fire sits southwest of Yosemite Valley and west of the main Sierra
ridges** — the domain contains some of the steepest terrain on the planet, but the
fire does not burn in it. Measured:

    Dome ZSF std:  whole domain 504.1 m
                   within 2 km of ignition  209.7 m
                   within 5 km              377.9 m
                   5-10 km                  452.3 m
                   >15 km                   633.5 m   <- the Valley and main ridges

**So §1's roughness column, which is domain-wide, overstates what the fire
experiences.** Across the cases:

    fire        domain std   <=5 km std   ratio
    Red Bank        31.5        21.5       1.47
    Cotton 2       184.4       162.3       1.14
    Silver         364.7       241.0       1.51
    Dome           504.1       377.9       1.33

The **ordering is preserved** and Dome is still roughest near the fire, so the
conclusions stand — but the overstatement is 14-51% and **not uniform**, so the
domain figure should not be used to compare cases. An earlier claim that Dome tested
mesh insensitivity at "nearly 3x Cotton 2's terrain forcing" is really 2.3x
(378 against 162). **Use near-fire roughness in future rows.**

### The perimeter is an early-state observation

`ngfs/perims/DOME_{AACFD673-...}.geojson`: `poly_PolygonDateTime` **2026-09-15
16:09**, 41 minutes after discovery and **two hours *before* the job's 18:11:53
ignition**. 100.07 acres = 40.5 ha; geometry exact (computed 40.5, 100.1% of stated).

It cannot score the forecast's end, but it pins the head-start deficit precisely:

    real fire     40.5 ha at 16:09
    model         ignites a point at 18:11, reaches 40.5 ha at 20:00
    -> the model is about 3.8 h behind the real fire at equal size

Centroid 608 m from the modelled ignition, so location is good.

**Shape compared at matched size rather than matched time**, which sidesteps the lag:
at 20:00 the modelled fire is 60.3 ha with **elongation 1.14 against the observed
2.23**, IoU 0.306. The shape deficit again, in the same direction and magnitude as
SINLAHEKIN, Dry River and Tabor — and IoU 0.306 is a success by the 0.25 standard.

**Matching on size rather than time is worth reusing** wherever detection lag makes
the clocks incomparable.

## 9. Scoring against observations — the ytd collections, and why timestamps are suspect

### The per-fire files are often the worse observation

`ngfs/perims` holds per-fire files **and** `perims_ytd_<date>.geojson`, each carrying
every fire's latest perimeter. The per-fire file is frequently unusable while the ytd
series has something better for the same fire:

    Dome, per-fire file   40.5 ha stamped 2 h BEFORE the job's ignition -> unscoreable
    Dome, ytd series      476.6 ha at 2026-09-16 13:31, INSIDE the forecast window

Working from the per-fire file cost a real validation. **Run
`ff_score_perim.py --irwin <id> --list` before concluding a fire cannot be scored.**

`extract_from_ytd`/`ytd_series` and a `--irwin`/`--ytd-dir` CLI were added
(`41ab9dc`). Two traps, both of which broke a first attempt: the files are a single
**~170 MB line** and `poly_IRWINID` is `"{AACFD673-...}"` — **literal braces inside a
string** — so naive brace matching starts in the wrong place and never closes; and
these files date in **RFC-822** (`"Wed, 16 Sep 2026 13:31:00 GMT"`) where the per-fire
files use `%Y/%m/%d`.

### The rescore

    fire       observed            modelled    obs_x    IoU    obs elong  mod elong
    Dome        476.6 ha @13:31Z    2861.6      6.00    0.165     2.77       1.16
    Red Bank     83.9 ha @00:58Z     410.3      4.89    0.204     1.83       1.76
    Tabor       158.9 ha            255-302   1.60-1.90 0.29-0.36 1.93    1.01-1.26

Not scoreable: **Silver** (only perimeter is 2026-09-05, fifteen days pre-forecast —
the weeks-old-fire pattern), **Cotton 2** and **Union** (no perimeter in the
collections at all).

**The overprediction is now confirmed against ground truth, not only against
WRF-SFIRE.** Dome 6.0x and Red Bank 4.9x, both **below the 0.25 IoU standard** — the
first genuine failures in this work. Only Tabor is near 1, and that is the suppressed
fire where the model also ignited 12.6 h late, so the two errors partly cancel.

### The shape deficit looks regime-dependent

**Red Bank: observed elongation 1.83, modelled 1.76 — essentially right.** Every other
case had observed 1.93-2.90 against modelled 1.01-1.76. So ForeFire is not producing a
fixed roundness; it tracks when the real fire *is* round and fails when the fire is
elongated. Red Bank is the flat plains case with the roundest observation.

That reframes the deficit from "ForeFire makes round fires" to **"ForeFire cannot
produce strong elongation"** — a different and more diagnosable problem, and
consistent with 09-18's finding that the wind source does not fix it.

### IR observation timestamps are often untrustworthy

**Per JH, fire progress datasets are problematic and the timestamps on IR observations
are frequently not to be trusted.** `src/ngfs/perim_cache.py` runs on cron and its
`load_geojson_with_processed_time` attaches a **`processed_utc` taken from the file's
mtime** — when the change was *noticed*. That is not stored in the ytd files
themselves (they carry none), so in practice the bound is **the ytd file's own
mtime**.

Using it as a sanity check:

    perims_ytd_2026-09-16.geojson  mtime 09-16 16:00  carries Dome's 13:31   (+2.5 h)
    perims_ytd_2026-09-17.geojson  mtime 09-17 10:00  carries Red Bank 00:58 (+33 h)

**Dome's claimed time is corroborated; Red Bank's is not.** A 33-hour gap means the
flight could be anywhere in that window, and since the modelled fire grows
monotonically, `obs_x` moves a long way with it. **Treat Red Bank's 4.89x as having a
wide error bar**, and check this bound before trusting any score.

Note also that the ytd **filename date is UTC while the mtime is local**, so a file
named `2026-09-20` can have an mtime of `2026-09-19 19:00`. Do not infer the notice
time from the filename.

## 10. Fuel moisture as a field — what is reachable and what is not

§7 item 12 (now §12 item 12) records that Md is a scalar. **Per JH, moisture as a field could be a
ForeFire input if a different spread model is used** — several Balbi variants exist —
though that would need a new fuels csv and possibly different fuel maps. Checked, and
the position is better than expected in one way and blocked in another.

### The models are all there

`libforefireL.so` inside the container carries **all eleven** propagation models:
`Rothermel`, `RothermelAndrews2018`, `Balbi2015`, `Balbi2020`, `BalbiNov2011`,
`BalbiNov2011Curv`, `BalbiNov2011TMdMl`, `BalbiUnsteady`, `Farsite`, `IsotropicFuel`,
`WindDriven`. The `forefire` executable is only 68 KB — a thin driver — so look in the
library, not the binary. Selecting one is a `setParameter[propagationModel=...]`
change in the templates, no rebuild.

### Moisture layers are already plumbed

`DataBroker.cpp` reads **`moisture` and `temperature` as XYZT layers straight out of
the netcdf**, by the same machinery as `windU`/`windV` (lines ~343, ~1285). So a
spatially *and* temporally varying field needs no new mechanism — `forefire_grib.py`
would just write another variable.

### But dead moisture is deliberately wired to the table

In `BalbiNov2011TMdMl.cpp`:

    deadMoisture = Md; //registerProperty("deadMoisture");

**Live moisture and temperature come from layers; dead moisture does not.** Someone
tied it back to the fuel table and commented out the layer registration. So today:

    live moisture, temperature   can be fields, no source change
    dead moisture                per-category table constant

### And the source change is effectively out of reach

**Per JH: the Singularity container was built on a separate machine from a Docker
image published by the ForeFire developers. Compiling ForeFire from source needs
specific compilers and a C++ netcdf library he could not get working in any
environment here.**

So uncommenting that one line is *not* a one-line change in practice — it requires a
rebuild that has not been achievable. Treat the container as fixed. The routes that
remain are to ask upstream, or to solve the build environment, and neither is a
afternoon's work.

### What a model switch would cost anyway

- **The fuel table.** `ff_fuels_behave13.csv` already carries Balbi parameters
  (`Rhod`, `Tau0`, `Deltah`, `r00`, `X0`, `Blai`), which is why one file serves a
  Rothermel run — but whether those values are *calibrated* for Balbi or merely
  inherited is unverified and must not be assumed.
- **The fuel maps.** Balbi variants separate dead and live loads (`Sigmad`/`Sigmal`)
  and depth, which Anderson 13 supplies only coarsely.
- **Every calibration measured here.** `windrf` and `windReductionFactor` were tuned
  against Rothermel behaviour, and 09-18 §7 warns such calibrations do not transfer.
  **A model switch would invalidate the `B/A` relationship measured across five fires
  rather than extending it.**

**So the practical near-term option is the variance-dependent radius of §7 item 12,
not a field.** A field is the better answer and it is blocked on a build problem, not
on a design one.

## 11. Hot Spring 226 — FMDA saw the rain and barely responded

`Hot_Spring_226_2026-09-21_16_00_00_F99FF1C2-AB69-44DD-A118-082DF3E7E09A`, Hot Spring
County, Arkansas. Ignition 2026-09-21 17:17:20Z at 34.4076, -92.6637. Flat (`ZSF` std
30.3 m), timber litter — cat 9 33%, cat 8 24%, cat 10 19%, `me` 0.25-0.30.

JH's question: the forecasts matched closely, but there was a jump in fuel moisture in
the window, **possibly precipitation FMDA was blind to and WRF was not.**

### FMDA was not blind to it

    rain event 22:00-00:00Z
      WRF   RAINNC 0.00 -> 0.23 mm accumulated    Md 0.0772 -> 0.1030   (+0.0258)
      FMDA  PRECIP 0.009 then 0.050               Md 0.0637 -> 0.0667   (+0.0030)

**FMDA's `PRECIP` registers the event in exactly the hours WRF does.** The difference
is the *response*: WRF's fuel moisture jumps **8.6x** more than FMDA's for the same
rain. Whether that is FMDA's Kalman update damping a small forcing against its RAWS
observations, or a genuinely different wetting response, is not separable from this.

Over the whole 27 hours the two agree well — **correlation 0.770, bias +0.0021, RMS
difference 0.0236**, FMDA peaking slightly higher (0.160 against 0.146). A localised
three-hour discrepancy inside a well-tracked series, not systematic blindness.

### Why the forecasts agreed anyway

Neither source ever exceeds 0.16 against `me` of 0.25-0.30, so **moisture never
approaches extinction here and cannot drive a stall in either model.** Same regime as
Dome; the opposite of Cotton 2, where grass at `me` 0.15 was crossed for 13 hours.

### Re-run with per-step moisture

    growth rate        ForeFire (Md vary)   WRF-SFIRE
      22:00                 85.4 ha/h          11.8
      23:00                 93.5               23.0
      00:00                 72.5  <- dip       17.3  <- dip
      01:00                 80.0               26.7

**Both dip at 00:00Z and both recover** — ForeFire to 78% of its prior rate,
WRF-SFIRE to 75%. Per-step moisture *does* reproduce the feature at comparable
relative depth, which **falsifies** a prediction that the 8.6x weaker wetting response
would give a much shallower dip.

**The absolute overprediction is untouched:** ForeFire 1551.8 ha over 24.5 h against
WRF-SFIRE's 440.9 ha, **3.5x** — consistent with Dome (6.0x against observation), Red
Bank (4.9x) and Cotton 2 (5.3x).

### The premise rested on a truncated run

**"ForeFire matched WRF-SFIRE closely" came from a run that stopped early.** The old
constant-Md result covered 17:30-05:00 and reached 630.6 ha, which against 440.9 ha
looked reasonable. The complete run reaches 1551.8 ha. The agreement was the run
stopping, not the models agreeing. The preserved `forefire_md_const` is **not** a
clean const-vs-vary comparison for the same reason — 24 perimeters against 50.

## 12. Most recent runs were truncated, and why

Surveying the last 36 workspaces with ForeFire output: **21 are truncated**, last
perimeter earlier than last wrfout.

    Calhoun_259    last wrfout 09-23 09:00    last FF perimeter 09-22 21:00
    BOON           last wrfout 09-23 00:00    last FF perimeter 09-22 12:00

**The cause is scheduling, not failure.** The cron fires at 00:10 while WRF is still
writing wrfouts, so `make_timing_table` globs a partial set and the run stops there.
Re-running once WRF has finished picks up the full set.

**Caution on the survey method:** runs made *before* the 09-21 clock fix carry
`valid_at` stamps that run ahead of real time, so they appear "complete" against the
last wrfout when they are not. The test is only meaningful for post-fix runs —
identifiable by having many distinct `fuelsTableFile` entries rather than one.

**Per JH, all recent forecasts are worth re-running with fuel moisture where
possible.** `ff.run_days(days2run=8, overwrite=True)` does this and was launched.

**The nightly cron will keep producing truncated runs** unless it is scheduled after
WRF reliably finishes, or made to skip fires whose WRF is still going. Nothing does
either yet.

## 13. Recovering WRF-SFIRE's fire from TIGN_G, and the Dome triangulation

### Cleaned workspaces are not lost

Per JH, a cleaned workspace has its wrfouts replaced by `saveout_d01_*`, and **fire
progression survives in `TIGN_G` in the retained final wrfout**. TIGN_G is each
fire-grid cell's ignition time in seconds since the simulation start, with a large
sentinel for cells that never burned, so **area at time t is the count of cells at or
below t** — the whole progression from one file.

`ff_growth.wrf_track` does this, and the CLI takes `wrf:<label>=<workspace>` so
WRF-SFIRE's own fire appears as a column beside the ForeFire forecasts.

**Validated where both variables survive.** On Hot Spring 226, TIGN_G gives 442.8 ha
against `FIRE_AREA`'s 440.9 at the final time and tracks at every hour. It runs
**0.5-1% high** because TIGN_G marks a cell burned at its ignition instant while
FIRE_AREA ramps fractionally — a consistent bias, not noise.

**It is better than FIRE_AREA for this purpose:** one file instead of 61, it works on
cleaned and uncleaned workspaces alike, and it evaluates at any instant rather than
only at output times.

**An earlier note in `wrfout_status` claimed cleaning made WRF-SFIRE's growth curve
unrecoverable and model-against-model comparison impossible. That was wrong** and is
corrected. The WRF-SFIRE reference does **not** have a shelf life; the archive stays
usable.

### Dome: observation, WRF-SFIRE and ForeFire together

Dome's `FIRE_AREA` was already gone, so this comparison existed only because TIGN_G
recovered it. At the observation time **2026-09-16 13:31Z**:

    observed (NIFC IR, ytd)     476.6 ha    1.00x
    WRF-SFIRE (from TIGN_G)     322.3 ha    0.68x
    ForeFire (HRRR, 25 m)      2861.6 ha    6.00x

    ForeFire / WRF-SFIRE = 8.88x
    final areas: WRF-SFIRE 613.6 ha, ForeFire ensemble mean 5139.5 -> 8.4x

**The two models bracket the observation.** WRF-SFIRE underpredicts by a third;
ForeFire overpredicts by six. That rules out the possibility that ForeFire only looks
bad because WRF-SFIRE is a soft target — **WRF-SFIRE is under the truth and ForeFire
is far over it**, so ForeFire's error against reality (6.0x) is real and its error
against WRF-SFIRE (8.9x) overstates it only slightly.

Hour by hour the divergence is stark: at 14:00Z ForeFire grows at **247 ha/h against
WRF-SFIRE's 29**, widening to 347 against 42 by 20:00Z.

This is the strongest evidence yet for the overprediction being ForeFire's, and it
sharpens §7 — winds, mesh and moisture together account for tens of percent
against a discrepancy approaching an order of magnitude.

## 14. IN FLIGHT AT HANDOFF — the batch re-run

`ff.run_days(days2run=8, overwrite=True)` was launched to re-run every recent
workspace with per-step moisture, fixing the truncation of §12.

    process   detached (PPID 1), so it should survive the session ending
    log       <session scratchpad>/rerun_all.log
    progress  6 fires when this was written

**The log lives in the session scratchpad and may not persist.** The durable way to
check is the workspaces themselves — re-run the survey in §12: a completed
re-run has its last ForeFire perimeter at or after its last wrfout, and many distinct
`fuelsTableFile` entries rather than one.

**Two caveats for reading the results.** The batch started *before* the
completeness guard was committed, so the already-running process is using the old
code and may still produce truncated output for fires whose WRF was unfinished; the
guard takes effect on the next invocation, including the nightly cron. And it started
before the saveout fix, so it may have skipped nothing but will have wasted time on
cleaned fires it could not previously recognise.

## 15. Silver-specific caveats that do not affect the ratio

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

## 16. A perimeter can be much older than the forecast

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

## 17. The two fuel paths differ cell by cell

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

## 18. Open items

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
12. **Md is a scalar per step, not a field, and the averaging radius is a real
    choice.** An earlier draft of this item called it a "domain-mean" — it is not.
    `read_fmda_hourly` averages within `moisture.radius_km`, default **10 km**, around
    the ignition on a 30 km domain. Measured on Cotton 2:

        FMDA grid spacing near the fire: 1.66 km  (so r=1 km catches no cells at all)

        06Z   r<= 5 km  mean 0.0960  sd 0.0107  range 0.076-0.113
              r<=10 km  mean 0.1133  sd 0.0179  range 0.076-0.148
              r<=21 km  mean 0.1226  sd 0.0173  range 0.076-0.191   straddles me
        12Z   r<= 5 km  mean 0.1762  sd 0.0421  range 0.132-0.256   straddles me
              r<=21 km  mean 0.1963  sd 0.0443  range 0.104-0.266   straddles me

    **The radius shifts Md by ~20%** (0.096 at 5 km against 0.123 at 21 km, 06Z).
    **The domain straddles cat 2's `me` of 0.15** — at 12Z even cells within 5 km span
    0.132-0.256, so part of the domain is arrested and part is burning freely while
    the model applies one number everywhere. And **the variance is itself
    time-varying**: sd 0.011 at 06Z against 0.042 at 12Z, a 4x change, so no fixed
    radius is right for both.

    **Per JH:** derive the mean from locations nearer the ignition, *depending on the
    variance across the domain* — fuel conditions away from the fire may differ from
    those where it is likely to spread. That argues for a variance-dependent radius
    rather than a tuned constant: widen while sd stays low, tighten toward the
    ignition when it does not. **The floor is the grid** — at 1.66 km spacing, ~5 km
    is the tightest radius averaging more than a handful of cells, and below that it
    reads one or two RTMA cells with whatever noise they carry.

    The deeper fix is Md as a **field** — see §8, which finds the layer machinery
    already exists but that dead moisture is wired to the table with the layer
    registration commented out, and that rebuilding the container to change it is not
    currently achievable.

## 19. NEXT SESSION

0. **Check the batch re-run first** (§14). It was still going at handoff and is
   detached, so it should have continued. Re-run the survey in §12 to see how many
   workspaces are now complete, and note that the batch predates both the
   completeness guard and the saveout fix.

1. (see §12 item 9). The cache should now cover its wind reversal,
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
