# ForeFire session handoff — 2026-09-24

Continues `SESSION_HANDOFF_forefire_2026-09-22.md`. Three strands, in the order they
were worked: calibrating `propagationSpeedAdjustmentFactor` across the batch, giving
each fuel category its own dead fuel moisture, and fixing a guard that was blocking
fires permanently. The pSAF thread produced the day's main correction and it is the
one to read first.

## 1. TL;DR

- **pSAF 0.175 overshot.** Measured over 59 fires: geometric mean **0.58**, 90% below
  1.0x. ForeFire now underpredicts.
- **The exponent was right; the baseline was wrong.** The 6.74x overprediction that set
  0.175 came from comparing **final areas at different valid times**. Scored at a
  common valid time the pSAF-0.6 overprediction was **3.85x**, not 6.74x. The widely
  quoted "ForeFire runs about 8x too much area" is substantially a timing artifact.
- **Now running at pSAF 0.337**, which puts the 25th percentile at 1.0 — about three
  quarters of fires overpredict, per JH's stated preference.
- **Md is now per fuel category**, weighted by WRF-SFIRE's own `fmc_gw01..03`. Effect on
  total area is small (+3.8% / -3.4%) but it is the physically right treatment.
- **The completeness guard has an age escape.** Two fires had been blocked for four days.

## 2. pSAF calibration — the main result and the main correction

### What was applied

Both templates in `/home/jhaley/forefire/tests` (`ignition_template.ff`,
`restart_template.ff`) went from `propagationSpeedAdjustmentFactor=0.6` to `0.175`, and
`run_days(days2run=10, overwrite=True)` re-ran the whole batch — 62 workspaces, ~2 h.

### What came back

Scored at the **last valid time both models reach**, ForeFire against WRF-SFIRE area:

| | n | geo mean | median | below 1.0x | within 2x |
|---|---|---|---|---|---|
| as first reported | 59 | 0.577 | 0.58 | 90% | 53% |
| stale pair removed (see §4) | 57 | **0.531** | 0.56 | 93% | — |

### Why it overshot

**Not the exponent.** The two-fire ladder (Dome 1.65, Hot Spring 226 1.42, mean 1.54)
predicted both calibration fires almost exactly:

| fire | predicted | measured |
|------|-----------|----------|
| DOME | 0.90 | **0.89** |
| Hot Spring 226 | 0.61 | **0.62** |

**The population baseline.** 0.175 was derived from a 6.74x median overprediction
measured by comparing each run's **final** area. But ForeFire chains to the last GRIB
while WRF-SFIRE stops at its own last output, so final-vs-final quietly credits ForeFire
with extra hours of growth. Under common-valid-time scoring the implied pSAF-0.6
baseline is **3.85x**. Correcting for that is the whole of the overshoot.

**This invalidates the framing that motivated the exercise.** "ForeFire is consistently
producing about 8 times too much fire area" was measured the same way and carries the
same inflation. The honest number is closer to 3.9x, and a good part of what looked like
a rate-of-spread error was a clock mismatch.

### Where it goes next

From the clean set (n=57) with exponent 1.54:

| target | pSAF | % below 1.0x |
|--------|------|--------------|
| centre the median | 0.255 | 51% |
| **25th pct at 1.0** | **0.340** | **25%** |
| 10th pct at 1.0 | 0.382 | 9% |

JH chose the 25th-percentile posture. **The batch is running at 0.337** (computed before
the stale pair was removed; the clean figure is 0.340 and the difference is far below
the noise). Predicted outcome: geo mean 1.46, median 1.54, 26% below 1.0x.

### Scoring tool

`src/ingest/ff_batch_ratio.py` (added today). It walks
`wksp/wfc-*`, reads ForeFire's chained step perimeters and WRF-SFIRE's TIGN_G curve,
intersects their valid times and reports the ratio at the latest common one. Fires with
no overlap are listed, not dropped.

## 3. Md is now per fuel category, using WRF-SFIRE's own weights

Until now every Anderson category got the **same** Md. `read_fmda_hourly` collapsed
FMC_GC classes 0/1/2 into one scalar with fixed weights (2000, 109, 30) and
`write_fuel_table` wrote it into the `Md` column of every row, so short grass and heavy
slash were handed identical moisture. Per JH the split should come from how much of
each fuel is 1-h, 10-h and 100-h.

### Where the weights come from

JH pointed at `etc/nlists/default.fire_cawfe_13`, which carries `fmc_gw01..fmc_gw05` —
one array per moisture class (1-h, 10-h, 100-h, 1000-h, live), 14 entries each. This is
**WRF-SFIRE's own weighting**, which settles the question: every ForeFire result here is
scored against WRF-SFIRE, so both should read a moisture field the same way.

Two facts from that table shape the implementation. `fmc_gw04` (1000-h) is **all zeros**,
so dropping the 1000-h class costs nothing. `fmc_gw05` (live) is non-zero for fuels 2, 4,
5, 7 and 10 and FMDA has no live class to fill it, so the dead weights are
**renormalised** over the three classes we have. Live stays in the table's own `Ml`
column rather than being folded into Md.

FMDA's `FMC_GC` has 6 slots but only 0-3 are moisture: `fuel_moisture_model.py` writes
`FMC_GC[:,:,:-2] = m_ext[:,:,:-2]` and the trailing two are Kalman extended state
(measured negative everywhere). `fuel_moisture_da.py` passes only `[:,:,:3]` to WRF, so
1-h/10-h/100-h is what both models actually consume.

### CAWFE weighting is load weighting, not surface-area weighting

Worth recording because the first implementation got this wrong. Rothermel weights a
size class by **surface area** (load x SAV); the plainer reading is **load fraction**.
They are not close for heavy fuels. Comparing CAWFE against both reconstructions:

| fuel | CAWFE 1h/10h/100h | area (Rothermel) | load (Anderson) |
|------|-------------------|------------------|-----------------|
| 1 grass       | 1.000 0.000 0.000 | 1.000 0.000 0.000 | 1.000 0.000 0.000 |
| 9 hdwd litter | 0.066 0.930 0.003 | 0.993 0.006 0.001 | 0.839 0.118 0.043 |
| 10 timber     | 0.300 0.200 0.500 | 0.942 0.034 0.024 | 0.300 0.200 0.500 |
| 12 slash      | 0.116 0.406 0.478 | 0.748 0.190 0.062 | 0.116 0.406 0.478 |

**CAWFE matches load weighting exactly for 11 of 13 fuels.** Area weighting stays 1-h
dominated everywhere (never below 0.748) because 1-h SAV is 18x the 10-h and 67x the
100-h, and it turned out to be nearly a no-op: 0.34% of cells crossed `me` against
8.14% for CAWFE, a 24x difference.

**Two fuels disagree with Anderson and are table errors.** Fuel 2 is 0.250/0.125/0.625
where Anderson loads give 0.571/0.286/0.143, and fuel 9 is 0.066/**0.930**/0.003 where
Anderson gives 0.839/0.118/0.043 — 1-h and 10-h transposed. **Per JH both are corrected
to the Anderson values**: fuel 9 is hardwood litter and is mostly a 1-hour fuel, as the
Anderson weighting shows, and a 93% 10-h weight would make eastern hardwood fires nearly
unresponsive to diurnal drying. Fuel 2 is timber grass and understory, which 0.625 on
100-h does not describe either.

This is not a marginal correction. **Cotton 2 is 63.3% fuels 2+9 by area** (55.95% fuel 2,
7.36% fuel 9) and **Hot Spring 226 is 48.7%** (33.32% fuel 9, 15.37% fuel 2), so the two
overridden categories cover roughly half to two thirds of both domains. The divergence
from WRF-SFIRE is deliberate and should be raised upstream — `class_weight_anderson_fuels`
set to `[]` follows the namelist exactly.

The Anderson reconstruction was validated on the way through and is kept as a fallback:
surface-area weighting built from those loads reproduces the `sd` column of
`etc/ff_fuels_behave13.csv` to 4 significant figures for all 9 fuels with no live load,
and the 4 that deviate order monotonically by live/dead ratio (5 > 4 > 10 > 2 > 7).

### It matters most when conditions move fast, as JH predicted

1-h fuel equilibrates within the hour, 100-h lags by days, so the gap should open when
temperature and humidity are moving. Measured over 131 fire-hours on six fires, the
grass-to-timber Md spread is **1.4x wider** in the fastest-changing third of hours
(0.089) than the slowest third (0.062). The *instantaneous* correlation is weak (+0.06)
because 100-h carries days of history, so the gap persists well after the change that
opened it.

Hot Spring 226 shows it directly. Across the wetting event the 1-h class climbs
0.069 -> 0.224 (3.2x) while 100-h moves 0.146 -> 0.151 (1.03x):

| hr | grass (f1) | hdwd litter (f9) | timber (f10) | slash (f13) |
|----|-----------|------------------|--------------|-------------|
| 00 | 0.074 | 0.089 | 0.113 | 0.115 |
| 09 | **0.224** | 0.149 | 0.170 | 0.156 |
| 21 | 0.069 | 0.112 | 0.117 | 0.125 |

The **ordering inverts** — grass is the driest fuel at hour 0 and the wettest at hour 9 —
which a single scalar cannot represent at all. At hour 9 grass sits above its `me` of
0.12 and is shut down while timber and slash sit at 0.156-0.170 against `me` 0.25 and
keep burning; the old scalar handed all of them 0.219, and fuels 4, 5 and 12 flip from
blocked to spreading.

### Switches

`etc/forefire.json`, under `moisture`:

    "weighting": "per_fuel",            # was "sav"; "single" and "sav" still work
    "class_weighting": "cawfe",         # or "load", or "area"
    "class_weight_namelist": "/data/jhaley/wrfxpy/etc/nlists/default.fire_cawfe_13",
    "class_weight_anderson_fuels": [2, 9]   # [] to follow the namelist exactly

`read_fmda_hourly(..., per_class=True)` returns the class tuple instead of the collapsed
scalar; `write_fuel_table` accepts a float (old behaviour, every row) or a tuple (per
row). A fuel index outside 1-13 — the no-fuel row — falls back to the collapsed scalar,
and a class FMDA did not supply that hour falls back rather than dropping the step.

## 4. The completeness guard needed an age escape

The guard from §12 skips a fire whose wrfout count is short, on the grounds that WRF is
still writing and the forecast would be truncated. Correct while WRF is running, wrong
once it has stopped: **a run that ended early never gains a file, so the fire is blocked
for ever.** Union_400257 sat at 52 of 55 and Ouachita_259 at 54 of 55 for four days,
skipped by every nightly run — three hours short of a full forecast, and getting none.

This surfaced as a scoring artifact, not as a complaint. Both fires showed up as wild
outliers in the pSAF 0.175 table (4.73x and 7.77x against a set centred near 0.58x)
because the guard skipped them, their `forefire` dirs still held the previous night's
**pSAF 0.6** output, and a 24 h recency filter let 1438-minute-old files through. Nothing
about the fires needed explaining. Per JH the fuel moisture at both tracks closely and
rises only near the end of the window, which is consistent: moisture was never the cause.

`wrfout_status` now takes cfg and adds two knobs:

    wrfout_stale_hours   default 6    newest output unchanged this long -> WRF is
                                      finished-short, not running; proceed
    wrfout_min_fraction  default 0.5  floor: below this, WRF died early and a badly
                                      truncated forecast is worse than none

Measured on the three fires the old rule blocked:

| fire | count | idle | outcome |
|------|-------|------|---------|
| Union_400257 | 52/55 | 81 h | **RUN** — finished short, 3 missing |
| Ouachita_259 | 54/55 | 81 h | **RUN** — 1 missing |
| Castle_Peak_RX | 11/55 | 81 h | **SKIP** — died early, under 50% |

Complete runs and cleaned workspaces (saveout path) are unaffected. The shortfall is
stated in the message rather than hidden, so a truncated forecast is identifiable
afterwards.

**Scoring lesson worth keeping.** A recency filter for batch scoring must key off the
batch's own start time, not a guessed number of hours — `ff_batch_ratio.py` now uses the
batch script's mtime. Two stale entries moved the reported geometric mean from 0.531 to
0.577 and would have silently biased the pSAF recommendation.


## 5. What per-fuel moisture actually did to the two test fires

Cotton 2 and Hot Spring 226 were re-run with per-fuel Md at pSAF 0.175, against their
own pSAF-0.175 scalar-Md output preserved first as `forefire_psaf175_scalarmd`. Same
pSAF, same winds, same mesh — the moisture weighting is the only variable.

Compared at a common valid time (Cotton 2's rerun picked up 53 steps against the
baseline's 35 as more GRIBs arrived, so raw finals are **not** comparable):

| fire | scalar Md | per-fuel Md | change | ratio vs WRF |
|------|-----------|-------------|--------|--------------|
| Cotton 2 | 149.7 ha | 155.3 ha | **+3.8%** | 0.35 -> 0.36 |
| Hot Spring 226 | 274.2 ha | 264.7 ha | **-3.4%** | 0.62 -> 0.60 |

Small, and **opposite in sign**. The prediction going in was that Hot Spring would grow
*more*, because fuels 4, 5 and 12 drop below their `me` at the humid peak where the
scalar had them blocked. It shrank 3.4% instead, and the reason is worth keeping:

**Per-fuel weighting damps the diurnal amplitude rather than shifting the level.** Heavy
fuels get *wetter during dry hours* (100-h stays wet when 1-h crashes) and *drier during
humid spikes* (100-h does not spike). Which phase dominates a given run decides the sign
of the net area change.

The intended effect is real and visible hour by hour — through Hot Spring's humid window
the per-fuel run grows consistently faster, 22.2 vs 16.6 ha/h at 06:00Z and 16.7 vs 14.5
at 10:00Z — it simply does not survive into the total.

**Read this as a null result on area and a correctness win on mechanism.** Do not expect
per-fuel Md to move the ForeFire/WRF-SFIRE ratio; expect it to make growth respond to the
right fuel at the right hour.

## 6. Corrections made during the session

Recorded because each was asserted before it was checked.

1. **"ForeFire runs ~8x too much area."** Measured final-vs-final. Common-valid-time is
   3.85x. This one cost a whole batch re-run at the wrong pSAF.
2. **"The two big outliers are a WRF-side disagreement."** Backwards. WRF-SFIRE's 103-157
   ha sat squarely inside its regional population (196, 258, 443 ha); ForeFire's 732 and
   802 ha were 3-6x its own regional norm — and in the end neither was real (§4).
3. **"Rothermel surface-area weighting is the correct physics for Md."** True of
   Rothermel, wrong here: WRF-SFIRE uses **load** weighting, which is what JH described.
   Area weighting would have been near a no-op (0.34% of cells crossing `me` vs 8.14%).
4. **"Per-fuel Md should make Hot Spring grow more."** It shrank 3.4% (§5).
5. Two defects in same-session code, found by testing rather than by review: an
   unreadable namelist raised instead of falling back, and that fallback then landed on
   area weighting when it should land on load. `time` was also missing from the imports
   when the guard escape was added.

## 7. Code state

Modified, **uncommitted**, all in `src/ingest/forefire.py` unless noted:

    DEAD_LOAD_13, SAV_1H_13, SAV_10H/100H, CAWFE_NAMELIST, CAWFE_ANDERSON_FUELS
    dead_class_weights(fuel, mode, namelist_path, anderson_fuels)   new
    cawfe_class_weights(namelist_path, n_dead, anderson_fuels)      new
    read_fmda_hourly(..., per_class=False)                          extended
    write_fuel_table(md, cfg)          md may now be a class tuple
    md_series                          per-class path + per-class range check
    make_script_set                    cache key handles a tuple
    wrfout_status(wksp_dir, cfg=None)  age escape + min-fraction floor
    import time                        was missing

`etc/forefire.json` (untracked) gained:

    "weighting": "per_fuel",
    "class_weighting": "cawfe",
    "class_weight_namelist": ".../etc/nlists/default.fire_cawfe_13",
    "class_weight_anderson_fuels": [2, 9]

Backups: `/tmp/psaf_work/forefire.json.bak-scalarmd`,
`/tmp/psaf_work/{ignition,restart}_template.ff.bak-psaf0.6`.

**Added today:** `src/ingest/ff_batch_ratio.py` (§2). It is the only tool that scores
the batch at a common valid time, and every number in §2 depends on it.

23 commits sit unpushed on `james_ngfs`; JH pushes these manually (SSH password).

## 8. In flight at handoff

- **pSAF 0.337 batch** — `run_days(days2run=10, overwrite=True)`, 62 workspaces, ~2 h,
  log `/tmp/psaf_work/psaf337_batch.log`, ends with `PSAF337 BATCH COMPLETE`.
  This is the first batch carrying **both** pSAF 0.337 and per-fuel moisture.
- **Union_400257 and Ouachita_259** — queued behind it on the same `flock`, log
  `/tmp/psaf_work/unblocked.log`. The 0.337 batch imported `forefire.py` before the age
  escape existed, so it cannot pick them up; this run brings them onto the same footing.

**Early partial read — 7 fires in, treat as a smoke test, not a result.** Geometric mean
1.98, median 1.65 against a prediction of 1.46 / 1.54. The sample is almost all
prescribed burns and includes a 0.7 ha fire, so it is not representative. Per-fire the
exponent is holding loosely (Eagle 0.50 -> 1.24 against 1.37 predicted, Dalton_RX
1.10 -> 2.80 against 3.02) with one large miss (RX_Stray_Creek 1.10 -> 9.63 against
3.02). Extra scatter is expected because these runs also carry per-fuel moisture.

To score when both finish:

    cd /data/jhaley/wrfxpy && PYTHONPATH=src:src/ingest \
      python src/ingest/ff_batch_ratio.py /tmp/psaf_work/psaf337_batch.py

## 9. Open items

1. **Verify 0.337 against the prediction** (geo mean 1.46, median 1.54, 26% below 1.0x).
   A large miss means the exponent does not hold over a 1.9x pSAF move, which would be
   new information — it held exactly over the 0.6 -> 0.175 move.
2. **pSAF probably belongs per fire, not global.** The size dependence is still there:
   `ratio ~ area^-0.40`. JH's standing direction is to predict it from NGFS detection
   fields (`frp`, `feature_frp`, `feature_detection_duration`, `fuel`, `land_cover`) over
   a large CONUS batch. Exclude connectivity outliers like Mud_Bayou.
3. **Raise fuels 2 and 9 upstream.** `default.fire_cawfe_13` has fuel 9 at
   0.066/0.930/0.003 where Anderson gives 0.839/0.118/0.043 — 1-h and 10-h transposed.
   Per JH hardwood litter is mostly a 1-hour fuel. We now diverge from WRF-SFIRE here
   deliberately; WRF-SFIRE itself is still wrong and its own forecasts carry it.
4. **`moisture_params` prints "unknown moisture source 'fmda_hourly'" then "FMDA moisture
   unusable"** on every run. Harmless — `md_series` does the per-step work — but the
   message is false and will mislead someone reading a log.
5. **Score against observations, not only WRF-SFIRE.** WRF-SFIRE is not truth: on Dome it
   ran 0.68x an IR perimeter while ForeFire ran 6.0x. Calibrating pSAF model-to-model
   inherits WRF-SFIRE's own bias.
6. **Two fires still unscored** — Nethery and STEEL_PASS had no overlapping valid time.

## 10. NEXT SESSION

Start by scoring the 0.337 batch (§8). If it lands near the prediction, the global-pSAF
question is closed at "0.34, with a known size dependence" and the next real move is
item 2 — per-fire pSAF from detection data, which is the only thing that addresses the
`area^-0.40` slope rather than centring around it.
