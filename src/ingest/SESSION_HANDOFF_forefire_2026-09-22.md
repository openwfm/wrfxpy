# ForeFire session handoff — 2026-09-22

Branch `james_ngfs`. Commits this session: `af212f4` (the wind comparison tool),
`1539a22` (the two-ensemble driver), `e5e5cb3` and this update. Follows
`SESSION_HANDOFF_forefire_2026-09-21.md`, whose §10 set the goal: **test the
GRIB-driven path in more complicated terrain.** Done, on two fires.

**Headline:** **ForeFire tolerates WindNinja's direction error and is sensitive to
its speed bias**, and the speed bias does not go one way in complex terrain. Silver
(Selkirk Mtns) ran 24% slow and produced a fire 28% *smaller* than the coupled run;
Dome (Yosemite), rougher still, ran 18% fast and produced one 16% *larger*. Both
times the area ratio tracked the speed ratio, and both times a 36-49 deg direction
RMS moved the ensemble centroid under 520 m on a multi-kilometre fire. **Terrain
roughness predicts none of it** (§1).

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

    * Union: relief 21-93 m; std not computed.
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

This closes 09-18 §10 item 6 for the wind half. The perimeter-scoring half
(`dry_ff_compare.py`) is still missing.

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

## 6. Silver-specific caveats that do not affect the ratio

- **`fire_init` was not invoked for this fire**, confirmed by measurement: cells
  inside the 2026-08-27 perimeter are 9.6% `NFUEL_CAT == 14` against 13.3% outside.
  A burn mask would have made inside ~100%. Per JH the module can acquire IR
  perimeters and mask consumed fuel to unburnable, but it did not run here, so both
  models burn through the August scar. The domain is sub-alpine with considerable
  talus near treeline, which is why there is a natural cat-14 background at all.
- **Absolute areas are therefore not realistic for this fire.** The ratio is,
  because both models use the identical fuel map.
- **The only perimeter available predates the forecast by 24 days** (§6), so nothing
  here is scored against observation.

## 7. A perimeter can be much older than the forecast

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

## 8. The two fuel paths differ cell by cell

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

## 9. Open items

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
7. `dry_ff_compare.py`, the perimeter-scoring half of 09-18 §10 item 6, is still
   unwritten.
8. **Ensemble spread does not generalise** (§5). Three fires give three orderings:
   HRRR more dispersed on Red Bank and Dome, tighter on Silver. Do not read a
   spread difference as meaningful without more cases.
9. **Cotton 2 is deferred, deliberately, until the cache covers its wind shift.**
   `COTTON_2_2026-09-21_20_00_00_283A99DF-AFD8-4F6A-8CA4-F827DA0AFBCD`, Napa County,
   ignition 2026-09-21 21:26:54Z at 38.63331, -122.06932. `ZSF` std **184.4 m**,
   which fills the gap between Red Bank (31.5) and Silver (365) — the mid-range
   regime §5 implicates as the awkward one. Second California case.

   Its interest is a **~164 deg wind reversal**: steady SW 194-227 deg from 19:00
   09-21 through 11:00 09-22, swinging to 360 deg at 12:00 and then round through
   N, NE, E to SE by 00:00 09-23. **That reversal begins on the last hour the f03
   cache reaches**, so the whole post-shift regime is out of range today. Per JH,
   held until the cache is complete rather than run on the pre-shift window.

   Two details for whoever picks it up. The WRF field goes **incoherent through the
   transition** — R falls to 0.161 at 11:00 and 0.118 at 12:00, below the 0.3 floor
   — so direction statistics *across* the shift are meaningless regardless of data;
   compare the steady regimes either side instead. And the f03 cache advances about
   an hour per hour, so the full window to 2026-09-23 00:00Z became reachable around
   00:00Z on 09-23. Nothing prunes `ingest/HRRRA`, so the early hours will still be
   there.

   **14 netcdfs of the pre-shift window are already built** in
   `/home/jhaley/forefire/tests/ffwksp_cotton2_hrrr`. `forefire_grib` skips existing
   netcdfs, so resuming costs only the new hours — raise `--steps` and re-run the
   same command.

## 10. NEXT SESSION

1. **Finish Cotton 2** (§9 item 9). The cache should now cover its wind reversal,
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
   `poly_PolygonDateTime` against the forecast window first (§7).
