# ForeFire session handoff — 2026-09-28

The ForeFire pipeline now runs a forecast from a wrfxpy job file with no WRF-SFIRE
anywhere in the loop. Winds and fuel moisture both come from HRRR forecasts.

## 1. TL;DR

- **`src/ingest/forefire_forecast.py` is the missing orchestration** (`089fb97`). Job
  file in, perimeters out. Verified end to end on two fires.
- **All three inputs are HRRR now.** Winds from the forecast cache via WindNinja, fuel
  moisture from FMDA's HRRR region, fuels from LANDFIRE.
- **Moisture varies per step and per fuel category.** On GALENA_RX all nine steps got
  their own fuel table, with the 1-hour class recovering 0.081 -> 0.117 overnight.
- **`etc/forefire.json` moved from `wksp_fmda/CONUS` to `wksp_fmda/CONUS_HRRR`.** CONUS
  is RTMA-driven, and RTMA is an analysis with no hours past the present, so a forward
  forecast silently reverted to the default fuel table partway through its own window.
- **A genuinely forward run works**: ignition at T+30, every wind step from HRRR leads
  valid in the future, every moisture step from FMDA `-fNN` forecast files.
- **The moisture switch does not break the pSAF calibration (§5).** RTMA and HRRR differ
  by 0.2-0.4 percentage points in the quantity that reaches the fuel table, well inside
  the run-to-run noise. **But near extinction they disagree 9% of the time** in the
  Oklahoma panhandle, and pSAF cannot compensate for a fire that never starts.
- **Piece 4 is done: keep 0.337 (§6).** On HRRR moisture the batch gives geo mean 1.39
  against 1.54 paired, but the 25th percentile — the actual calibration target — moves
  only 1.11 -> 1.05, implying 0.326. Inside the noise.
- **Neither moisture product reproduces WRF-SFIRE's stalls better (§7).** Over the seven
  biggest movers RTMA matches on two, HRRR on one, neither on a fourth. BOON is the
  extreme case, not a representative one. The sign of the RTMA/HRRR difference predicts
  the direction of every area change, and that difference runs -0.064 to +0.042 per
  fire -- 10-20x §5's estimate, and **regionally signed**, so averaging domains hides it.
- **BOON alone looked like a counterexample to "HRRR is better" (§7).** RTMA put fuel 2 above extinction
  overnight and reproduced a stall WRF-SFIRE independently shows; HRRR missed it by 0.003
  and ran 2.94x against RTMA's 1.83x. At BOON HRRR is 0.03 drier than RTMA -- ten times
  the difference at either §5 sample domain, so §5 was under-sampled.
- **Caution: a regional moisture comparison under-predicted the impact.** §5's 0.2-0.4
  point difference implied the calibration was safe by a wide margin; the real aggregate
  move was 10-15% in area, because area does not respond proportionally to moisture near
  extinction. Re-run rather than infer.

## 2. What was built

`forefire_forecast.py`:

    PYTHONPATH=src:src/ingest python src/ingest/forefire_forecast.py <job.json> \
        [--steps 8] [--issue-utc ...] [--domain-km N] [--overwrite]

It reads the ignition, builds the per-step netcdfs through `forefire_grib`, writes the
`.ff` scripts through `make_script_set` (which applies the per-step moisture), runs
them, and collects the outputs. Every one of those pieces already existed; the only
thing joining them was a hand-written scratch script per fire.

**Job time formats.** Three spellings are in circulation and all three appear in real
files: `2026-09-27_21:08:00` in job files and workspace `input.json`, ISO with a `Z`,
and the relative `T+30` that the templates under `tests/` use. Accepting only the first
is how a template silently becomes a 1970 ignition.

**Domain.** Derived from the job's own grid — `cell_size / subgrid_ratio` is the fire
mesh WRF-SFIRE would have used, `cell_size * domain_size` its extent — so a ForeFire run
and a WRF-SFIRE run of the same job are comparable rather than accidentally different.
When that exceeds WindNinja's measured cell ceiling the **domain** is trimmed, not the
mesh: coarsening the solve degrades exactly the ridge-top winds that drive the fastest
spread, and an 8 h fire does not reach the edge of a 50 km box. A 61 km job at 25 m
(5.95 M cells) trims to 50 km.

## 3. Verified end to end

**GALENA_RX** (Oregon, 2026-09-27 21:08, 31 km at 25 m, 8 steps):

- winds from the 18z cycle leads 3-5, then the 00z cycle leads 0-4 — `find_gribs`
  switching to the fresher cycle as soon as it exists, which is the piece-2 rule working
- nine distinct fuel tables, 1-hour class 0.0807 -> 0.1166 overnight
- **3.1 -> 76.0 ha over 8 hours**, steady growth, plausible for a prescribed burn

**A forward run** from `tests/default_simple_forecast.json`, whose `T+30` resolved to 30
minutes ahead of the clock:

| | |
|---|---|
| wind steps | 06z cycle, leads 6-12, valid 12:00-18:00 — **steps 3-7 ahead of the clock** |
| moisture steps | 10z cycle, leads 2-8 — **every hour from a `-fNN` forecast file** |

That is the whole point of the exercise: nothing in the run waited on a WRF-SFIRE
forecast, and nothing was read from a file that only exists in hindsight.

**The forward fire does not spread, and that is correct.** Fuel 2 sits at Md
0.164-0.171 against its moisture of extinction of 0.150 for the first five steps — too
wet. The sixth step reaches 0.1426, below extinction but only barely, and Rothermel's
rate of spread goes to zero as Md approaches `me`, so no measurable growth in one hour
is the right answer rather than a stalled front. Colorado front range, late September,
1862-2190 m. **I first read this as the front failing to restart after a stall; the
measurements do not support that** — the front is present and valid at every step
through 18:00, and only one step was ever below extinction.

## 4. A correction carried over

Last session's §7 traced a bug where an out-of-range-dry 1-hour moisture class was
replaced by roughly the 100-hour value. The clamp fixed it. Worth keeping straight:
**the fix prevents the stall; it was never shown that a stalled front can restart.**
That question is still open and needs a case where moisture crosses `me` decisively
mid-run, not marginally.

## 5. RTMA vs HRRR FMDA — does the pSAF calibration travel?

Per JH: two 30x30 km domains, two weeks of hourly analyses, box means. The grids differ
(RTMA 1258x2145, HRRR 1038x1799) so there is no cell-by-cell comparison to make; a box
mean is also what ForeFire consumes, since `md_series` takes a `radius_km` mean around
the ignition. `src/ingest/fmda_compare_regions.py` (`c67eaf5`), 339 common hours,
2026-09-14 to 09-28.

| domain | quantity | RTMA | HRRR | HRRR-RTMA | RMS | corr |
|---|---|---|---|---|---|---|
| Sierra / DOME | 1-h | 0.0841 | 0.0887 | +0.0047 | 0.0097 | 0.938 |
| (37.57, -119.62) | **fuel 2 Md** | **0.0863** | **0.0882** | **+0.0019** | 0.0060 | 0.951 |
| Oklahoma panhandle | 1-h | 0.1244 | 0.1259 | +0.0014 | 0.0272 | 0.906 |
| (36.75, -101.60) | 100-h | 0.1821 | 0.1597 | **-0.0224** | 0.0374 | 0.931 |
| | **fuel 2 Md** | **0.1343** | **0.1307** | **-0.0036** | 0.0239 | 0.904 |

**The answer to JH's question is yes, with a caveat.** In the quantity that actually
reaches the fuel table the two products differ by **0.2 to 0.4 percentage points**, with
correlations of 0.90-0.95. That is far inside the ~3% per-fire reproducibility noise
(09-25 handoff §4), so the pSAF calibration travels and 0.337 should survive the switch.

**The caveat is that moisture does not act proportionally near extinction.** Counting
hours where fuel 2 could burn at all (Md < me = 0.15):

| domain | RTMA | HRRR | disagree |
|---|---|---|---|
| Sierra / DOME | 339/339 | 339/339 | **0 (0%)** |
| Oklahoma panhandle | 218/339 | 233/339 | **29 (9%)** |

The Sierra sits at Md ~0.087, far below extinction, so the products are interchangeable
there. The panhandle sits at ~0.13 against a threshold of 0.15 — close enough that the
same 0.4-point difference flips the verdict in 9% of hours, and HRRR burns in **more**
of them. **pSAF cannot compensate for that**: it scales the rate of a fire that spreads,
and has no effect on one that never starts. So the calibration travels for fire *size*
while the burn/no-burn boundary shifts slightly toward more fire in marginal regimes.

**Where they differ most is the 100-hour class** (HRRR 2.2 points drier in the
panhandle, against 0.1 for the 1-hour). Fuel 2 weights 100-h at only 0.1429, which is
why the derived Md difference stays small. Fuels that weight the coarse classes more
heavily would show a larger divergence, and that is worth checking before trusting this
result across all 13 categories.

## 6. Piece 4 done — the pSAF re-check on HRRR moisture

Batch re-run with `fmda_dir` on `CONUS_HRRR`, scored against the preserved RTMA baseline
(`src/ingest/reference/psaf337_rtma_baseline_2026-09-25.txt`). **Paired on the 42 fires
present in both:**

| | RTMA | HRRR | change |
|---|---|---|---|
| geometric mean | 1.54 | 1.39 | **-10.0%** |
| median | 1.55 | 1.31 | -15.2% |
| below 1.0x | 7 | 8 | |
| within factor 2 | 34 | 33 | |

69% of fires produce less fire under HRRR moisture, with a per-fire log-change sd of
**0.300** -- individual fires move about +/-35%, extremes from -64% (HOP-PATTERSON_PILE_RX)
to +93% (Merit_Creek). That is an order of magnitude above the ~3% reproducibility noise,
so it is signal.

**Keep 0.337.** The calibration posture was the 25th percentile at 1.0, and it is still
there: 1.11 -> 1.05, implying pSAF **0.326**, a 3% move inside the noise. Restoring the
*median* would want 0.371, but the median was never the target.

**A prediction that was too optimistic, worth recording.** From the two-domain moisture
comparison in §5 I said 0.337 "should survive" on the strength of a 0.2-0.4 point mean
difference. The actual aggregate move is 10-15% in area -- the right conclusion, but by a
narrower margin than that comparison implied. The reason is the same non-linearity the
Oklahoma panhandle exposed: **area does not respond proportionally to moisture near
extinction**, so a small mean difference flips fires across the threshold and moves area
far more than the mean difference suggests. A regional moisture comparison is therefore a
weak predictor of a calibration's stability, and is not a substitute for re-running.

The HRRR set also has slightly better factor-of-2 coverage (79% vs 76%) on the full sets
despite producing less fire, which is mildly in favour of JH's view that the HRRR product
is the better one.

**Overwriting, and where the old perimeters are.** `run_forecasts` rewrites
`<wksp>/forefire` in place, so the 09-25 RTMA-moisture perimeters were overwritten for
the 50 re-run workspaces. **JH preserved them**: `/home/jhaley/tmp/old_<grid_code>_*.kml`,
350 `_final` and 352 `_steps`, covering every fire that moved most in the table above.
The `_steps` files are the useful ones -- they show *when* the two runs diverge, not just
the endpoint, which is what distinguishes a moisture-threshold stall from a uniform
scaling.

    /home/jhaley/tmp/old_BOON_2026-09-21_21_00_00_703D4319-...-DFFDCFDF888A_final.kml

The matching new perimeters are `<wksp>/forefire/*_final.kml` and `*_steps.kml`. The
scored table is committed under `src/ingest/reference/` so the next A/B does not depend
on a scratchpad. For side-by-side runs in future use `sweep_ff_params`, which takes
`dest_dir` and `suffix` and keeps both.

## 7. BOON — a counterexample, and a correction to §5

JH noticed the RTMA-driven BOON run nearly stalled around step 45 with growth only in
patches, then made big jumps that looked questionable. Both halves check out, and the
conclusion runs against the general prior.

**The stall is fuel 2 crossing extinction.** Reconstructing both moisture series for
BOON's 52 steps (Mendocino County, mixed fuels: 2, 5 and 10 near the ignition):

| | fuel 2 (me 0.15) | fuel 5 (me 0.20) | fuel 10 (me 0.25) |
|---|---|---|---|
| RTMA | 0.086-0.185, **blocked 22 of 52 steps** | 0.075-0.191, never | 0.114-0.170, never |
| HRRR | 0.059-**0.147**, never blocked | 0.046-0.149, never | 0.089-0.139, never |

The blocked steps are overnight local (PDT = UTC-7). Fuels 5 and 10 keep burning
throughout, which is exactly JH's "only small growth in some areas": the grass shuts off
and the fire creeps on in the heavier fuels.

**The jumps are not an artifact.** WRF-SFIRE does the same thing on its own: growth falls
from 21.7 ha/h at 20:00 local to 2.2-3.8 ha/h from 03:00 to 10:00, then goes 3.4 -> 29.9
-> 63 -> 81 -> 96 -> **110.8 ha/h** through the afternoon. That is the burning period, and
the RTMA run reproduced it. The HRRR run, which never stalls, is smooth throughout and
shows no such structure.

**On this fire RTMA is the better product.**

| | final | vs WRF-SFIRE |
|---|---|---|
| WRF-SFIRE | 603 ha | — |
| ForeFire, RTMA moisture | 1103 ha | 1.83x |
| ForeFire, HRRR moisture | 1776 ha | 2.94x |

HRRR misses the stall by **0.003** -- its fuel 2 peaks at 0.147 against a threshold of
0.150. RTMA did not get closer by luck; it reproduced a physical overnight stall that
WRF-SFIRE independently shows.

**The correction to §5.** At BOON, HRRR runs **0.027-0.038 drier than RTMA** throughout,
with a similar diurnal amplitude -- a regional dry bias, roughly **ten times** the
difference measured at either sample domain (+0.002 Sierra, -0.004 panhandle). So §5's
"the products differ by 0.2-0.4 points" was under-sampled: the difference is regionally
variable and two domains did not capture it. That also explains why the batch-wide area
change (§6, 10-15%) so far exceeded what §5 implied -- it is not only the extinction
non-linearity, the underlying moisture difference is itself much larger in some regions.

**What this does not settle.** §6's aggregate still says keep pSAF 0.337, and the 25th
percentile is a robust statistic. But "HRRR is the better product" is not established
fire by fire, and BOON is a clear case against it. The useful next test is not another
regional moisture comparison but a **per-fire** one: for the fires that moved most in §6,
which product's Md series matches the stall structure WRF-SFIRE shows? That is a
falsifiable question and BOON answers it against HRRR.


### The seven big movers — BOON was not representative (2026-09-29)

Applying the same reconstruction to every fire that moved most in §6
(`src/ingest/ff_moisture_stall_test.py`, `7135030`):

| fire | dArea | fuel | Md RTMA | Md HRRR | HRRR-RTMA | mechanism |
|---|---|---|---|---|---|---|
| Merit_Creek | **+93%** | 10 | 0.2560 | 0.1925 | **-0.064** | threshold, 27 steps |
| BOON | **+61%** | 2 | 0.1367 | 0.0992 | **-0.038** | threshold, 22 steps |
| Menogyn_Piles_RX | **+45%** | 10 | 0.1783 | 0.1669 | -0.012 | smooth |
| HOP-PATTERSON_PILE_RX | -64% | 8 | 0.1802 | 0.1830 | +0.003 | smooth |
| D6_RX_Treasured_Landscapes | -51% | 5 | 0.0983 | 0.1401 | **+0.042** | threshold, 9 steps |
| Little_Blue_Stem_RX | -49% | 10 | 0.1485 | 0.1880 | **+0.039** | smooth |
| Silver | -37% | 5 | 0.0821 | 0.1205 | **+0.039** | smooth |

**The sign of the moisture difference predicts the direction of the area change in all
seven.** HRRR drier -> more fire; HRRR wetter -> less fire.

**§5 was worse than under-sampled.** Per-fire differences run -0.064 to +0.042 against
the +/-0.004 measured at the two sample domains -- 10 to 20x larger -- and, more
importantly, the difference **changes sign regionally**. Averaging over domains would
have hidden it completely. A two-domain comparison is not a valid way to assess this.

**Only 3 of 7 are threshold crossings** (Merit_Creek, BOON, D6). The other four are
smooth offsets, so their perimeters should scale rather than change shape -- a concrete
prediction for the KML comparison.

**Correction to the §7 lean.** §7 suggested RTMA might be the better product. Over seven
fires that does not hold. Of the four where WRF-SFIRE genuinely stalls (>30% of hours
below a quarter of peak growth): RTMA matches on **Merit_Creek** and **BOON**, HRRR
matches on **D6**, and neither explains **HOP-PATTERSON** (WRF stalls 41%, neither product
blocks). Each also blocks spuriously once where WRF shows no stall -- RTMA on Menogyn (16
steps), HRRR on Little_Blue_Stem (14). That is 2-1 and 1-1 on seven fires: noise. **BOON
is the extreme case** -- WRF stalls 77% of its hours there against 0-41% elsewhere --
**not a representative one.** Neither product is systematically better at reproducing
stalls, and §6's aggregate remains the only solid basis for preferring HRRR.

**One case unexplained.** HOP-PATTERSON_PILE_RX moved -64% on a +0.003 moisture
difference in its dominant fuel. Every other fire's magnitude tracks its moisture
difference; this one does not. Either a non-dominant fuel carries it or something other
than moisture is involved.


### Corrected: area-weight the domain, and use the moisture *ratio* (2026-09-29)

The table above picked each fire's fuel from a 50x50 cell window at the **ignition**.
That is not the fuel the fire spreads through. On HOP-PATTERSON it selected fuel 8 -- 28%
of the domain and the one fuel where the two products agree to 0.003 -- while **fuel 5
covers 48% and differs by 2.8x**. That is the whole of the "unexplained" case.

Area-weighting over the whole domain instead:

| fire | dArea | Md RTMA | Md HRRR | **ratio** | crossings |
|---|---|---|---|---|---|
| BOON | +61% | 0.1383 | 0.1032 | **0.75** | 22 |
| Merit_Creek | +93% | 0.2058 | 0.1661 | **0.81** | 39 |
| Menogyn_Piles_RX | +45% | 0.1703 | 0.1556 | **0.91** | 10 |
| Little_Blue_Stem_RX | -49% | 0.1346 | 0.1789 | **1.33** | 40 |
| Silver | -37% | 0.0961 | 0.1298 | **1.35** | 0 |
| D6_RX_Treasured_Landscapes | -51% | 0.1034 | 0.1412 | **1.37** | 12 |
| HOP-PATTERSON_PILE_RX | -64% | 0.0821 | 0.1195 | **1.46** | 0 |

**Sorted by ratio the seven separate perfectly by sign**, and HOP-PATTERSON -- the
previous anomaly -- is now the extreme point at both ends: largest ratio, largest change.

**Two corrections that matter beyond this table.**

1. **Use the ratio, area-weighted over the domain**, not an absolute difference at the
   ignition. Md sampled at a point says little about what the fire runs through.
2. **Threshold crossings are not the main mechanism.** HOP-PATTERSON and Silver have
   **zero** crossings and still move -64% and -37%. Rothermel's rate sensitivity in dry
   fuel does the work; crossings are an additional effect on top, not the driver. The
   earlier "3 of 7 are threshold crossings" was wrong on both counts.

**JH's reading of HOP-PATTERSON holds in substance** -- a moisture difference large
enough to stop the HRRR fire reaching terrain where it could run uphill -- but the gate
is spread *rate*, not fuel gaps. `ff_domain_connectivity.py` measures the latter and
clears this domain: 97.7% burnable, 100% of it reachable from the ignition, fragility
0.08. It does flag the genuine connectivity cases -- **Mud_Bayou at 47.6% burnable and
fragility 0.72**, the known bridge crossing, and **Silver at fragility 1.00** -- so it is
worth running before fitting anything, since those fires change area for reasons that
have nothing to do with spread rate.

(Its relief column reads 11986 m on Mud_Bayou, which is a fill value in the altitude
field, not terrain. Do not trust that column without checking.)


### Connectivity over the whole batch (2026-09-29)

`ff_domain_connectivity.py` run over all **361** retained fire domains.

| | min | 10th | median | 90th | max |
|---|---|---|---|---|---|
| burnable | 0.347 | 0.880 | 0.981 | 0.998 | 1.000 |
| reach | 0.000 | 0.776 | 0.999 | 1.000 | 1.000 |
| fragility | 0.008 | 0.015 | 0.053 | 0.350 | 1.000 |

**70 of 361 flagged (19%)** at burnable<0.80, reach<0.90 or fragility>0.50. The worst:
STEEL_PASS (43.5% burnable, **0.0% reachable**), Silver (fragility 1.00), Foster (0.96),
Bear_River_MBR_Unit_2c_RX (0.3% reachable), **Mud_Bayou (47.6% burnable, 0.72)** -- the
known bridge-crossing case, correctly picked out without being told about it.

**STEEL_PASS explains itself.** It is one of the fires logged as "no common valid time"
and therefore unscorable; with none of its burnable fuel reachable from the ignition
there was never going to be a fire to score.

### Do flagged fires distort the calibration? Noise, not bias.

| | flagged (n=14) | clean (n=38) |
|---|---|---|
| geometric mean | 1.398 | 1.458 |
| median | 1.315 | 1.335 |
| **sd of log ratio** | **0.698** | **0.484** |
| within factor 2 | 71% | 79% |

Central tendency matches; **scatter is 45% larger**. Excluding them moves the implied
pSAF from 0.328 to 0.323, against 0.337 in use -- negligible.

**So: keep them in the global calibration, drop them from anything fitted per fire.**
Excluding 27% of the scored set to move a constant by 1.5% is a bad trade. But for the
NGFS-predicted per-fire pSAF (09-24 handoff §9 item 2) the extra scatter is precisely
what would corrupt the fit, and that is the same reasoning already recorded for
Mud_Bayou. Run this first and exclude the flagged fires there.

**A bug in the metric, fixed.** The first run reported fragility down to **-0.583**, which
is impossible by its definition. When the ignition cell fell out of every component after
erosion, the nearest-surviving-component fallback could select a *different and larger*
component, giving reach2 > reach. Now intersected with the original reachable component
first; the range is [0.008, 1.000] and STEEL_PASS moved to the top of the ranking.


### Merit_Creek was a failed run, not a moisture effect (2026-09-29)

JH noticed Merit_Creek's `_steps.kml` showed only a few steps. It is worse than a small
fire -- the RTMA run **died two hours in**:

| | perimeters | ForeFire ha | WRF ha at the common time | ratio |
|---|---|---|---|---|
| RTMA | **3** of 53 step files | 0.7 | 1.3 | 0.55 |
| HRRR | 53 | **207.4** | 194.8 | **1.06** |

The scorer takes the latest time **both** models reach, so the RTMA run was compared
against WRF-SFIRE two hours in: small against small, giving a ratio of 0.55 that looked
entirely ordinary. Comparing it with the healthy run's 1.06 produced a "+93% moisture
effect" that put Merit_Creek at the top of the movers list and into the per-fire stall
analysis. **None of that was real.** JH confirms the HRRR run tracks WRF-SFIRE closely,
which the numbers agree with -- 207.4 against 194.8 is one of the closest matches in the
batch.

**This does not touch §6.** Only 1 of the 42 paired fires was scored at a different valid
time. Excluding it the effect is slightly *larger*: geo mean -11.7% against -10.0%,
median -15.9% against -15.2%. Keep pSAF 0.337.

**It does invalidate Merit_Creek's place in §7.** It was one of the two fires counted as
"RTMA reproduces the stall". RTMA did not reproduce a stall there; it produced a
degenerate run. The §7 tally becomes RTMA 1 (BOON), HRRR 1 (D6), neither 1
(HOP-PATTERSON) -- which strengthens rather than changes the conclusion that **neither
product is better**.

**Guard added** to `ff_batch_ratio.py`: skip a fire whose ForeFire track covers under
half the WRF window, or under half of whose step files yield a perimeter. It fires on
nothing in the current batch -- the genuinely truncated runs are already caught as "no
common valid time", and Merit_Creek's RTMA output has been overwritten -- but the failure
is silent and the ratio it produces is plausible, which is exactly why it needs a guard
rather than an eye.

**A bug in that guard, caught before it did damage.** The step-file test globbed
`*.geojson`, but a run writes one `final` geojson per step and `track()` skips those, so
the denominator was doubled. It flagged Bear_River, COTTON_2 and Huckleberry as truncated
at "47 of 108" when 47 of 54 is healthy. Now counts only what `track()` reads.

## 8. Open items

1. ~~**Piece 4, the pSAF re-check.**~~ **DONE — see §6. Keep 0.337.**
2. ~~**FMDA forecast mode has no schedule.**~~ **DONE** — cron entry added:

       30 2,8,14,20 * * *   ./hrrr_cycler.sh f CONUS_HRRR

   One run per 6-hourly HRRR cycle. At 02/08/14/20 local (MDT) `cycle_start` lands on
   00/06/12/18Z respectively, ~2 h after the cycle so the leads have published, giving
   11 leads and ~11 h of forward moisture per run.

   **Correction to what this item first said.** It advised a *separate* flock. That is
   wrong and would have been a data race: forecast mode is **not read-only**. Lead 1
   seeds from lead 0, which is the cycle's *analysis* file, and it writes that file
   itself when it is missing — the 2026-09-28 run ran DA three times and created
   `fmda-CONUS_HRRR-20260928-10.nc`. It therefore shares `hrrr_cycler_fmda.lock` with
   the hourly analysis. Minute 30 because 11 leads take ~35 min and the analysis fires
   at :27, leaving 57 minutes of clear air.
3. **The last step of a run loses its front** — writes an empty geometry with the
   timestamp reset to the anchor. Seen on Eagle_Springs and again here. Harmless to the
   forecast, but it is the same signature as a genuinely broken chain, so it makes real
   breakage harder to spot.
4. **`ff_growth.track` collapsed seven step geojsons into two rows** on the forward run.
   Did not affect the result — the areas were read from ForeFire's own logs — but it
   would corrupt scoring, so it should be fixed before the next batch is scored.
5. Items 3-6 of the 09-24 handoff §9 remain untouched: the fuel 2/9 transposition
   upstream, the false "unknown moisture source" log line, scoring against observations
   rather than only WRF-SFIRE, and the two fires with no overlapping valid time.

## 9. NEXT SESSION

**Look at the preserved perimeters (§6).** JH's own next step. The old RTMA-moisture
runs are in `/home/jhaley/tmp/old_<grid_code>_*.kml` and the new HRRR ones in
`<wksp>/forefire/`. Start with the fires that moved most in §6 -- Merit_Creek (+93%),
BOON (+61%), Menogyn_Piles_RX (+45%), HOP-PATTERSON_PILE_RX (-64%),
D6_RX_Treasured_Landscapes (-51%), Little_Blue_Stem_RX (-49%), Silver (-37%) -- all of
which are present in both sets.

**What to look for, and a falsifiable prediction.** A smooth moisture offset scales a
perimeter roughly isotropically; a *threshold* crossing changes its **shape**, by
switching off one fuel category while others keep burning. By the corrected §7 analysis
**BOON, Merit_Creek, Little_Blue_Stem_RX, D6 and Menogyn cross thresholds; HOP-PATTERSON
and Silver do not** -- so those last two should scale rather than change shape. If they
show shape changes instead, the smooth-offset reading is wrong. The `_steps` files are
what distinguish the two, because they show when the divergence starts.

**The per-fire moisture test is done** (§7, 2026-09-29) -- do not redo it. Over all seven
movers **neither product reproduces WRF-SFIRE's stalls better**, and BOON was the
extreme case rather than a representative one. What came out of it instead: the
HRRR/RTMA moisture **ratio**, area-weighted over the whole domain, separates all seven
fires by sign, and threshold crossings are *not* the main mechanism -- HOP-PATTERSON and
Silver have zero crossings and still move -64% and -37%.

**Run `ff_domain_connectivity.py` before any per-fire fit.** Over all 361 domains it
flags 70 (19%) as connectivity-gated -- fires whose area changes because the fire reaches
a region or does not, which has nothing to do with spread rate. It finds Mud_Bayou
unprompted (47.6% burnable, fragility 0.72) and explains STEEL_PASS, long logged as
unscorable, as having **0.0% of its burnable fuel reachable from the ignition**.

Those fires are **noise, not bias**: the same central tendency as clean fires but 45%
more scatter in the log ratio (sd 0.698 against 0.484). So **keep them in the global
calibration** -- excluding 27% of the scored set to move pSAF from 0.328 to 0.323 is a
bad trade -- and **exclude them from the NGFS per-fire pSAF fit** (09-24 handoff §9 item
2), where scatter is exactly what corrupts a fitted relationship. That is the same
reasoning already standing for Mud_Bayou, now measurable across the batch rather than
case by case.

**Re-run `fmda_fcst_skill.py` once the cron has accumulated cycles.** This session's
forecast-vs-analysis result came from 2 cycles with n=1 at leads 3-9, and JH's read of it
was right: the apparent wet drift was the period, not the system. The two cycles go
opposite ways (+0.061 and -0.007 at the 10-h class), and the archive confirms 09-25 was
drying 1.1 sd faster than normal while 09-28 was moistening. The real signature is that
the forecast **lags the analysis tendency in either direction** -- which matters more
than a bias, because it cannot be corrected away and is largest during exactly the
rapid-drying events that drive fire. At 4 cycles/day x 32 leads the cron yields ~128
pairs/day. **Score against the analysis tendency, not lead alone.**

**The full-duration comparison JH asked for**: `forefire_forecast.py --steps 27` covers
the median WRF-SFIRE run, and both the GRIB cache (+33 h) and the 33 h moisture span it.

**Watch the FMDA forecast cron.** Its first firing was silently cancelled -- `flock -n`
lost the lock to a 12-minute analysis under load 37 -- and it now uses `-w 1200` with a
skip log, but that has not yet been seen to work.

Commits remain unpushed on `james_ngfs` and `release-fmda-fixes`; JH pushes them manually
(SSH password).
