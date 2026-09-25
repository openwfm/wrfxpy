# ForeFire session handoff — 2026-09-25

A short session: the scoring pass that the 09-24 handoff's §10 asked for, plus writing the
result down where it will survive. One bug fixed along the way (§6).

## 1. TL;DR

- **pSAF 0.337 is verified and the global-pSAF question is closed.** 65 fires: geometric
  mean **1.55**, median **1.53**, 17% below observed, 75% within a factor of two, against
  a prediction of 1.46 / 1.54 / 26%.
- **The `pSAF^1.54` exponent holds over a 1.9x move**, as it did over the 0.6 -> 0.175
  move. It can be trusted for future adjustments.
- **The templates now explain the number.** Both `.ff` templates carry a comment block
  recording the calibration, so 0.337 does not read as a magic constant.
- **Correction: the size dependence has not vanished.** It is `area^-0.198` over fires of
  10 ha or more — half the old `-0.40`, still real. An all-fires fit says `-0.050`, and
  that is an artifact of two tiny fires.
- **`ff_batch_ratio.py` could silently score the whole archive** when given a bare number
  of hours that matched a filename. Fixed (§6). It did not affect any number here.
- **Reruns are not reproducible to better than ~3%**, so differences below about 5%
  between ForeFire configurations are noise (§4).
- **Nothing was re-run.** The scored population is the 0.337 batch from 09-24 plus the
  nightly cron; the only run today was a 27 ha smoke test to prove the templates parse.

## 2. The result

Scored with `src/ingest/ff_batch_ratio.py`, ForeFire against WRF-SFIRE burned area at the
latest valid time both models reach:

    cd /data/jhaley/wrfxpy && PYTHONPATH=src:src/ingest \
      python src/ingest/ff_batch_ratio.py 18

The argument is hours. `ff_batch_ratio.py` normally takes the batch script and uses its
mtime to bound the window, but `/tmp/psaf_work/psaf337_batch.py` had been cleaned up. The
templates' own mtime is the same boundary — the moment pSAF became 0.337 — so hours from
that works and is not a guess. **If a batch script still exists, pass it**; the reason
that option exists is that a guessed window let two stale pSAF-0.6 fires into the 0.175
table and moved its geometric mean from 0.531 to 0.577.

| | predicted (09-24 §2) | measured |
|---|---|---|
| geometric mean | 1.46 | **1.55** |
| median | 1.54 | **1.53** |
| below 1.0x | 26% | **17%** |
| within a factor of 2 | — | **75%** |
| n | 57 | **65** |

The median is essentially exact. The batch runs slightly more conservative than aimed
for, which is the harmless direction given JH's preference for overprediction.

### Why stop at 0.337

Rescaling the scored ratios by the exponent shows it sitting at the peak of the
within-2x curve:

| pSAF | geo mean | median | below 1.0x | within 2x |
|------|----------|--------|-----------|-----------|
| 0.300 | 1.30 | 1.29 | 26% | 72% |
| **0.337** | **1.55** | **1.53** | **15%** | **75%** |
| 0.380 | 1.85 | 1.83 | 11% | 60% |
| 0.420 | 2.15 | 2.13 |  9% | 43% |

Moving up to 0.38 buys 6 points fewer underpredictions and costs 15 points of accuracy.
There is no argument for further global tuning.

## 3. Correction — the size dependence is halved, not gone

I reported mid-session that it had collapsed (`ratio ~ area^0.02`). **That was wrong.**
The all-fires fit is flat because two sub-10-ha fires (0.7 and 1.3 ha) sit at the far
left of a log axis and act as leverage points. Over the 63 fires of 10 ha or more:

    ratio ~ area^-0.198     log-log r = -0.478

against `area^-0.40` at pSAF 0.6. Real, and about half as steep.

**The residual is not a clean power law**, which matters before anyone fits one:

| WRF area | n | geo mean ratio |
|---|---|---|
| <200 ha | 20 | 1.71 |
| 200-600 ha | 25 | 1.77 |
| >600 ha | 20 | 1.18 |

Non-monotonic. Small and mid-size fires behave alike; it is the **large** fires that are
already close to right. A single exponent describes that shape badly.

**Method note worth keeping.** Both errors this session — this one and the 0.175/0.337
mix-up in §5 — came from trusting a summary statistic without checking what was in the
sample. A log-log slope over a range that spans 0.7 ha to 6000 ha is dominated by its
endpoints; bin the data before believing the fit.

## 4. Template comment

`/home/jhaley/forefire/tests/ignition_template.ff` and `restart_template.ff` — **outside
the git repo**, referenced from `etc/forefire.json` under `templates.ignition` and
`templates.restart`. This handoff is therefore the only version-controlled record of what
they say. Both now carry, directly above the parameter:

    # Calibrated against WRF-SFIRE burned area over a 65-fire CONUS batch, 2026-09-24/25.
    # 0.6 (ForeFire default) overpredicted 3.85x; 0.175 overcorrected to 0.58x.  0.337 puts
    # the 25th percentile at 1.0x: geo mean 1.55, median 1.53, 17% of fires below observed,
    # 75% within a factor of two.  Area responds as pSAF^1.5, so changes here are steep.
    # Deliberately biased high -- per JH, overprediction is the safer error for a forecast.
    # Caveat: this is model-vs-model.  WRF-SFIRE is not truth; on Dome it ran 0.68x an IR
    # perimeter.  Re-derive against observations before treating 0.337 as physical.

Inserted with a `sed` address matched on the parameter *name*, not a line number, so it
stays put if the file is reordered. The two templates still differ only in
`minSpeed=0.001` and the ignition-vs-restart block, which is how they should differ.

Verified by re-running `RX_Stray_Creek__LCRD` (27 ha) under the cron `flock`: 51 steps,
all "Success: Process completed", no "wrong input file, check your settings" — the
generated `.ff` scripts parse and the step chain completes.

**The rerun did not reproduce exactly, and that is worth knowing.** Same pSAF, same
inputs, same valid time (09-24 21:00): 270.8 ha against the 261.8 ha scored a few hours
earlier, a **3.4%** difference, ratio 9.96 against 9.63. Not diagnosed. The likely
candidate is the hourly FMDA files — `md_series` reads them per step and the FMDA cycler
rewrites them, so a rerun can pick up revised moisture for hours already simulated;
front-tracking nondeterminism is the other possibility.

Either way it sets a noise floor: **differences below about 5% between two ForeFire
configurations are not evidence of anything.** The mesh-insensitivity result from 09-22
(<2% across 25/50/100 m) sits below that floor and should be re-read as "no detectable
effect", not "a 2% effect".

## 5. Correction — which value was actually running

Worth recording because it cost time. JH's instruction this morning was "leave it set at
0.175", but the templates had been at **0.337** since 11:04 on 09-24 — 0.175 was scored,
found to underpredict (geo mean 0.577, 90% of fires below 1.0x) and raised the same day.
The cron forecasts JH looked at and called "a big improvement" were the 0.337 ones, so
the instruction and the observation agreed; only the label was off.

I also mis-scored the batch twice before getting it right, both times by reading the
wrong files:

- the per-step geojsons in the **workspace root** are bare GeoJSON Features, not
  FeatureCollections, so a reader expecting `features` finds nothing;
- the combined `<name>.geojson` in the root is **WRF-SFIRE's own** perimeter export, not
  ForeFire's. Scoring against it compares WRF-SFIRE with itself and returns a suspiciously
  tidy 1.01-1.19 — the ~5% gap between TIGN_G and the polygonised perimeter.

ForeFire's perimeters are the geojsons in `<wksp>/forefire/`. `ff_batch_ratio.py` already
reads the right ones; the lesson is to use it rather than writing another ad-hoc scorer.

## 6. Bug fixed — `ff_batch_ratio.py` could silently score the whole archive

Its argument is "a batch script, or a number of hours", and it decided which by calling
`osp.exists()` **first**. A stray executable named `1` has been sitting in the repo root
since 2025-12-04, so `ff_batch_ratio.py 1` resolved to that file, took its mtime as the
window start, and scored **125 fires spanning three different pSAF values** — reported,
with a straight face, as the last hour. Geometric mean 1.52, median 1.26: plausible
enough to be believed.

Fixed by trying `float()` first and only falling back to a path, with an explicit error
if the argument is neither. A bare number can no longer be captured by a file that
happens to share its name.

The `18.01`-hour run that produced every number in §2 was never affected — nothing is
named `18.01` — and re-running the 1 h check after the fix returns exactly the one fire
that was re-run. **Worth deleting `/data/jhaley/wrfxpy/1`**; it is an old scratch script,
but it is not this tool's job to be robust against it and now it is.

## 7. Open items

Items 3, 4, 5 and 6 of the 09-24 handoff's §9 are untouched and still open: the fuel 2/9
transposition upstream, the false "unknown moisture source" log message, scoring against
observations rather than only WRF-SFIRE, and the three fires with no overlapping valid
time (Nethery, STEEL_PASS, Eagle_Springs).

Item 1 is closed. Item 2 — per-fire pSAF from NGFS detection fields — is the next real
move, now with a concrete target (§3).

## 8. NEXT SESSION

Item 2, but weigh item 5 first. The entire calibration is model-vs-model and WRF-SFIRE is
not truth; on Dome it ran 0.68x an IR perimeter while ForeFire ran 6.0x. A **global**
constant fitted against WRF-SFIRE inherits its bias on average, which is easy to revise
later. A **per-fire** correction fitted the same way would bake that bias in fire by
fire, which is not. If there are enough usable IR perimeters to fit against instead, that
ordering is worth the delay.

Commits remain unpushed on `james_ngfs` and `release-fmda-fixes`; JH pushes them manually
(SSH password).
