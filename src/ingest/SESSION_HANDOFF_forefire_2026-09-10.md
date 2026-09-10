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

## 7. Open items

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
