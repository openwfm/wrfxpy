# ForeFire session handoff — 2026-10-07

Developed land in LANDFIRE was being handed burnable fuel, so fires crossed interstates
and parking lots. This builds a remapped fuel raster, wires it into both the WRF-SFIRE
and ForeFire paths, and tests it on PICO.

**Section numbers continue from the 2026-10-06 handoff**, which continues 10-05 and 10-01.
Read this as the 10-07 continuation.

## 5s. Developed land remapped to no-fuel (code 91 -> 14)

### The defect

`forefire_grib.remap_fuel` gives LANDFIRE code **91 (urban/developed)** the value of its
nearest burnable neighbour, as geogrid's `var_wisdom` fill does, and 91 is **not** in
`BARRIER_CODES`, so the `-r max` pass does not restore it either. Every road and parking
lot therefore burns. §4 recorded this as deliberate ("highways are deleted on purpose")
and §5q showed its cost directly: on Dallas, **100% of the observed barrier was burnable
in the model**.

### What the data supports

An **NLCD 2024 land cover** raster is on disk
(`landfire/Annual_NLCD_LndCov_2024_CU_C1V1`), 30 m Albers like FBFM13 and grid-aligned
with it to an integer cell offset, so the two can be compared without resampling. It
separates developed intensity, which FBFM13's single code 91 cannot.

| domain | FBFM13 code 91 | NLCD 23+24 (med+high) | NLCD 21-24 (any dev) |
|---|---|---|---|
| **Pulaski** (industrial Little Rock) | **51.8%** | 30.8% | 66.5% |
| **Dallas** (ranch land) | 0.04% | **0.00%** | 0.76% |

Where NLCD calls a cell medium-or-high-intensity development, FBFM13 calls it 91 in
**98.1%** of cases, so code 91 is a reliable developed indicator and the only question is
which developed cells should be non-burnable.

**NLCD would not have fixed Dallas.** It finds *no* medium/high development there -- the
barrier is an unimproved ranch track, too narrow and too unpaved for any 30 m land-cover
product. **This work fixes interstates and parking lots, not rural roads**, and Dallas
stays unexplained.

### What was built, per JH: all of code 91 -> 14

- **`landfire/LFmosaic_FBFM13_CONUS_noburn91/Tif/LFmosaic_FBFM13_CONUS_noburn91.tif`**
  1.68 GB. **458,359,208 cells changed, 2.887% of CONUS.** Strip-wise copy preserving
  geotransform, projection, Int16, LZW, 256x256 tiling and NoData 32767, so WRFXPY's usual
  routine reads it unchanged. Built in 3.7 min by `scratchpad/remap91.py`.
- **`etc/vtables/geo_vars.json_mosaic_noburn91`** -- the geogrid path. Identical to
  `_mosaic` but for `NFUEL_CAT`.
- **`fuel_tif` config key** in `etc/forefire.json`, threaded through both
  `build_step_ncs` call sites, defaulting to `FUEL_TIF`. **Needed because ForeFire reads
  the raster directly and does not go through `geo_vars.json`** -- without it the geo_vars
  change would alter WRF-SFIRE while ForeFire kept the old fuels, confounding the
  comparison. It also makes the change a toggle rather than an edit.

### Verified, both paths

| | burnable | non-burnable |
|---|---|---|
| ForeFire `remap_fuel`, original | 97.73% | 2.27% |
| ForeFire `remap_fuel`, remapped | **69.31%** | **30.69%** |
| geogrid `geo_em` NFUEL_CAT, original | 96.9% | 3.1% |
| geogrid `geo_em` NFUEL_CAT, remapped | **70.9%** | **29.1%** |

The two paths agree to within a point, so both models now see the same landscape. Over
the PICO domain 303,232 cells changed, **all were 91, all are now 14, none remain**; what
they had been is fuel 5 (78.3%), fuel 8 (14.1%) and fuel 2 (6.8%).

### A prediction that was wrong, and why that matters

`ff_domain_connectivity` on the remapped PICO domain:

| | burnable | reach | fragility |
|---|---|---|---|
| original | 97.7% | 99.9% | 0.05 |
| remapped | 69.3% | **0.7%** | **1.00** |

Reach of 0.7% says only 0.7% of burnable ground connects to the ignition -- worse than
PIUTE's 79.8%, the most gated domain in the 391-fire archive -- and predicts near-total
confinement. **WRF-SFIRE contradicts it**: at 10-07 01:00 the remapped run is at 90.4 ha
and still growing.

So the erosion-based fragility measure **overstates confinement for a level-set model**,
which can round the ends of a linear barrier in a way the measure does not anticipate.
Worth knowing before that measure is used to exclude domains from calibration (10-01 §6
item 5).

### The barrier signature

WRF-SFIRE, same fire, same weather, fuel raster the only difference:

| valid | original ha | no_burn91 ha | ratio |
|---|---|---|---|
| 10-06 20:00 | 9.3 | 8.9 | 0.95 |
| 10-06 21:00 | 19.2 | 18.9 | 0.98 |
| 10-06 22:00 | 50.2 | 49.6 | **0.99** |
| 10-06 23:00 | 86.1 | 72.0 | 0.84 |
| 10-07 00:00 | 124.1 | 84.4 | 0.68 |
| 10-07 01:00 | 166.1 | 90.4 | **0.54** |

**The runs track within 1% for four hours and then split.** That is a barrier, not a fuel
reduction -- a general loss of burnable area would separate them from the first hour. The
fire reached a road around 22:00 and stopped against it; by 01:00 the road has removed
46% of the burned area and the gap is still widening.

### CORRECTION: the first PICO ForeFire result was a defect of mine

An earlier draft of this section reported ForeFire removing **19%** of PICO's burned area
against WRF-SFIRE's 83%, and explained it as the two models handling barriers differently
-- front tracking against a level set. **That explanation was wrong and the number
measured a bug.**

`BARRIER_CODES` was `(92, 94, 95, 96, 97, 98, 99)`. The remapped raster writes old code 91
as **14**, which was not in that tuple, so the `-r max` continuity pass -- the whole
mechanism §5b built to stop barriers being sampled into dashes -- never applied to the new
road barriers. The `near` warp sampled most of them away, the roads came out dotted, and
the fire walked through the gaps. geogrid does not go through `remap_fuel`, so only the
ForeFire side carried the defect: there was no model disagreement to explain.

Fixed in `5c853c1` by adding 14 to `BARRIER_CODES`; 14 never occurs in stock FBFM13, so it
is inert on the unmodified raster.

### What the three configurations actually show

| configuration | orig ha | no_burn ha | removed |
|---|---|---|---|
| WRF-SFIRE (25 m) | 1722.4 | 288.3 | **83%** |
| ForeFire 50 m, broken mask | 2863.2 | 2330.7 | 19% |
| ForeFire 25 m, broken mask | 2909.9 | 381.6 | **87%** |
| ForeFire 50 m, **fixed mask** | 2063.7 | 266.9 | **87%** |

**The models agree.** With the mask working ForeFire removes 87% against WRF-SFIRE's 83%,
landing at 266.9 ha against 288.3 -- within 8% on a fire that burns 1700-2900 ha without
roads. Two models, different front representations and meshes, converging once they see
the same barriers.

**Resolution was not the issue either.** 25 m with the broken mask also recovers 87%,
because finer sampling catches the roads directly without needing the max pass. So the
mask fix buys the right answer at 50 m, which is a quarter of the WindNinja cost: **25 m
is not required for road barriers if the mask is right.**

JH's prediction, made before the runs, was 320 ha for WRF-SFIRE and ForeFire larger but
under 400. WRF-SFIRE came in at **288.3** and fixed ForeFire at **266.9**, both within
about 20%; the ordering claim held for the 25 m run (381.6) but the fixed 50 m run landed
slightly below WRF-SFIRE rather than above.

## 5t. The remap is right except under fire weather

### Against observed outcomes, not model-vs-model

| fire | conditions | observed | original fuels | remapped | correct |
|---|---|---|---|---|---|
| **Pulaski** (Sep) | normal | **3 acres / 1.2 ha** | 25.6 ha WRF, 168 ha FF | **~3 acres WRF, 0 ha FF** | remapped |
| **PICO** | normal | -- | 1722 / 2863 ha | **267 / 288 ha**, models agree | remapped |
| **Palisades** (2025-01-07) | **red flag, Santa Ana** | **9,490 ha** | **5,023.8 ha (0.53x)** | 587.5 ha (**0.06x**) | **original** |

Per JH's rerun, **Pulaski with the remapped fuels came in at about 3 acres against a
reported 3 acres.** That is the fire that was 6.57x too large in ForeFire and 20x too
large in WRF-SFIRE, now essentially exact -- and it is why "ForeFire gives 0 ha" is the
best of the three estimates in absolute error, not an over-correction as recorded earlier
in this session.

**Palisades is the counter-case and it is severe.** The remap removes 88% there -- the
same fraction as PICO's 87%, so the fuel change behaves identically -- but the original
run was already low at 0.53x and the remap makes it **9x worse**. Per JH the no_burn
variant did constrain the fire along roads and through housing developments; the fire
simply did not respect them.

### The rule, per JH

**Use the fuel map with urban areas unburnable, unless the NWS has issued a red flag
warning.** There are plenty of cases where urban development created no barrier at all.

### The switch already exists, and had a bug

`ngfs_incident.NGFS_red_flag()` reads **`nws_fire_wx_code`** from the NGFS detection data
-- per JH, the detection metadata carries a weather code telling whether fire weather
conditions exist -- and sets `self.red_flag` during incident setup, before the job json is
written. Codes **3 or 4** mean fire weather. The file already carries a comment that this
should replace the NWS zone lookup (`red_flag_incident`), which is commented out.

The condition read `(int(c) == 3 or int(c)) == 4`, which for code 3 evaluates `True == 4`
-- **false**. Code 3 never registered; code 4 worked only because `False or 4` is 4.
Fixed in `0099132`. Code 3 is 1.1% of detections in the file sampled against 98.9% for
code 4.

This is operationally the right source: it is per-detection and available **at detection
time**, which §5r established is the binding constraint, unlike the retrospective incident
fields.

### Two cautions before wiring it up

1. **Red flag is a coarse proxy.** Palisades burned through housing on ~100 mph winds and
   single-digit humidity; a red flag warning is issued at far lower thresholds and covers
   large zones for long periods, so it will be true for many fires that behave like
   Pulaski. Too permissive a switch reverts to burnable roads on ordinary fires -- the
   6.57x failure. A wind-speed or spread-component threshold may separate the two better,
   and the gridded winds are in hand at forecast time.
2. **It rests on one red-flag case.** Palisades is among the most extreme WUI fires on
   record. Two or three more extreme-wind comparisons would show whether the rule
   generalises.

### In flight

- JH's WRF-SFIRE `PICO_..._no_burn` run (job 18191), `geo_vars.json_mosaic_noburn91`.
- A ForeFire pair, `wksp_forefire/PICO_fuelcmp_{orig,noburn91}`, winds and moisture held
  identical so the raster is the only variable.

**JH's prediction, recorded before the results:** WRF-SFIRE no_burn91 finishes near
**320 ha**; ForeFire larger but **under 400 ha**.

That is a claim about barrier handling, not just size. WRF-SFIRE's level set can round the
end of a linear barrier; ForeFire's front tracking at `perimeterResolution 25` may not.
§5l found the two models diverge most on fragmented domains, so this is the same question
in a controlled form: identical fuel, identical weather, different front representation.

## 6. Open items

1. **Is blanket code 91 too aggressive?** It makes 51.8% of Pulaski and 29% of PICO
   non-burnable, including low-intensity and open-space development where WUI fires do
   spread. A 30 m "developed" cell is often a 6 m road plus 24 m of grass, and fires cross
   two-lane roads by radiant heat and spotting. **NLCD 23/24 is the ready refinement** --
   30.8% of Pulaski rather than 51.8% -- and the rasters are already grid-aligned.
2. **Rural roads remain invisible.** Neither FBFM13 nor NLCD sees the Dallas track. That
   needs a roads vector layer, not a land-cover raster.
3. **The fragility measure needs re-checking** against the level-set result above before
   it is used to exclude domains.
4. **Agriculture (code 93) is still filled burnable** by the same mechanism, 4.8% of the
   PICO domain. Defensible -- crops burn -- but it is the same code path and was not
   examined.
