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
