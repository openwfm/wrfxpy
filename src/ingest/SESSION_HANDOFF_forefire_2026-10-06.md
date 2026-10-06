# ForeFire session handoff — 2026-10-06

The perturbed ensemble, built and made fast; one real bug behind it; and the finding that
everything validated so far comes from a population we are not trying to forecast.

**Section numbers continue from the 2026-10-05 handoff on purpose**, as that one continues
2026-10-01: §5q and §5r are referenced by number from commit messages. Read this as the
10-06 continuation, not as a standalone.

## 1. TL;DR

- **A 49-member ensemble now costs about 8 minutes**, against ~18.5 for the WindNinja
  pass it shares. 49 of 49 members complete. Perturbs pSAF and fuel moisture; the
  moisture axis had to be built, since moisture is not a `setParameter` (§5q).
- **One real bug** (`b8e89c9`): sweep members were handed a fuel table path the container
  could not resolve, so ForeFire ran with **no fuel properties** and 42 of 49 fires died
  after two steps. Both moisture directions failing identically is what gave it away
  (§5q).
- **The container was half the cost of every run.** `squashfuse` is absent, so apptainer
  extracted a 298 MB SIF on **every** invocation. A pre-built sandbox takes a 28-step
  member from 102 s to 64 s, identical result, and is now live in `etc/forefire.json`
  (§5q).
- **`sweep_job` takes `workers`.** Threads, not processes, because each step is already an
  `apptainer exec`. Measured 1/4/8/12/16 with 8/8 complete and zero failures at every
  level (§5q).
- **Both ensembles scored against IR, and neither perimeter is a valid target.** Dallas's
  observation sits below all 49 members because a road the fuel map lacks stopped the real
  fire; ORCHARD_ROAD's is reached only at the 98th percentile because crews burned
  additional ground. 48% and 71% of those boundaries are straight runs (§5q).
- **0 of 34 forward-run fires came from an unassociated hotspot** -- none of 10,345
  workspaces did. `unknown.run` is `False` in `etc/ngfs.json`. Everything calibrated so
  far is a managed fire (§5r).

## 5q. The ensemble work: one real bug, and the container is the bottleneck (2026-10-06)

Per JH: a 49-member ensemble, 7 pSAF x 7 fuel-moisture offsets, scored against an
observed perimeter. The ensembles themselves still need re-running; what came out of the
attempt is more useful than the attempt.

### The moisture axis had to be built, and it is the stronger one

Moisture is not a `setParameter` -- it reaches ForeFire through the per-step
`fuelsTableFile` -- so `md_series` and `make_script_set` gained `md_offset`, applied to
every dead class **before** the range check so a dry member clamps at the 0.02 floor
rather than going negative. In `sweep_job` it travels as the reserved key `_md_offset`.

Measured on Dallas, over the ranges chosen:

| axis | range | final area |
|---|---|---|
| moisture | -0.040 to +0.040 | 4031.3 -> 388.0 ha, **10.4x** |
| pSAF | 0.221 to 0.513 | 125.4 -> 626.9 ha, **5.0x** |

**Moisture has twice the leverage of pSAF**, which inverts the weighting a 7x7 grid
implies. §5h argued for moisture because it reaches moisture-limited fires that pSAF
cannot; it turns out to be the dominant axis generally.

### The bug: sweep members got a fuel table path the container could not resolve

42 of 49 members died after two steps, in **both** moisture directions identically, with
no crash. A drier fire burning less is not physics, and that symmetry is what gave it
away. ForeFire's own log had it:

```
WARNING: could not load file for fuel properties /work/fueltables/fuels_Md_cawfe_0.0480_...csv
```

It ran with **no fuel properties at all** and the fire went out.

A sweep runs each member in `base/<tag>` while binding `base` as `/work`, because the
netcdf symlinks point back into base and must be followable. `make_script_set` resolved
the fuel table against the *member* directory, writing `/work/fueltables/<name>` for a
file at `/work/<tag>/fueltables/<name>`.

**Why it hid:** a sweep that varies only `setParameter` values leaves Md untouched, so
the filenames match tables the baseline run already left in `base/fueltables`, and the
wrong path happens to resolve. The zero-offset members were not working -- they were
silently reading the **baseline's** tables. The first sweep ever to perturb moisture
changed the filenames and everything failed. **§5f's pSAF sweep is unaffected** for the
same reason it hid.

Fixed in `b8e89c9`: `make_script_set` takes `bind_dir` and resolves against it, falling
back to `forefire_dir` so ordinary runs are unchanged. Verified -- offsets 0.000, +0.001,
-0.040, +0.040 all complete 28 steps, giving 1997.9 / 1963.8 / 4031.3 / 388.0 ha,
monotonic in moisture.

**This is also why a 0.001 offset was the right test.** A thousandth of a point of
moisture cannot extinguish a fire, so when it failed exactly like -0.040 the fault had to
be the code path rather than the perturbation.

### Member concurrency: `sweep_job` now takes `workers`

Threads, not processes: every ForeFire step is already an `apptainer exec` subprocess and
Python releases the GIL for its duration, so N threads put N containers in flight.
Checked rather than assumed -- no `os.chdir` in the path, every write target is the
per-member `run_dir`, and `_FMC_GW`/`_FMDA_GRID` are check-then-set with fully built
values.

8 members per level, pSAF only so `workers` is the single variable:

| workers | wall min | min/member | complete | speedup |
|---|---|---|---|---|
| 1 | 14.8 | 1.86 | 8/8 | 1.0x |
| 4 | 8.2 | 1.02 | 8/8 | 1.8x |
| 8 | 8.4 | 1.04 | 8/8 | 1.8x |
| 12 | 5.8 | 0.72 | 8/8 | 2.6x |
| 16 | 3.9 | 0.49 | 8/8 | 3.8x |

**Zero container failures at every level.** The 14 `exit code 255` failures seen earlier
at 25-way were on the broken code -- members dying for the fuel-table reason -- and were
wrongly attributed here to concurrency. 16x the workers buys only 3.8x, and 4 and 8 are
indistinguishable, which is a shared-resource plateau rather than CPU saturation.

### The shared resource is the container image, and it is cheap to fix

`squashfuse` is not installed, so apptainer cannot loopback-mount the 298 MB SIF and
**extracts it to a temporary sandbox on every invocation**:

| | SIF | pre-built sandbox |
|---|---|---|
| one `exec` | 1.45 s | **0.18 s** |
| 8 at once | 5.4 s | **0.2 s** |
| one real 28-step member | 102.1 s | **63.5 s** (identical result, 2564.2 ha) |

Building the sandbox once takes **0.89 s and 871 MB**. The 8-at-once figure is the
important one: the SIF path serialises, which is exactly the 4-to-8 plateau above.

**What it takes:** `apptainer build --sandbox <dir> forefire.sif`, then point
`apptainer.sif_image` at the directory. No code change -- `exec` takes a sandbox wherever
it takes a SIF. Better still would be installing `squashfuse`, which needs admin but
needs no duplicate on disk and fixes every apptainer user on the machine.

At 28 steps a member and 49 members an ensemble that is **1,372 extractions of a 298 MB
image**, about 30 minutes of pure overhead per ensemble. It also helps every other run in
the project -- the cron, the forward batch, the single deterministic forecast -- which
matters for the 8-hour latency target.

`etc/forefire.json` has **not** been changed; the sandbox built for the test is in the
session scratchpad and is temporary.

### Perimeter quality: trust vertex spacing, not metadata

Per JH, `poly_MapMethod` is unreliable; dense vertices mean IR, sparse means hand-drawn.
Measured over the five perimeters matching forward runs:

| fire | verts | median gap | metadata says | vertex test |
|---|---|---|---|---|
| ORCHARD_ROAD | 877 | **4.0 m** | Auto-generated | **IR** |
| Dallas | 173 | **13.9 m** | Mixed Methods | **IR** |
| Bowie_9066 | 14 | 28.0 m | Mixed Methods | too small to judge |
| Ice_Box | 28 | 114.2 m | Auto-generated | sparse (and `GISAcres` contradicts its own polygon) |
| Hatch_Grande | 226 | **164.8 m** | **IR Image Interpretation** | **hand-drawn** |

**The metadata is exactly backwards on the two extremes.** Hatch_Grande, the only one
labelled IR, is the sparsest; ORCHARD_ROAD, labelled auto-generated, is the densest by
40x. Hatch_Grande is unusable anyway -- 3650 ha two hours into a forecast started from a
point ignition, so it measures a pre-existing fire, not spread.

### Dallas: the missing road, and what an ensemble cannot fix

Per JH, a road appears to have stopped northward spread, the direction most growth took.
The perimeter carries a **323 m run straight to within 3.6 m**, north of the centroid --
a linear feature, not fire behaviour. Sampled in ForeFire's own fuel field, **100% of it
is burnable** (fuel 2 and 5), and the domain is only 0.9% non-burnable overall. That is
§4's "highways are deleted on purpose": `remap_fuel` fills LANDFIRE code 91.

So part of Dallas's 6x overprediction is **geometric**, and no pSAF or moisture member
can reproduce a boundary the fuel map does not contain. Worth keeping as the clearest
statement of what a parameter ensemble can and cannot diagnose.

## 5r. Everything validated so far is a known incident (2026-10-06)

Checked across the forward set and the whole archive: **0 of 34 forward-run fires, and
none of 10,345 workspaces, came from an unassociated hotspot.** Every one carries an
IRWIN GUID.

The reason is a switch, not an accident. `etc/ngfs.json`:

```
unknown:  run False        wfo_list ['KBOI']        red_flags False
          unknown_cfg_file jobs/base_ngfs_cfg_short.json   (exists)
```

`ngfs_start.py:1539` enters the unknown-fire branch only on `'unknown'` in argv or
`ngfs_cfg['unknown']['run']`. It is off, and scoped to one office, so the nightly cron
has only ever processed fires that already had an incident record.

**Correction to a test used here:** a `grid_code` is never missing an id string -- per JH,
an artificial one is constructed from the `feature_tracking_id` for hotspots with no
incident. So "has a GUID" does not by itself prove a fire was known; the config above is
the evidence that stands.

### Why this matters more than it sounds

A fire with an incident record is a fire with crews, control lines and burnout
operations. That is the same fact that showed up geometrically in §5q:

| | straight runs >100 m, as a share of boundary |
|---|---|
| Dallas | **48%** |
| ORCHARD_ROAD | **71%** |

Free-burning fire does not make straight lines. So the validation targets are compromised
**because of the selection**, not by bad luck -- and the two fires are biased in opposite
directions, Dallas's road holding the fire below free spread while ORCHARD_ROAD's burnout
adds area the fire never spread into. Neither is a test of spread physics.

**Consequence for §5p.** The 31-fire forward baseline -- geo mean 1.42, 58% within 2x --
is calibrated on managed fires against WRF-SFIRE. As model-vs-model that is sound. Its
transfer to unassociated hotspots is **untested**, and those are typically smaller and
earlier, where the `area^-0.198` size dependence bites hardest.

### The target population, and why the switch stays off for now

GOES-19 alone, last 24 h: **22,080 detections, 18 known incidents, 13,359 possible
wildland fires.** Per JH it is a good bet **99.99% of those are not wildland fires at
all**, and deciding which deserve a forecast is an open problem. NGFS already accounts for
FRP and signal persistence, so the useful additions are elsewhere:

- proximity to **populated areas**
- hotspots during **fire weather watches and warnings**
- **lightning detection** data, when available, to indicate probable ignition sources

Per JH the sequence is deliberate: **get the low-latency pipeline working, then turn the
switch on.** Turning it on first would only queue forecasts for things that are not fires.

That also explains why the 8-hour product matters for exactly the fires the switch would
admit -- new hotspots, no incident, no crews yet. For a known incident with resources
already dispatched, a forecast is far less decisive. ORCHARD_ROAD is the illustration:
resources were dispatched at 22:14 UTC, **32 minutes before** the satellite-estimated
ignition of 22:46.

### Also recorded: do not use retrospective incident fields operationally

`attr_FireDiscoveryDateTime` (100% populated) and `attr_DiscoveryAcres` (87%, median
0.40 ha, 38% already over 1 ha) are tempting as a better ignition time and a non-point
initial state. **They are not available for the use case** -- they exist only once a fire
has an incident record, which is exactly what the target hotspots lack. The operational
analogues at detection time are `frp`, `feature_frp`, `pixel_area`, `confidence` and
`feature_tracking_id`. Noted because the fields look useful and are not.

Per JH, newer NGFS versions carry a **platform-independent tracking id**, so GOES and
VIIRS detections will share one reference per incident -- which is what would let a
feature be followed across satellites rather than per-platform.

## 6. Open items

Carried from the 10-01 and 10-05 handoffs, plus:

1. **Two Dallas members did not complete** (47 of 49). Unchecked -- it may be two fires
   running out of the domain rather than a failure.
2. **Stale diagnostic directories** in `wksp_forefire/Dallas_*`: `d_*`, `f_*`, `w01_*`
   through `w16_*`, `img_*`, `live_check`. Harmless, and excluded from scoring by the
   `p*_m*` glob, but worth clearing.
3. **The sandbox goes stale silently** if `forefire.sif` is ever rebuilt. Nothing links
   them. Rebuild with `apptainer build --sandbox /home/jhaley/forefire/forefire_sandbox
   /home/jhaley/forefire/forefire.sif`, or install `squashfuse` and drop the duplicate.
4. **Concurrency above 16 is untested.** The trend was still improving there, and the
   earlier 25-way failures were the fuel-table bug rather than contention, so the real
   ceiling is unknown.

## 7. NEXT SESSION

**The validation route is closed; the latency route is open.** §5r is the reason to stop
looking for a perimeter to score against. Any fire with a usable perimeter had crews on
it, and the fires the pipeline exists for mostly never produce one. Model-vs-model stays
the primary calibration, not as a fallback but because the target population has no
observations by construction.

So the work that moves things forward:

**1. Latency, end to end.** §5q removed the container overhead and made members cheap;
WindNinja is now essentially the whole cost of a forecast (18.5 min of ~20). An 8-hour
forecast is ~6 min of WindNinja and under a minute of ForeFire. Driving that down is the
single thing that unblocks the unknown-fire switch, per JH.

**2. Triage, which has no answer yet.** 13,359 possible wildland fires a day on GOES-19,
~99.99% not fires. FRP and persistence are already in NGFS; the candidates JH named are
proximity to populated areas, fire weather watches and warnings, and lightning detection
when available. The platform-independent tracking id in newer NGFS versions is what would
let a feature be followed across GOES and VIIRS rather than per-platform.

**3. The ensemble is ready and has nothing worth running on.** It works, it is fast, and
both axes behave. What it needs is a fire where the answer is not predetermined by
suppression -- which, per §5r, means the switch, which means latency first.

Still open and untouched: isolating the mesh (09-29 §7), Finley_RX at 0.06 as evidence
that pSAF and moisture are not independent knobs (§5p), and excluding the 73
connectivity-flagged domains before any per-fire pSAF fit (10-01 §6).
