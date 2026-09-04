# NGFS automated forecasting

Draft 1 — working document. Facts here were read out of the code on 2026-09-04;
where the code and the intent disagree, this file describes the code.

This package turns NOAA **NGFS** satellite fire detections into WRF-SFIRE
forecast jobs, automatically, on a 30-minute cycle. It assumes you already know
WRF, WRF-SFIRE and wrfxpy; it is a layer on top of that stack, not a
replacement for any part of it.

    NGFS detections (GOES + VIIRS)
              │
        this package  ── decides *whether*, *where* and *when* a fire ignited
              │
        jobs/<grid_code>.json      ← a standard wrfxpy job description
              │
        forecast.sh → src/forecast.py → WRF-SFIRE
              │
        wksp/<grid_code>/          ← standard wrfxpy output
              │
        map_locations.py           ← web map + per-incident pages

The package's whole job is the second box: produce a defensible ignition
estimate and a correct job file. Everything downstream is ordinary wrfxpy.

---

## 1. What "automated" means here

A fire enters the system when NGFS assigns it an incident id (`known_incident_id`,
the IRWIN UUID). Every 30 minutes the system:

1. downloads the last ~48 h of GOES detections and the last ~48 h of VIIRS
   detections, plus a NASA FIRMS URT top-up for detections NGFS has not
   published yet;
2. groups detections by incident id into one `ngfs_incident` per fire;
3. for each unstarted incident younger than `lookback_time` (24 h), estimates
   the ignition point and time;
4. writes a job file, and starts the forecast — once per incident, ever.

"Once per incident, ever" is enforced by `started_inc_ids`, a list of incident
ids carried forward in the saved state. See §6.

### The 4-hour / VIIRS rule

An incident is not forecast the moment it appears. `ngfs_day.start_incidents`
starts it when **either**:

- VIIRS detections exist for it (`len(inc.viirs_data) > 0`), or
- it is more than 4 hours old.

The reason is resolution. GOES gives you a detection within minutes but its
pixel is ~2 km; VIIRS gives you ~375 m but only on an overpass. Waiting for
VIIRS buys a much better ignition location, so the system waits — but not
forever.

While it waits it issues **one** coarse forecast from the GOES/NIFC location,
with the incident id's last field replaced by `InitialFcast`:

    {A5F72521-C511-4FA9-AFDA-C400E719496D}   ← full forecast, VIIRS-refined
    {A5F72521-C511-4FA9-AFDA-InitialFcast}   ← initial forecast, GOES location

Both appear in `jobs/`. The pair is deliberate: it lets you compare a
GOES-only ignition against a VIIRS-refined one for the same fire.

### Unnamed fires

NGFS also reports detections it has *not* named — `Possible Wildland Fire`,
and the near-solar-farm / near-persistent-emitter variants.
`ngfs_day.unknown_incidents` synthesizes an incident id and name for these from
the NWS forecast office code plus the feature tracking id, but only inside the
WFO list in `ngfs_cfg['unknown']['wfo_list']` (currently `KBOU` alone). These
get a shorter, cheaper configuration: `jobs/base_ngfs_cfg_short.json`. An
incident id containing `.` marks one of these (`ngfs_incident.unknown`).

---

## 2. Two systems, and which one issues forecasts

There are two parallel installations sharing one copy of this package.

|                    | monolith (production)      | this package (testing)              |
|--------------------|----------------------------|-------------------------------------|
| entry point        | `src/ngfs_start.py`        | `src/ngfs/ngfs_start_2.py`          |
| runs from          | `/data/jhaley/wrfxpy`      | `/data/jhaley/new_wrfxpy/wrfxpy`    |
| state directory    | `<install>/ngfs`           | `<install>/ngfs`                    |
| issues forecasts?  | **yes**                    | no — submission commented out       |
| detection sources  | GOES                       | GOES + NGFS-processed VIIRS         |

**The monolith is what currently produces real forecasts.** `src/ngfs_start.py`
is a ~1700-line precursor to this package and is not authoritative, but it is
live, and 36 sibling scripts in `src/` import its classes, so it cannot simply
be removed.

**This package deliberately does not submit.** In
`ngfs_incident.start_forecast` the two lines that would launch the job are
commented out:

```python
def start_forecast(self, sleep_time=1):
    cmd = f'sleep {sleep_time}; {self.json_start_code}'
    #subprocess.Popen(cmd,shell=True)
    #os.system(f'(sleep {sleep_time}; {self.json_start_code})')
    print(f'Starting job after {sleep_time} delay: {self.json_start_code}')
    self.started = True
```

Everything else runs: the job file is written, `started_inc_ids` is updated, and
the incident is marked started. So the testing system behaves exactly as if it
had forecast, and its `jobs/` directory is a record of what it *would* have run.
That is the point — it is there to be diffed against the monolith.

Note the consequence before enabling submission: `started_inc_ids` has been
accumulating ids for forecasts that were never issued, including VIIRS-only
fires the monolith never saw. Flipping those two lines on without deciding what
the existing ledger means will cause those fires to be skipped silently.

> **The package directory is a symlink.** There is only one copy:
> `/data/jhaley/new_wrfxpy/wrfxpy/src/ngfs` → `/data/jhaley/wrfxpy/src/ngfs`.
> Editing `src/ngfs` in the main checkout changes what the live loop runs on its
> next execution. Also, `/home/jhaley/work` → `/data/jhaley`, so
> `/home/jhaley/work/wrfxpy` and `/data/jhaley/wrfxpy` are the same directory.

---

## 3. Running it

### Environment

The runtime environment is **`wrf_test`**, not `wrfx`. It is the only conda
environment with sklearn, Basemap and folium together (`wrfx` has no sklearn).
Python 3.7.7, pandas 1.3.5.

If you invoke the interpreter directly instead of activating the environment,
you must set `PROJ_LIB` — Basemap reads it at import and raises
`KeyError: 'PROJ_LIB'` without it:

    PROJ_LIB=/home/jhaley/anaconda3/envs/wrf_test/share/proj

### A normal run

    cd /data/jhaley/new_wrfxpy/wrfxpy
    conda activate wrf_test
    export PYTHONPATH=src
    python /data/jhaley/wrfxpy/src/ngfs/ngfs_start_2.py now unknown

**The working directory matters.** Almost every path in the configuration is
relative — `ngfs_directory` is the bare string `"ngfs"`, job files go to
`jobs/`, logs to `logs/` — all resolved against the process working directory.
Running from the wrong directory writes state into the wrong installation.

`PYTHONPATH=src` is needed because the package imports monolith-era modules
(`utils`, `simple_forecast`, `ngfs_dictionary`, `ngfs_helper`, `state_names`,
`ingest.*`) that live in `src/` rather than in the package.

### Command-line arguments

Arguments are matched by substring against the whole of `sys.argv`, in
`ngfs_day.sys_args_override` and `make_csv_date_str`:

| argument         | effect                                                        |
|------------------|---------------------------------------------------------------|
| `now`            | run for the current time rather than a fixed date; sets `today_forecasts` |
| `ftp`            | select the FTP CSV source                                     |
| `api`            | select the OGC API source                                     |
| a path ending `.csv` | read that file as the GOES source and derive the date from its name; also clears `now` |
| `full_process`   | reprocess incidents restored from state instead of reusing them |
| `behave`, `cawfe` | select the burn model namelist — *currently ineffective, see below* |
| `unknown`        | *(no effect — see below)*                                     |

`unknown` is passed by the cron script but nothing reads it; the unnamed-fire
path is controlled by `ngfs_cfg['unknown']['run']` in the config file.

`behave` and `cawfe` are matched correctly but assign
`ngfs_cfg["fire_namelist_path"]` at the top level, while
`make_incident_configuration` reads `ngfs_cfg['run_cfg']['fire_namelist_path']`.
So the argument is recognized and still has no effect. The burn model is in
practice set by `run_cfg.fire_namelist_path` in the config file. See §9.

### Cron

    0,30 * * * 0-6  cd /data/jhaley/new_wrfxpy/wrfxpy; ./cron_ngfs.sh > log_cron_ngfs.log

`cron_ngfs.sh` activates `wrf_test`, runs the entry point, waits 2 minutes, then
copies the log to `ngfs/cron_ngfs_<timestamp>.log`. A typical run takes about
35 seconds. `log_cron_ngfs.log` is overwritten every cycle, so the timestamped
copies in `ngfs/` are the real history.

The `PYTHONPATH=src` line in `cron_ngfs.sh` has no `export` and so never
reaches python. It works anyway because `ngfs_start_2.py` does its own
`sys.path.insert(1, 'src/')`.

---

## 4. The pipeline

`ngfs_start_2.py` is about 50 lines of orchestration over two classes. Read it
first; it is the outline of everything below.

```
config_manager.load_cfgs()            etc/ngfs.json + etc/conf.json
persistence.get_old_incidents()       restore state from newest pickles
ngfs_day(...)                          ← one object per run
  ├─ add_goes_data()                  NGFS GOES CSV (FTP) or OGC API
  ├─ add_viirs_data()                 NGFS VIIRS scene + FIRMS URT top-up
  ├─ add_red_flags()                  NGFS fire-wx codes, else NWS zones
  ├─ add_pop_data()                   county population table
  ├─ add_incidents()                  detections → ngfs_incident objects
  ├─ process_incidents()              per incident: ignition estimate + job file
  ├─ start_incidents()                apply the start rules, mark started
  └─ save_outputs()                   summary CSV, map PNG, ignition CSV, state
```

### `ngfs_day` — one run

Owns the run's detection data (`self.data` for GOES, `self.viirs_data` for
VIIRS), the incident list, and the started-id ledger. Key detail: `self.data`
is *carried forward* from the previous run's pickle and re-deduplicated on
`(latitude, longitude, acq_date_time)`, so a fire's detection history survives
across cycles rather than being re-downloaded from the start.

`add_incidents` is where old and new meet. For each incident id in the data it
either revives the `ngfs_incident` object from the previous run's state
(appending new detections to it) or constructs a fresh one. An id already in
`started_inc_ids` is constructed fresh but immediately marked `started`, which
is what prevents a second forecast for a fire that dropped out of the data and
came back.

VIIRS-only incidents — ids present in the VIIRS data but not in GOES — are
picked up here too. These are fires the GOES-driven monolith cannot see.

### `ngfs_incident` — one fire

`process_incident` is the scientific core:

- `find_old_features` walks backward through the full GOES data for earlier
  occurrences of the incident's `feature_tracking_id`, from before NGFS named
  it. This is what recovers detections predating the incident id.
- `set_ignition_point` estimates location and time. It averages detections from
  the first 3 hours, preferring the **SWIR (band 5) terrain-corrected**
  positions (`latitude_b5`/`longitude_b5`, with `-999` as fill) over the nominal
  pixel centers — but only if the two estimates agree to within 0.04°,
  otherwise it falls back to the nominal mean. Ignition time is the earliest of
  `acq_date_time`, `incident_start_time` and `pixel_date_time`.
- `polar_ign_estimate` refines that using VIIRS: if a VIIRS detection falls
  inside the GOES ignition pixel, `new_ign_latlon` / `new_ign_utc` are set and
  take priority in the job configuration.
- `incident_landcover`, `set_incident_demographics` and `NGFS_red_flag` attach
  fuel, county/state/population, and fire-weather context.

`make_incident_configuration` deep-copies the base job configuration and
specializes it: ignition point and time, domain center, `truelats`/`stand_lon`,
start and end times, GRIB source and Landfire tables by region (AK → NAM198,
HI → NAM196, PR/VI → GFSF, western states → the 2024 Landfire vtable), the
GRIB cycle start for cycle-based sources, and FMDA versus equilibrium fuel
moisture. Then `expand_domain_size` grows the domain if the detection footprint
would not fit inside it, and the result is written to
`jobs/<grid_code>.json`.

`grid_code` is `<incident_name>_<start_utc>_<uuid without braces>` with `#()  :`
replaced by underscores.

---

## 5. Configuration — `etc/ngfs.json`

Read by `config_manager.load_cfgs`, which also loads the wrfxpy configuration
named by `wrfxpy_cfg` (normally `etc/conf.json`). The base job description comes
from `jobs/base_ngfs_cfg.json` via `make_base_configuration`; if it is missing,
`simple_forecast.questionnaire()` runs interactively — which will hang a cron
job, so make sure it exists.

| key | meaning |
|---|---|
| `ngfs_directory` | state and output directory, **relative to cwd** (`"ngfs"`) |
| `base_cfg` | base wrfxpy job description all forecasts are derived from |
| `pop_data` | county population table (tab-separated) |
| `run_cfg.fire_namelist_path` | selects the burn model; `behave_13` also forces `time_step = 6` |
| `run_cfg.lookback_time` | hours; incidents older than this are never started (24) |
| `run_cfg.num_starts` | maximum forecasts per run (50) |
| `run_cfg.job_sleep` | seconds staggered between job starts, to avoid metgrid collisions (150) |
| `run_cfg.today_forecasts` | equivalent to passing `now` |
| `goes_cfg` / `viirs_cfg` | source (`ftp` or `api`), host, remote and ingest directories, `days_to_get`, satellites/sectors |
| `firms_cfg` | NASA FIRMS URT satellites, ingest directory, `days_to_get`, API token path |
| `fmda_cfg.use_fmda` | fuel-moisture data assimilation; CONUS only — AK/HI fall back to equilibrium FMC |
| `region_cfg` | per-region GRIB source and Landfire vtable overrides |
| `unknown` | unnamed-fire handling: `run`, `wfo_list`, and the short job configuration |
| `perims_cfg.perim_dir` | NIFC perimeter directory for `add_nifc_perims` |

`region_cfg` holds only the regions that need their own fuels, elevation **and**
GRIB source: `alaska`, `hawaii`, `prvi`. CONUS is deliberately absent — see §5a.
Each entry needs `state`, `geo_vars_path`, `grib_source` and `msg`, and a state
must appear in at most one entry, since the last match would otherwise win.

One config caveat: `perims_cfg.perim_dir` is `ngfs/perims/`, which does not
exist in the testing install, and `add_nifc_perims` has no callers.

`etc/ngfs.json` is per-installation and gitignored. The tracked example is
`etc/ngfs.json.initial`, following the wrfxpy convention (`conf.json.initial`,
`tokens.json.initial`); copy it and edit for your system.

---

## 5a. Fuels

Fuel category (`NFUEL_CAT`) and elevation (`ZSF`) come from a small JSON naming
two GeoTIFFs, resolved in this order by `forecast.py`:

1. `geo_vars_path` in the job configuration, which the NGFS base configs set to
   `etc/vtables/geo_vars.json_mosaic`, overridden per region by `region_cfg`;
2. failing that, `etc/vtables/geo_vars.json`.

**CONUS fuels are not chosen by state.** LANDFIRE refreshes its products in a
rolling fashion, so a release covers only the areas updated that cycle and the
rest of the raster is fill. `LF2025_FBFM13_CONUS` carries valid codes over
roughly the western third of CONUS and nothing elsewhere, despite its name.
Those update regions cut across state lines — New Mexico is covered in its
western third only — so no state list can select the right release, and an
earlier attempt to do so would have handed fill-value fuel to fires in
uncovered areas.

Instead, `src/ingest/landfire_mosaic.py` stacks the staged releases
newest-first and fills each pixel from the newest release that has data there.
The current CONUS mosaic takes 21.4 % of its pixels from LF2025 and 40.4 % from
LF2024, giving full CONUS land coverage with the newest available fuels
everywhere.

Two fill values matter and both are treated as gaps: `32767`, the declared
nodata, meaning *inside the footprint but not updated this cycle*; and
`-9999`, which is **not** declared in the GeoTIFF metadata and means *outside
the product footprint*. Because `-9999` is undeclared, any gap test based on a
raster's own nodata value silently accepts it as a fuel category. The mosaic
collapses both to a single declared nodata so nothing downstream inherits that
trap.

Each mosaic is written with a `.manifest.json` sidecar recording the exact
releases, their file sizes and mtimes, and the pixel count each contributed —
so a forecast's fuels are traceable to specific LANDFIRE releases rather than
to "some mosaic".

**When a new release is staged:**

    export PYTHONPATH=src
    python -m ingest.landfire_mosaic --list
    python -m ingest.landfire_mosaic --product FBFM13 --region CONUS
    python -m ingest.landfire_mosaic --product FBFM13 --region CONUS --build

The third command is a dry run reporting each layer's estimated contribution;
`--build` writes the raster and its manifest. Expect roughly half an hour for
CONUS. No code or state lists need to change. The tool refuses to proceed if
the layers disagree on CRS, pixel size or grid alignment, rather than
resampling — pasting misaligned layers would put fuel categories on the wrong
ground.

Alaska, Hawaii and PRVI are single releases rather than mosaics: AK is a
complete raster (LF2025, verified 100 % valid), HI uses LF2023, and PRVI has
not been updated since the 2020 release.

---

## 6. State and the incident ledger

State lives in the `ngfs_directory` as one pickle of the whole `ngfs_day`
object per run:

    ngfs/pkl_ngfs_day_<YYYY-MM-DD>_<YYYYMMDD_HH_MM>_testing.pkl.gz

`persistence.get_old_incidents` globs that directory, sorts by mtime, keeps the
last 7 days, and merges forward: the union of all `started_inc_ids`, the
incident objects, and the newest run's detection frame. Facts that are easy to
get wrong:

- **`save_pickle` only runs when `start_count > 0`** (see
  `ngfs_day.save_outputs`). State is written only by runs that started a
  forecast, so pickle timestamps are *not* a log of when the loop ran.
- **mtime is load-bearing.** Both the age filter and the "newest state" choice
  sort by mtime. Any tool that rewrites a state file must preserve its
  timestamp, or it promotes stale state to newest and corrupts the ledger.
- **Pickles are gzip level 6**, written to a pid-qualified temp name and
  `os.replace`d into position, so an interrupted save cannot leave a truncated
  file where the next run looks. Reads accept `.pkl`, `.pkl.gz` and `.pkl.xz`,
  so older uncompressed files stay readable.
- **Class paths differ by writer.** The monolith runs as a script, so its
  pickles record `__main__.ngfs_day`; package-written ones record
  `ngfs.ngfs_day.ngfs_day`. Reading a monolith pickle from the package raises
  `AttributeError: module '__main__' has no attribute 'ngfs_day'`.
  `persistence.pickle_writer` determines ownership from the first 128 bytes —
  do not rely on the `_testing` filename suffix.
- Only the 48 most recent hours of detections are kept in each pickle
  (`save_pickle` trims before writing). Since a pickle is written every 30
  minutes, each detection is stored roughly 96 times across the series.
- `started_inc_ids` grows without bound — currently ~3,790 entries, each new
  incident id compared against the whole list. It lives inside the pickle, so
  there is nowhere to record a first-seen timestamp; bounding it properly means
  extracting an explicit ledger file. Open.

### Compressing old state

`persistence` has a CLI for the backlog. It gzips at the byte level, verifies by
SHA-256 of the decompressed bytes, and deletes the original only on a match.
Without `--delete` it only reports. Files newer than `--min-age-hours`
(default 6) are never touched, so it is safe to run while the loop is live, and
an interrupted sweep can simply be re-run.

    cd /data/jhaley/new_wrfxpy/wrfxpy
    PROJ_LIB=/home/jhaley/anaconda3/envs/wrf_test/share/proj \
    PYTHONPATH=/data/jhaley/wrfxpy/src \
    /home/jhaley/anaconda3/envs/wrf_test/bin/python -m ngfs.persistence ngfs [--delete] [--limit N]

It **refuses the monolith's own state directory**: `src/ngfs_start.py` globs
only `*.pkl`, so compressing there would hide state from it and cause fires to
be forecast a second time.

---

## 7. Modules

### The pipeline

| module | role |
|---|---|
| `ngfs_start_2.py` | entry point; ~50 lines of orchestration |
| `ngfs_day.py` | one run: data acquisition, incident assembly, start rules |
| `ngfs_incident.py` | one fire: ignition estimation, job configuration |
| `persistence.py` | state save/restore, pickle compression, run outputs |
| `config_manager.py` | loads `etc/ngfs.json` and the base job configuration |
| `constants.py` | thresholds, column groups, region tables |

### Detection ingest

| module | role |
|---|---|
| `ngfs_ftp.py` | NGFS CSV over FTP/HTTP from `bin.ssec.wisc.edu`; the default source |
| `ngfs_api.py` | the same detections over NGFS's OGC API; selected with `api` |
| `firms_data.py` | NASA FIRMS URT/NRT VIIRS, for detections NGFS has not published yet |
| `ngfs_dictionary.py` | column dtypes and the v1↔v2 NGFS schema translation |

The v1/v2 split matters when reading old data: v1 used `lat`/`lon`/`incident_id_string`,
v2 uses `latitude`/`longitude`/`known_incident_id`. `ngfs_day.data_to_v2`
detects v1 by the presence of a `lon_tc` column and renames through the
dictionary. Much of the code still carries `try`/`except KeyError` pairs for
both spellings.

### Context data

| module | role |
|---|---|
| `nifc.py` | current NIFC incident locations and perimeters (ArcGIS) |
| `perim_cache.py` | caches NIFC perimeters, storing only geometry changes |
| `watch_duty.py` | scrapes watchduty.org reports for cross-reference |

### Visualization and web

These are separate consumers. They read `wksp/` output and live detection
feeds; they are not part of the forecasting pipeline.

| module | role |
|---|---|
| `map_locations.py` | the main folium web map; run from cron every 10 minutes |
| `incident_webpage.py` | per-incident page: forecast perimeters, area series |
| `map_utils.py` | shared map building blocks, popups, legends, `scp` upload |
| `sat_webpage.py` | JPSS satellite ground tracks and overpass prediction |

`map_locations.py` runs via `ngfs_map.sh` from `/data/jhaley/wrfxpy` (the
monolith install), writes `wildfire_map.html`, and `scp`s it to
`~/NGFS/` on the CSU engineering servers, trying each `linux*.engr.colostate.edu`
in turn.

### Tools

| module | role |
|---|---|
| `latency.py` | job-file time versus estimated ignition time |
| `upload_forecasts.py` | pushes forecast directories to a remote host |
| `watchduty_nifc_table.py` | maps IRWIN ids to watchduty catalog numbers |
| `utilities.py` | one timestamp helper *(broken — see §9)* |

Outside the package but part of this workflow:

| module | role |
|---|---|
| `src/ingest/landfire_mosaic.py` | builds the LANDFIRE fuels mosaic and its manifest (§5a) |

Deliberately untracked in git: `latency.py`, `upload_forecasts.py`,
`watchduty_nifc_table.py`, `utilities.py` (tooling), `tmp.txt` (a stray diff
dump), `ngfs_api.json` (a downloaded OpenAPI spec nothing reads).

---

## 8. Outputs

Per run, in `ngfs_directory` unless noted:

| file | written by | content |
|---|---|---|
| `jobs/<grid_code>.json` | `make_incident_configuration` | the wrfxpy job description |
| `logs/<grid_code>.log` | `forecast.sh`, once submission is enabled | forecast log |
| `detection_summary_testing.csv` | `detection_summary` | per-satellite counts and total FRP, latest run |
| `detection_summary_<date>_testing.csv` | `detection_summary` | same, kept per day |
| `NGFS__testing<date>.png` | `print_base_map` | CONUS map of incidents by status; Alaska panel spliced on when a fire is north of 54° |
| `forecast_ignition_pixels_<date>_<HH_MM>.csv` | `save_incident_text` | the ignition estimate for each newly started incident, including whether VIIRS refined it |
| `pkl_ngfs_day_..._testing.pkl.gz` | `save_pickle` | run state; **only when a forecast started** |
| `cron_ngfs_<timestamp>.log` | `cron_ngfs.sh` | that cycle's full stdout |

`forecast_ignition_pixels_*.csv` is the file to look at when checking ignition
quality: it carries both the original GOES estimate and the VIIRS-refined
`forecast_ign_lat`/`forecast_ign_lon`/`forecast_ign_UTC`, plus the
`viirs_pixel_ign` flag.

---

## 9. Rough edges

Recorded so nobody documents them as working, or debugs them twice. None of
these are currently causing failures.

**Dead code in the package** (verified zero callers):

- `ngfs_day.cluster_data` — would raise `AttributeError` if called;
  `db.labels_c_cloc` should be `db.labels_`. Also uses `DataFrame.append`,
  removed in later pandas.
- `ngfs_day.add_polar_data` — superseded by `add_viirs_data`; commented out in
  the entry point.
- `ngfs_day.prioritize_incidents` — population-based ordering, commented out in
  the entry point.
- `ngfs_incident.add_nifc_perims` — no callers, and its configured directory
  does not exist here.
- `ngfs_incident.find_old_detections` — superseded by `find_old_features`.
- `set_cmd_str` is defined *inside* `process_incident`, so it is a local
  function that is never called, not a method.

**Bugs that silently do nothing:**

- `sys_args_override`: the `behave`/`cawfe` branches write
  `ngfs_cfg["fire_namelist_path"]` at the top level, while
  `make_incident_configuration` reads `ngfs_cfg['run_cfg']['fire_namelist_path']`.
  The arguments are matched but have no effect. The
  `==`-for-`=` and `is`-for-`in` typos in this method were fixed on 2026-09-04.
- `make_incident_configuration` used to test for `region_cfg_REMOVE_THIS`, so
  the package never consulted `region_cfg` and fell through to a hardcoded
  state list. Fixed on 2026-09-04: fuels now come from the mosaic (§5a) and
  `region_cfg` is read for real. Note that the **monolith always did** read
  `region_cfg` — its guard is `'region_cfg' in ngfs_cfg.keys()` — so this was
  only ever a package defect, and production was not mis-selecting fuels.
- `unknown_incidents` reports `unknown_count - skipped`, but `skipped` is reset
  inside the per-WFO loop, so with more than one WFO the count is wrong.
- `utilities.timestamp_from_string` takes `date_str` and uses `csv_date_str` —
  an immediate `NameError`. The working copy is `ngfs_day.timestamp_from_string`.
- `ngfs_day.py`'s `__main__` block calls `ngfs_day(ngfs_cfg, sys_args=...)` and
  `nd.add_data(...)`; neither the keyword nor the method exists. The module's
  self-test does not run.

**Deferred, from the current state of the installation:**

- ~201 GB of monolith-format pickles (about 4,201 files) sit in the testing
  install's state directory. Nothing appears to read them; they are excluded
  from `get_old_incidents` only by the 7-day age filter.
- Two pickle files match neither class marker. Inspect before assuming junk.

---

## 10. Troubleshooting

**`KeyError: 'PROJ_LIB'`** — Basemap, at import. Activate `wrf_test` or set
`PROJ_LIB` (§3).

**`ModuleNotFoundError` for `utils`, `simple_forecast`, `ngfs_helper`,
`state_names`** — `PYTHONPATH` does not include `src/`. These are monolith-era
modules outside the package.

**`AttributeError: module '__main__' has no attribute 'ngfs_day'`** — you are
reading a monolith-written pickle from the package. Expected; since the
`get_old_incidents` hardening it is logged and skipped rather than fatal.

**No state pickle from a run** — normal if no forecast started; `save_pickle`
is gated on `start_count > 0`.

**A run hangs with no output** — `jobs/base_ngfs_cfg.json` is probably missing,
so `make_base_configuration` fell through to the interactive questionnaire.

**A fire is not being forecast** — check, in order: is its id already in
`started_inc_ids`; is its estimated ignition older than `lookback_time` (24 h);
has `num_starts` (50) been reached this run; is it younger than 4 hours with no
VIIRS detections yet (in which case only the `InitialFcast` job is expected).

**Nothing at all happened** — remember submission is commented out in this
package (§2). Look in `jobs/` for what it *would* have run.

**Fuels look wrong, or geogrid finds no NFUEL_CAT** — check what the job file's
`geo_vars_path` resolves to *from the working directory the forecast runs in*;
these paths are relative and the two installations have separate
`etc/vtables`. Then confirm the GeoTIFFs it names exist, and that the ignition
point actually falls inside the fuels raster with a value that is neither
`32767` nor `-9999` (§5a). A CONUS raster does not cover PRVI at all.

---

## 11. Git

Work happens on `james_ngfs`; `master` has no `src/ngfs` at all. A new worktree
defaults to branching from `origin/master`, so create them explicitly:

    git worktree add -b <name> .claude/worktrees/<dir> james_ngfs

A fresh worktree also lacks the ~74 untracked `src/*.py` modules, the
per-installation configs (`.gitignore` names `etc/conf.json`,
`etc/tokens.json`, `etc/clusters.json` and `etc/ngfs.json` individually — there
is no `etc/*.json` wildcard), and the concrete `etc/vtables/geo_vars.json_*`
files, which are untracked rather than ignored. To make a worktree importable
for testing, symlink those in from the main checkout, and remove the symlinks
before committing.

`core.fileMode` is set to `false` in this repo; without it, 181 files show as
modified purely from a past `chmod`.

The main checkout holds roughly 12,000 untracked files and about 1 TB of data,
including a 45 GB pickle and multi-GB LANDFIRE archives, and `.gitignore` does
not cover `landfire/`, `*.pkl`, `*.tif`, `*.geojson` or the root data dumps.
**Never `git add -A` here.**
