# NGFS session handoff — 2026-09-04

Follows `SESSION_HANDOFF_2026-09-03.md`, which still holds for system topology,
the `wrf_test`/`PROJ_LIB` gotchas and the state-pickle work. Read
`/data/jhaley/wrfxpy/AGENTS.md` first for the approval gate and commit policy.

Today's session was documentation plus a fuels correctness fix. Most of what
was learned is now in `src/ngfs/README.md`, which is tracked — this file
records what happened, what is verified, and what is still open.

---

## 1. Commits — five, all on `james_ngfs`, all UNPUSHED

    325ec0b  Fix comparison and identity operators in sys_args_override
    996e35d  Select fuels by data coverage, not by state
    ee864cd  Add landfire_mosaic: build a fuels mosaic from staged releases
    a381c53  Add ngfs.json.initial example config and ignore the live one
    d55a931  Document the NGFS package: README and module docstrings

`origin/james_ngfs` is still at `9e35dda`, so **b70ca1e plus these five are
unpushed** — six commits total. Pushing has to happen from a separate terminal;
it does not work from inside Claude Code (no ssh-agent, passphrase-encrypted
key, no askpass). This was raised and declined before; do not propose changing
the ssh setup.

`325ec0b` is the user's own edit, committed separately at his request. Both it
and `996e35d` shared a file with the documentation commit, so `ngfs_day.py` and
`ngfs_incident.py` were each split by hand: save the file, `git checkout HEAD --`
it, re-apply only the functional hunk, `git commit -- <path>` (pathspec, which
bypasses the index and leaves other staged paths alone), then restore. Verified
by md5 both times. **Note `git commit -- <path> -m "msg"` does not work** —
everything after `--` is a pathspec; put `-m`/`-F` before it.

---

## 2. The fuels problem and what was done

**The problem.** LANDFIRE refreshes its products by region, so a release covers
only the areas updated that cycle and the rest of the raster is fill.
`LF2025_FBFM13_CONUS` has valid fuel codes over roughly the western third of
CONUS and nothing elsewhere, despite the name and a CONUS-wide bounding box:

| product | valid |
|---|---|
| LF2022_FBFM13_220_CONUS | 61.8 % |
| LF2023_FBFM13_240_CONUS | 61.7 % |
| LF2023_FBFM13_240_CONUS_sw_update | 11.3 % |
| LF2024_FBFM13_250_CONUS | 61.7 % |
| LF2025_FBFM13_CONUS | **21.9 %** |

Coverage does not follow state lines — New Mexico is covered in its western
third only. Measured against 400 recent ignition points from job files: LF2024
has valid fuel at 400/400, LF2025 at 269/400. Of the 132 fires the old code
tagged as "updated states", six would have received fill-value fuel from
LF2025, all in eastern and central New Mexico. So simply repointing at `_2025`
was not safe, and **no state list can be correct** for this.

**Two fill values, one undeclared.** `32767` is the declared nodata and in a
rolling release means *inside the footprint but not updated*. `-9999` is **not**
declared in the GeoTIFF metadata and means *outside the product footprint*
(~38 % of every CONUS product). Any gap test using a raster's own nodata value
silently accepts `-9999` as a fuel category. This caught out a first attempt at
measuring coverage during the session; do not repeat it.

**What was built.** `src/ingest/landfire_mosaic.py` stacks the staged releases
newest-first and fills each pixel from the newest release that has data there.
The CONUS mosaic:

    /data/jhaley/wrfxpy/landfire/LFmosaic_FBFM13_CONUS/Tif/LFmosaic_FBFM13_CONUS.tif
    2.61 GB, 156336 x 101538, EPSG:5070, 30 m, nodata 32767, overviews to /128
    LF2025 21.4 %, LF2024 40.4 %, older three 0.0 %, gaps 38.3 % (outside CONUS)
    + LFmosaic_FBFM13_CONUS.tif.manifest.json

The manifest records the releases used, their sizes and mtimes, and the pixels
each contributed. It was backfilled from the build log for this first raster
(the counts reconcile exactly with the grid size); future builds write it
automatically.

**Wiring.** `region_cfg` now holds only `alaska`, `hawaii`, `prvi`, and CONUS
fuels come from `geo_vars_path` in the base job configurations. See README §5a.

---

## 3. Files changed outside git — the part a rebuild would lose

None of this is in version control. Backups of everything touched are in the
session scratchpad, which **will not survive**, so re-create from here if needed.

**`etc/vtables/geo_vars.json` → `geo_vars.json_mosaic`** (was `geo_vars.json_2024`).
Revert with one `ln -sfn`. `geo_vars.json_2024` is untouched.

**`geo_vars.json_mosaic`** created in the main install; symlinked into the
testing install.

**Both `ngfs.json` files** — `region_cfg` reduced to `alaska`/`hawaii`/`prvi`,
with `prvi` newly added (it had been hardcoded in the package and *absent* from
`region_cfg`, so the monolith, whose guard was always live, had no PRVI handling
at all). The live config is now gitignored; `etc/ngfs.json.initial` is the
tracked example.

    /data/jhaley/wrfxpy/etc/ngfs.json                    (real file)
    /data/jhaley/new_wrfxpy/wrfxpy/etc/ngfs.json  ->  ngfs.json_testing

**`geo_vars_path` added to the base job configurations**, for provenance — each
job file now records the fuels it used. Four paths, **three** distinct files:

    /data/jhaley/wrfxpy/jobs/base_ngfs_cfg.json                    (real)
    /data/jhaley/wrfxpy/jobs/base_ngfs_cfg_short.json              (real)
    /data/jhaley/new_wrfxpy/.../jobs/base_ngfs_cfg.json        ->  base_ngfs_cfg.json_quick
    /data/jhaley/new_wrfxpy/.../jobs/base_ngfs_cfg_short.json  ->  main install's copy

**Consequence:** with `geo_vars_path` pinned in the base configs, the
`geo_vars.json` symlink no longer governs NGFS forecasts — it only affects
non-NGFS wrfxpy runs. To change NGFS fuels, edit `geo_vars.json_mosaic` or
rebuild the mosaic in place.

**Testing install's AK and HI vtables symlinked to the main install's.** They
had drifted badly: testing was on LF2022 for both while production used LF2025
(AK) and LF2023 (HI), so the two systems configured AK/HI fires differently —
a confound in exactly the comparison the testing install exists to make. Now
all four region vtables resolve identically from both installs.

**Deleted** `etc/ngfs.json_initial` (untracked, obsolete schema, backed up).

---

## 4. Corrections to earlier beliefs

- **The monolith always did read `region_cfg`.** Its guard is
  `'region_cfg' in ngfs_cfg.keys()` (`src/ngfs_start.py:426`), which is true.
  Only the *package* had the dead `region_cfg_REMOVE_THIS` key. An earlier
  claim in this session that `OR/WA/ID/MT/WY/CO` were getting the older vtable
  was therefore wrong for production — the monolith honoured the full 11-state
  list. Production was never mis-selecting CONUS fuels.
- **There is no `etc/*.json` gitignore wildcard.** `.gitignore` names
  `conf.json`, `tokens.json`, `clusters.json` individually, and now `ngfs.json`.
  The concrete `geo_vars.json_*` vtables are untracked, not ignored.
- **NAM and NAM218 are treated equivalently** by the system (legacy, from before
  the full menu of weather inputs existed), per the user. Dropping the
  `updates` region entry therefore cost the monolith nothing meteorologically.
- **PRVI was latently broken, not actually broken.** No PR/VI fire has ever
  been forecast; the four job files with `grib_source: GFSF` are hand-made
  Alaska tests from 2024.
- The repo convention for example configs is `<name>.initial` with a **dot**
  (`conf.json.initial`, `tokens.json.initial`), which is what the tracked files
  use.

---

## 5. What is verified, and how

- **Mosaic correctness.** All six New Mexico ignition points that LF2025 alone
  cannot serve now return real fuel codes; 400/400 ignition points have valid
  fuel, matching LF2024. In a 3000x3000 window in northern California, 437,407
  pixels genuinely differ between LF2024 and LF2025 and the mosaic follows
  LF2025 on **every one**, so it delivers newer fuels rather than a copy of
  2024. `-9999` is entirely absent from the output (0.00 %).
- **Mosaic mechanics.** `build()` was tested on small synthetic rasters with
  known fills, checking pixel offsets, fill precedence where both layers are
  fill, fall-through where only the newer is fill, and exact contribution
  counts across a multi-strip loop. That test found two real defects in the
  tool (fixed 256-px tiling and fixed overview levels both fail on rasters
  smaller than a block, which would have bitten HI and PRVI).
- **Region selection.** Eight state cases driven through the real
  `make_incident_configuration`: CA/NM/TX/OR record the mosaic; AK, HI, PR, VI
  each get their own fuels and grib source. Also checked that no CONUS state
  matches a region entry and no state matches two.
- **Live loop.** Cron runs at 08:32, 09:00, 09:32, 10:02 completed with zero
  tracebacks and no `KeyError`/`AttributeError`, 36–42 s, unchanged from
  baseline. Three ran during the mosaic build.
- **Not yet verified end to end:** no cron run since the change has actually
  processed an incident (they have all reported `Started 0 of 3x`), so the new
  region block has not been exercised *in the live loop* — only through the
  direct eight-case test. **First thing to check next session:** a job file
  newer than 2026-09-04 10:01 in either `jobs/` directory, and whether its
  `geo_vars_path` is `etc/vtables/geo_vars.json_mosaic` for a CONUS fire. A
  background watcher was left running for this but dies with the session.

---

## 6. Open items

**Ready to do:**

- **Push the six commits.** From a separate terminal.
- **AK mosaic is unnecessary.** `LF2025_FBFM13_AK` is 100 % valid with zero fill
  of either kind, confirmed at seven Alaska locations. The user judged this
  correctly; it needs nothing.
- **`LF2020_PRVI_IA` does not parse** in `landfire_mosaic --list`: its archive
  and directory names disagree and the TIFs sit a level deeper
  (`LF2020_PRVI_IA/LF2020_FBFM13_PRVI/Tif/`). PRVI works today via
  `geo_vars.json_prvi`, so this is cosmetic unless the tool should cover it.
- **The mosaic build reads all five layers** though three contribute nothing,
  costing ~7 of ~13 minutes of reads. Left alone deliberately: keeping every
  layer stays correct for regions where older releases do contribute. A
  `--layers N` flag would fix it if it ever matters.

**Carried over from 2026-09-03, still open:**

- `started_inc_ids` grows without bound (now ~3,790). The ids live inside the
  pickle so there is nowhere for a first-seen timestamp; bounding it properly
  means extracting an explicit ledger file. Still the user's design call.
- **Ledger semantics at cutover.** `start_forecast`'s two submission lines are
  still commented out, so the package records forecasts it never issued,
  including for VIIRS-only fires the monolith never saw. Decide what the
  accumulated ledger means before enabling submission, or those fires will be
  silently skipped.
- ~201 GB of monolith-format pickles (4,201 files) in the testing state
  directory, expected to be deleted eventually.
- Two pickle files matching neither class marker.
- `add_nifc_perims` has no callers and its configured `ngfs/perims/` does not
  exist in the testing install.

**Rough edges recorded in README §9** (found while documenting, none causing
failures): dead code in the package (`cluster_data` — which also has a
`db.labels_c_cloc` typo — `add_polar_data`, `prioritize_incidents`,
`find_old_detections`, and `set_cmd_str`, which is defined inside
`process_incident` and so is a never-called local); `unknown_incidents` resets
`skipped` inside the per-WFO loop so its count is wrong for more than one WFO
and would `NameError` on an empty list; `utilities.timestamp_from_string`
raises `NameError` on any call; `ngfs_day.py`'s `__main__` self-test calls a
keyword and a method that do not exist.

**One live bug not yet fixed:** the `behave`/`cawfe` command-line arguments are
matched correctly since `325ec0b`, but they assign
`ngfs_cfg["fire_namelist_path"]` at the top level while
`make_incident_configuration` reads `ngfs_cfg['run_cfg']['fire_namelist_path']`.
They are recognised and still have no effect. One-word fix, but it makes those
arguments live for the first time, so it is a behaviour change and was left for
the user.

---

## 7. Local working state — leave alone

Unchanged from the previous handoff and still true:

- `src/fmda/fuel_moisture_model.py` has a **staged** change that is the user's
  to commit. It survived every commit today; all commits were pathspec-limited
  so it was never swept in. Keep doing that.
- Nine other tracked files are modified by the user, including
  `src/ngfs_start.py`. The committed `ngfs_start.py` does not define
  `make_geo_folder`; it exists only in that working-copy modification.
- `core.fileMode` is `false` here; without it 181 files show as modified from a
  past `chmod`.
- ~12,000 untracked files and ~1 TB of data. **Never `git add -A`** — there is a
  45 GB pickle and multi-GB LANDFIRE archives, and `.gitignore` does not cover
  `landfire/`, `*.pkl`, `*.tif`, `*.geojson` or the root data dumps.

---

## 8. Rebuilding the mosaic next LANDFIRE cycle

    cd /data/jhaley/wrfxpy
    export PYTHONPATH=src
    export PROJ_LIB=/home/jhaley/anaconda3/envs/wrf_test/share/proj
    python -m ingest.landfire_mosaic --list
    python -m ingest.landfire_mosaic --product FBFM13 --region CONUS
    python -m ingest.landfire_mosaic --product FBFM13 --region CONUS --build

The third command is a dry run reporting each layer's estimated contribution;
`--build` writes the raster and manifest, about half an hour for CONUS. Stage
the new release as a directory named `LF<year>_FBFM13[_<version>]_CONUS` with
its GeoTIFF under `Tif/`, and nothing else needs to change — no code, no state
lists. The tool refuses to proceed if layers disagree on CRS, pixel size or
lattice alignment rather than resampling.
