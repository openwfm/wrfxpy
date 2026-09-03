# NGFS session handoff — 2026-09-03

Context for the next session. Read `/data/jhaley/wrfxpy/AGENTS.md` first: it
defines the approval gate, commit policy and development rules that govern work
here. Notably, outline changes and get approval **before** editing, and never
commit without explicit instruction.

---

## 1. System topology — the things that are not obvious

**There is only ONE copy of the `ngfs` package.**

    /data/jhaley/new_wrfxpy/wrfxpy/src/ngfs  ->  /data/jhaley/wrfxpy/src/ngfs   (symlink)

Editing `src/ngfs` in the main checkout changes the code the live loop runs, on
its next execution. Always work in a git worktree and merge deliberately; the
merge is the deploy.

**`/home/jhaley/work` is a symlink to `/data/jhaley`.** So `/home/jhaley/work/wrfxpy`
and `/data/jhaley/wrfxpy` are the same directory (verified by inode). Not a
second installation.

**Two systems run in parallel:**

| | monolith (production) | new package (testing) |
|---|---|---|
| entry point | `/data/jhaley/wrfxpy/src/ngfs_start.py` | `/data/jhaley/wrfxpy/src/ngfs/ngfs_start_2.py` |
| runs from | `/data/jhaley/wrfxpy` | `/data/jhaley/new_wrfxpy/wrfxpy` |
| state dir | `/data/jhaley/wrfxpy/ngfs` (266 GB) | `/data/jhaley/new_wrfxpy/wrfxpy/ngfs` (220 GB) |
| issues forecasts? | yes | no — submission deliberately disabled |
| driven by | GOES | GOES + NGFS-processed VIIRS (sees fires the monolith cannot) |

The new system runs every 30 min via `/data/jhaley/new_wrfxpy/wrfxpy/cron_ngfs.sh`,
logging to `log_cron_ngfs.log` (overwritten each run, then copied to
`ngfs/cron_ngfs_<timestamp>.log`).

`src/ngfs_start.py` is the non-authoritative precursor but is the system
actually producing forecasts. ~35 sibling scripts in `src/` import its classes,
so it cannot simply be deleted.

---

## 2. Environment and invocation

- **Runtime env is `wrf_test`**, not `wrfx`. It is the only conda env with
  sklearn + Basemap + folium together (`wrfx` has no sklearn). Python 3.7.7,
  pandas 1.3.5.
- **`PROJ_LIB` must be set** when invoking the interpreter directly instead of
  `conda activate wrf_test` — Basemap reads it at import and raises
  `KeyError: 'PROJ_LIB'`:
  `PROJ_LIB=/home/jhaley/anaconda3/envs/wrf_test/share/proj`
- **No pyarrow, fastparquet or pytables** in any env, so Parquet and HDF5 are
  unavailable without installing into a scientific env (which AGENTS.md forbids
  without approval). gzip/bz2/lzma are stdlib and were used instead.
- The `PYTHONPATH=src` line in `cron_ngfs.sh` does nothing — no `export`, so it
  never reaches python. It works only because `ngfs_start_2.py` does its own
  `sys.path.insert(1, 'src/')`.
- **Pushing to origin does not work from inside Claude Code.** No ssh-agent
  runs, `~/.ssh/id_rsa_github` is passphrase-encrypted, no askpass helper
  exists. The user pushes from a separate terminal. He explicitly declined
  changing the ssh setup — do not propose it again.
- A worktree defaults to branching from `origin/master`, which has no
  `src/ngfs`. Create them explicitly from `james_ngfs`:
  `git worktree add -b <name> .claude/worktrees/<dir> james_ngfs`
- A worktree lacks the ~74 untracked `src/*.py` modules and gitignored
  `etc/*.json`. To make one importable for testing, symlink them in from the
  main checkout, then remove them before committing.
- The user dislikes branch proliferation (the repo has dozens). Land work on
  `james_ngfs` by fast-forward and delete the temporary branch.

---

## 3. What changed on 2026-09-03

Ten commits, `2e7aa35` → `9e35dda`, all on `james_ngfs` and **pushed**.

    dbcb2c5  Add AGENTS.md defining project instructions and agent operating rules
    ad3207a  Add src/ngfs package: modules required to run ngfs_start_2.py
    b64fc5b  Add NGFS visualization modules: map_locations, incident_webpage, watch_duty
    4134939  Add perim_cache.py: NIFC fire perimeter caching
    3afa99c  Compress state pickles with gzip and write them atomically
    193f099  Stop config_manager importing make_geo_folder from the monolith
    bab9052  Remove duplicate get_fmda_path from config_manager
    daae3e7  Add helper to compress previously saved state pickles
    6ccd78d  Make get_old_incidents survive missing or unreadable state
    9e35dda  Read the perimeter config key the deployed config defines

**Version control.** `src/ngfs` was entirely untracked; 19 of 25 files now
tracked. Deliberately still untracked: `latency.py`, `upload_forecasts.py`,
`watchduty_nifc_table.py`, `utilities.py` (tooling), `tmp.txt` (stray diff
dump), `ngfs_api.json` (downloaded OpenAPI spec nothing reads).

**Pickle compression.** `save_pickle` writes gzip level 6 via `pd.to_pickle` to
a pid-qualified temp name, then `os.replace`. Constants in `constants.py`:
`PICKLE_COMPRESSION`, `PICKLE_SUFFIX`, `PICKLE_PATTERNS`.
`persistence.state_pickle_files` accepts `.pkl`, `.pkl.gz`, `.pkl.xz`, so old
files stay readable and no migration was needed. Measured on 70.9 MB: gzip-6 =
7.1x at 2.11 s write / 0.59 s read; level 9 = 7.23 s for 3 % more; xz = 12.1x
at 14.6 s. Older pickles reach 8.6–8.8x.

**Backlog sweep.** `persistence.compress_state_pickles` gzips at the byte level,
verifies by SHA-256 of decompressed bytes, deletes the original only on a match.
Result: 1,996 files, 202.50 GB → 25.50 GB, zero failures; directory 397 → 220 GB.

**Import hygiene.** `config_manager` no longer imports from `ngfs_start`, and
`get_fmda_path` is no longer duplicated. `ngfs_start` is now absent from the
package import graph — before this, a clean checkout of `james_ngfs` could not
import the package at all.

---

## 4. Verified working in production

- Compressed pickles are **written** by `save_pickle` (13 overnight, 9.2–10 MB
  each) and **read back** to restore state on cron runs.
- Cron runs complete with zero tracebacks.
- `origin/james_ngfs` == local == `9e35dda`.

---

## 5. Facts worth not rediscovering

- **Pickle class paths differ by writer.** The monolith runs as a script, so its
  pickles record `__main__.ngfs_day`; package-written ones record
  `ngfs.ngfs_day.ngfs_day`. Reading a `__main__` pickle from the package raises
  `AttributeError: module '__main__' has no attribute 'ngfs_day'`. Determine
  ownership by reading the first 128 bytes — `persistence.pickle_writer` does
  this. Do **not** rely on the `_testing` filename suffix; `set_save_name` is
  annotated for removal.
- **`ngfs_start.py` globs only `*.pkl`.** Compressing anything in the monolith's
  own state directory would hide it and cause re-forecasting.
  `compress_state_pickles` refuses such directories.
- **mtime is load-bearing.** Both `get_old_incidents` implementations sort by
  mtime and apply a 7-day filter. Any tool rewriting a state file must preserve
  timestamps or it promotes stale state to "newest" and corrupts the ledger.
- **`save_pickle` is gated on `start_count > 0`** (`ngfs_day.save_outputs`), so
  state is written only by runs that started a forecast. Pickle timestamps are
  therefore *not* a record of when the loop ran.
- The 4,201 legacy monolith pickles in the new install's state dir are excluded
  from `get_old_incidents` only by the 7-day age filter. Since `6ccd78d` an
  unreadable pickle is logged and skipped rather than fatal.
- `ngfs_directory` in config is the relative string `"ngfs"`, resolved against
  the process working directory.

---

## 6. Open items

**User is still thinking about:**

- `started_inc_ids` grows without bound — **3,790 entries**, and every new
  incident id is compared against the whole list. His idea: store a first-seen
  timestamp per id and compare only against roughly the last two months. Note
  the entries live inside the pickle, so there is nowhere to put a timestamp;
  that argues for extracting a small explicit ledger file (id + first-seen date)
  as authoritative, with the pickle demoted to a data cache.

**Known, deferred:**

- **201 GB of monolith-format pickles** (4,201 files) in the new install's state
  directory. Nothing appears to read them; he expects to delete them later.
  `compress_state_pickles` skips them by design.
- **2 pickle files** whose header matched neither class marker — inspect before
  assuming they are junk.
- `add_nifc_perims` has **no callers** and its configured `ngfs/perims/` does
  not exist in the new install. `9e35dda` fixed the config key it reads; wiring
  it up is a separate decision.
- **Ledger semantics at cutover.** `start_forecast` deliberately has its two
  submission lines commented out so the system behaves as if it forecast — job
  files are written for comparison against the monolith. Consequence:
  `started_inc_ids` records forecasts never issued, including for VIIRS-only
  fires the monolith never saw. When submission is enabled, decide what that
  accumulated ledger should mean, or those fires will be silently skipped.
- Structural: each pickle stores a 48-hour rolling window written every 30 min,
  so each detection is stored ~96 times across the series. The redundancy is
  *between* files, not inside them. That is the real case for storing references
  to cached ingest files instead of the data — but re-deriving frames is not
  byte-reproducible (order-dependent `drop_duplicates`, wall-clock windows,
  hourly FIRMS URT snapshots), so any such scheme needs a byte-exact sidecar
  rather than replaying the pipeline.

---

## 7. Local working state — leave alone

- `src/fmda/fuel_moisture_model.py` has a **staged** change, blob
  `01673b43289815e738806e70a9f69e6f3cb3d283` (+21/−33: generalises FMDA from 3
  fuel-moisture classes to arbitrary `n`, plus a `cPickle`→`pickle` fix and a
  Kalman-gain index correction). It survived every merge this session and is the
  user's to commit. Use pathspec-limited commits (`git commit -- <paths>`) so it
  is never swept in.
- 8 other tracked files are modified, including `src/ngfs_start.py` (+249/−89).
  The **committed** `ngfs_start.py` does not define `make_geo_folder` — it
  exists only in that working-copy modification.
- `core.fileMode` was set to `false` in this repo; without it 181 files show as
  modified purely from a past `chmod`. Revert with
  `git config --unset core.fileMode`.
- ~12,000 untracked files and ~1 TB of data live here. **Never `git add -A`** —
  there is a 45 GB pickle and multi-GB LANDFIRE archives. `.gitignore` does not
  cover `landfire/`, `*.pkl`, `*.tif`, `*.geojson` or the root data dumps.

---

## 8. Running the compression helper

    cd /data/jhaley/new_wrfxpy/wrfxpy
    PROJ_LIB=/home/jhaley/anaconda3/envs/wrf_test/share/proj \
    PYTHONPATH=/data/jhaley/wrfxpy/src \
    /home/jhaley/anaconda3/envs/wrf_test/bin/python -m ngfs.persistence ngfs [--delete] [--limit N]

Without `--delete` it reports only. `--limit 0` with `--delete` removes verified
originals for already-compressed files without compressing anything new. Files
newer than `--min-age-hours` (default 6) are never touched, so a run in progress
is safe, and interrupted sweeps can be re-run.
