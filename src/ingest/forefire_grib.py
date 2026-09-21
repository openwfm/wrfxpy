#!/usr/bin/env python
"""Build ForeFire input netcdfs from GRIB winds, with no WRF-SFIRE in the loop.

`forefire.make_FF_nc` takes its four ingredients from a completed WRF-SFIRE run:
`NFUEL_CAT` and `ZSF` from `wrfinput_d01`, `UF`/`VF` from the wrfout.  That run is
the latency blocker -- a ForeFire forecast cannot start until WRF finishes.  Every
ingredient has a non-WRF source already staged on this machine:

    fuel      LFmosaic_FBFM13_CONUS.tif   (the raster geogrid itself uses)
    altitude  LF2020_Elev_220_CONUS       (matches ZSF to RMS 2.2 m)
    windU/V   WindNinja downscaling a raw HRRR wrfprs GRIB
    domain    synthesised for the grid built around the ignition

so this module builds the same netcdfs from an ignition time and location alone.

It writes files named exactly as `forefire.ff_nc_name` expects, into the ForeFire
run directory.  `forefire.make_script_set` skips its own builder for any step whose
netcdf already exists, so the existing sequence machinery -- `make_script_set`,
`run_timing_table`, the `.ff` templates -- runs unmodified on these.  That
machinery already produces the time-lagged ensemble: each step advances one hour,
saves state with `print[sim_..._END_STEP.ff]`, then runs on to the common
`END_TIME` and writes `final_<step>.geojson`.

Environments: this runs under `wrf_test` (Python 3.7) like the rest of the ForeFire
work, and shells out to the separate `windninja` env for WindNinja_cli and GDAL.
That env must be *activated*, not merely put on PATH, or PROJ_DATA is unset and
gdalwarp fails inside PROJ with a misleading proj.db error.
"""

from __future__ import print_function

import os
import os.path as osp
import argparse
import glob
import logging
import re
import subprocess

import numpy as np
import pandas as pd


# ---------------------------------------------------------------- configuration

WINDNINJA_ENV = 'windninja'
CONDA_SH = '/home/jhaley/anaconda3/etc/profile.d/conda.sh'

#NFUEL_CAT points at the same mosaic as etc/vtables/geo_vars.json, so a GRIB-driven
#run and a geogrid run see identical fuels.
FUEL_TIF = '/data/jhaley/wrfxpy/landfire/LFmosaic_FBFM13_CONUS/Tif/LFmosaic_FBFM13_CONUS.tif'
ELEV_TIF = '/data/jhaley/wrfxpy/landfire/LF2020_Elev_220_CONUS/Tif/LC20_Elev_220.tif'

#hourly f03 HRRR cache filled by the FMDA HRRR cycler.  f03 only, so this supports
#an ignition in the recent past but cannot reach ahead of now; when forecast cycles
#are cached, point --grib-dir at those and widen `_grib_valid_time`.
HRRR_CACHE = '/data/jhaley/clean_wrfxpy/wrfxpy/ingest/HRRRA'

#per-fuel-category wind reduction, from etc/nlists/default.fire.  UF/VF in a wrfout
#already carry this, so an external field written into the slot they fill must have
#it applied or the fire is driven ~3.3x too fast.
WINDRF_13 = np.array([0.36, 0.36, 0.44, 0.55, 0.42, 0.44, 0.44,
                      0.36, 0.36, 0.36, 0.36, 0.43, 0.46, 1e-7])
NO_FUEL_CAT = 14
#category 14 cannot burn and its 1e-7 is a sentinel, not a reduction.  Applying it
#would punch near-zero wind holes into the field, so non-burnable cells take the
#modal burnable factor instead, keeping the field continuous across them.
NO_FUEL_WINDRF = 0.36

#WindNinja's measured ceiling lies between 4.0 M cells (runs) and 5.76 M (fails);
#the failure is an int index overflow reported as an out-of-memory error, so more
#RAM does not move it.
MAX_CELLS = 4000000


class ForeFireGribError(Exception):
    pass


# --------------------------------------------------------------- shelling out

def _bash(script, what):
    """Run a bash snippet with the windninja env activated; return stdout."""
    cmd = '. {} && conda activate {} && {}'.format(CONDA_SH, WINDNINJA_ENV, script)
    proc = subprocess.Popen(['bash', '-lc', cmd],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out, _ = proc.communicate()
    out = out.decode('utf-8', 'replace')
    if proc.returncode != 0:
        raise ForeFireGribError('{} failed (exit {}):\n{}'.format(what, proc.returncode, out))
    return out


def _transform(x, y, s_epsg, t_epsg):
    out = _bash('echo "{} {}" | gdaltransform -s_srs EPSG:{} -t_srs EPSG:{} -output_xy'.format(
        x, y, s_epsg, t_epsg), 'gdaltransform')
    parts = out.split()
    if len(parts) < 2:
        raise ForeFireGribError('gdaltransform returned no coordinate:\n{}'.format(out))
    return float(parts[0]), float(parts[1])


# -------------------------------------------------------------------- domain

def utm_epsg(lat, lon):
    """EPSG code of the UTM zone containing a point."""
    zone = int((lon + 180.0) / 6.0) + 1
    return (32600 if lat >= 0 else 32700) + zone


def build_domain(lat, lon, domain_km=30.0, fire_dx=30.0):
    """Define the square fire grid centred on an ignition point.

    The grid is in the local UTM zone at fire_dx metres.  The 30 m default is
    LANDFIRE's native resolution, so fuels and terrain are cut rather than
    resampled.
    """
    nx = int(round(domain_km * 1000.0 / fire_dx))
    ny = nx
    if nx * ny > MAX_CELLS:
        raise ForeFireGribError(
            '{} x {} = {:.2f} M cells is past the measured WindNinja ceiling of '
            '{:.1f} M cells; reduce --domain-km or coarsen --fire-dx'.format(
                nx, ny, nx * ny / 1e6, MAX_CELLS / 1e6))

    epsg = utm_epsg(lat, lon)
    cx, cy = _transform(lon, lat, 4326, epsg)
    half_x = nx * fire_dx / 2.0
    half_y = ny * fire_dx / 2.0
    x0, y0 = cx - half_x, cy - half_y
    x1, y1 = x0 + nx * fire_dx, y0 + ny * fire_dx

    #corner lon/lat for the netcdf BBoxWSEN.  Taking SW and NE corners is the same
    #simplification make_FF_nc makes for the LCC case.
    west, south = _transform(x0, y0, epsg, 4326)
    east, north = _transform(x1, y1, epsg, 4326)

    return {
        'lat': lat, 'lon': lon, 'epsg': epsg,
        'nx': nx, 'ny': ny, 'dx': float(fire_dx),
        'x0': x0, 'y0': y0, 'x1': x1, 'y1': y1,
        'west': west, 'south': south, 'east': east, 'north': north,
        'Lx': float(nx * fire_dx), 'Ly': float(ny * fire_dx),
    }


# -------------------------------------------------------------- raster cuts

def _warp(src, dst, domain, resample, fmt='AAIGrid', extra=''):
    """Cut and reproject a raster onto the exact fire grid."""
    _bash('gdalwarp -q -overwrite -t_srs EPSG:{epsg} -te {x0} {y0} {x1} {y1} '
          '-ts {nx} {ny} -r {r} -of {fmt} {extra} "{src}" "{dst}"'.format(
              epsg=domain['epsg'], x0=domain['x0'], y0=domain['y0'],
              x1=domain['x1'], y1=domain['y1'], nx=domain['nx'], ny=domain['ny'],
              r=resample, fmt=fmt, extra=extra, src=src, dst=dst),
          'gdalwarp of {}'.format(osp.basename(src)))
    return dst


def read_asc(path):
    """Read an AAIGrid into a south-to-north array.

    AAIGrid writes row 0 at the *north* edge; WRF fields and ForeFire's domain
    (SWx/SWy at the southwest corner) are indexed south-to-north, so the array is
    flipped.  Getting this wrong mirrors the fire north-south and is invisible on
    a symmetric domain.
    """
    a = np.loadtxt(path, skiprows=6)
    return np.flipud(a)


def remap_fuel(raw):
    """LANDFIRE FBFM13 codes -> NFUEL_CAT, following etc/../var_wisdom.py.

    Anderson categories 1-13 pass through.  15-90, 92 and 94-99 are non-burnable
    and become 14.  0, 91 (urban) and 93 (agriculture) take their nearest valid
    neighbour, as geogrid's 'nearest' fill does, so a fire is not stopped by a
    road or a field boundary that LANDFIRE happens to label unburnable.
    """
    f = np.rint(raw).astype(np.int32)
    out = np.full(f.shape, NO_FUEL_CAT, dtype=np.int16)
    valid = (f >= 1) & (f <= 13)
    out[valid] = f[valid].astype(np.int16)

    nearest = (f == 0) | (f == 91) | (f == 93)
    if nearest.any() and valid.any():
        from scipy import ndimage
        #indices of the nearest zero element of ~valid, i.e. the nearest valid cell
        idx = ndimage.distance_transform_edt(~valid, return_distances=False,
                                             return_indices=True)
        out[nearest] = out[idx[0][nearest], idx[1][nearest]]
    return out


# ---------------------------------------------------------------- WindNinja

def _wn_output(out_dir, mesh_m):
    """Find the vel/ang pair WindNinja just wrote."""
    vel = sorted(glob.glob(osp.join(out_dir, '*_vel.asc')))
    ang = sorted(glob.glob(osp.join(out_dir, '*_ang.asc')))
    if not vel or not ang:
        raise ForeFireGribError(
            'WindNinja wrote no vel/ang pair into {} -- it exits 0 on inputs it '
            'cannot read, so check the run log above'.format(out_dir))
    return vel[-1], ang[-1]


def windninja_uv(grib, dem_tif, domain, work_dir, time_zone, mesh_m=None,
                 wind_height=6.096, num_threads=32, vegetation='grass'):
    """Downscale one GRIB and return (u, v) on the fire grid, in m/s.

    Returns the *unreduced* wind at wind_height; windrf is applied separately so
    the reduction stays visible and per-fire tunable.
    """
    mesh_m = float(mesh_m or domain['dx'])
    out_dir = osp.join(work_dir, 'wn_' + osp.basename(grib).replace('.', '_'))
    if not osp.isdir(out_dir):
        os.makedirs(out_dir)

    #--output_speed_units defaults to mph, which would inflate the field 2.24x with
    #nothing looking wrong.  --time_zone is mandatory whenever --forecast_filename
    #is given, and the timestamp WindNinja puts in the output filename is in that
    #zone rather than UTC.
    log = _bash(
        'WindNinja_cli --num_threads {nt} --elevation_file "{dem}" '
        '--initialization_method wxModelInitialization --forecast_filename "{grib}" '
        '--time_zone {tz} --mesh_resolution {mesh} --units_mesh_resolution m '
        '--output_speed_units mps --output_wind_height {h} '
        '--units_output_wind_height m --vegetation {veg} '
        '--write_ascii_output true --output_path "{out}"'.format(
            nt=num_threads, dem=dem_tif, grib=grib, tz=time_zone, mesh=mesh_m,
            h=wind_height, veg=vegetation, out=out_dir),
        'WindNinja on {}'.format(osp.basename(grib)))
    logging.debug(log)

    vel_asc, ang_asc = _wn_output(out_dir, mesh_m)
    #WindNinja solves on its own mesh and its extent can be inset from the DEM, so
    #both fields are warped onto the exact fire grid before use.
    vel = read_asc(_warp(vel_asc, osp.join(out_dir, 'vel_grid.asc'), domain, 'bilinear'))
    ang = read_asc(_warp(ang_asc, osp.join(out_dir, 'ang_grid.asc'), domain, 'near'))

    #_ang is the direction the wind blows *from*, in degrees; convert to vector
    #components blowing towards.
    rad = np.deg2rad(ang)
    u = -vel * np.sin(rad)
    v = -vel * np.cos(rad)

    _check_wind(vel, grib)
    return u.astype(np.float32), v.astype(np.float32)


def _check_wind(vel, grib):
    """Reject the silent-garbage failure mode.

    WindNinja exits 0 on a raw NAM218 file and produces speeds of 25.8 to 45,632
    m/s, so the exit code proves nothing and the values must be checked.
    """
    finite = vel[np.isfinite(vel)]
    if finite.size == 0:
        raise ForeFireGribError('WindNinja produced no finite speeds for {}'.format(grib))
    hi = float(np.nanmax(finite))
    if hi > 100.0:
        raise ForeFireGribError(
            'WindNinja produced a peak speed of {:.1f} m/s for {} -- this is the '
            'signature of an input it cannot actually read (it exits 0 on those). '
            'Confirm the file is an HRRR wrfprs GRIB.'.format(hi, grib))


def apply_windrf(u, v, fuel, windrf=None):
    """Apply the per-fuel-category wind reduction UF/VF already carry."""
    table = np.asarray(WINDRF_13 if windrf is None else windrf, dtype=np.float64)
    if table.size != 14:
        raise ForeFireGribError('windrf table must have 14 entries, got {}'.format(table.size))
    factor = np.full(fuel.shape, NO_FUEL_WINDRF, dtype=np.float64)
    burnable = (fuel >= 1) & (fuel <= 13)
    factor[burnable] = table[fuel[burnable] - 1]
    return (u * factor).astype(np.float32), (v * factor).astype(np.float32)


# --------------------------------------------------------------- the netcdf

def write_ff_nc(out_path, domain, fuel, altitude, u, v):
    """Write one ForeFire input netcdf, matching make_FF_nc's structure exactly."""
    import netCDF4 as nc

    ny, nx = fuel.shape
    for name, arr in (('altitude', altitude), ('windU', u), ('windV', v)):
        if arr.shape != (ny, nx):
            raise ForeFireGribError('{} is {} but fuel is {}'.format(name, arr.shape, (ny, nx)))

    bbox_wsen = '{},{},{},{}'.format(domain['west'], domain['south'],
                                     domain['east'], domain['north'])
    wsenlbrt = '{},{},{},{},0.0,0.0,{},{}'.format(
        domain['west'], domain['south'], domain['east'], domain['north'],
        domain['Lx'], domain['Ly'])

    with nc.Dataset(out_path, 'w', format='NETCDF4') as ff:
        ff.createDimension('nx', nx)
        ff.createDimension('ny', ny)
        ff.createDimension('nz', 1)
        ff.createDimension('nt', 1)
        ff.createDimension('fx', nx)
        ff.createDimension('fy', ny)
        ff.createDimension('fz', 1)
        ff.createDimension('ft', 1)
        ff.createDimension('wind_rows', ny)
        ff.createDimension('wind_columns', nx)
        ff.createDimension('wind_dimensions', 1)
        ff.createDimension('wind_directions', 1)

        dom = ff.createVariable('domain', str)
        dom.type = 'domain'
        dom.BBoxWSEN = bbox_wsen
        dom.WSENLBRT = wsenlbrt
        dom.SWx = np.float32(0.0)
        dom.SWy = np.float32(0.0)
        dom.Lx = np.float32(domain['Lx'])
        dom.Ly = np.float32(domain['Ly'])
        dom.Lz = np.float32(0.0)
        dom.t0 = np.float32(0.0)
        dom.Lt = np.float32(np.inf)
        dom.SWz = np.float32(0.0)

        fv = ff.createVariable('fuel', 'i2', ('ft', 'fz', 'fy', 'fx'), fill_value=False)
        fv.type = 'fuel'
        fv[:] = np.ascontiguousarray(fuel[np.newaxis, np.newaxis, :, :].astype(np.int16))

        av = ff.createVariable('altitude', 'i2', ('nt', 'nz', 'ny', 'nx'), fill_value=False)
        av.type = 'data'
        av[:] = np.ascontiguousarray(altitude[np.newaxis, np.newaxis, :, :].astype(np.int16))

        uv = ff.createVariable('windU', 'f4',
                               ('wind_dimensions', 'wind_directions', 'wind_rows', 'wind_columns'),
                               fill_value=np.nan)
        uv.type = 'wind'
        uv[:] = np.ascontiguousarray(u[np.newaxis, np.newaxis, :, :])

        vv = ff.createVariable('windV', 'f4',
                               ('wind_dimensions', 'wind_directions', 'wind_rows', 'wind_columns'),
                               fill_value=np.nan)
        vv.type = 'wind'
        vv[:] = np.ascontiguousarray(v[np.newaxis, np.newaxis, :, :])

        #provenance, so a netcdf built this way is never mistaken for a WRF one
        ff.source = 'forefire_grib.py: LANDFIRE fuels/terrain + WindNinja winds'
        ff.epsg = domain['epsg']
        ff.fire_dx = domain['dx']


# ------------------------------------------------------------- timing table

_GRIB_RE = re.compile(r'hrrr\.t(\d{2})z\.wrfprsf(\d{2})\.grib2$')


def as_naive_utc(when):
    """Timestamp as naive UTC.

    GRIB valid times are built naive from filenames, while an ignition string
    carrying a 'Z' parses tz-aware, and comparing the two raises.  The rest of the
    ForeFire code is naive throughout, so everything is normalised to that here
    rather than each comparison growing its own guard.
    """
    ts = pd.Timestamp(when)
    if ts.tzinfo is not None:
        ts = ts.tz_convert('UTC').tz_localize(None)
    return ts


def _grib_valid_time(path):
    """Valid time of a cached HRRR wrfprs file, from its cycle and lead hour."""
    base = osp.basename(path)
    m = _GRIB_RE.search(base)
    if not m:
        return None
    day = re.search(r'hrrr\.(\d{8})', path)
    if not day:
        return None
    cycle = pd.Timestamp('{}T{}:00:00'.format(day.group(1), m.group(1)))
    return cycle + pd.Timedelta(hours=int(m.group(2)))


def find_gribs(start_utc, end_utc, grib_dir=None):
    """Cached HRRR files whose valid time falls in [start, end], in time order."""
    grib_dir = grib_dir or HRRR_CACHE
    start_utc, end_utc = as_naive_utc(start_utc), as_naive_utc(end_utc)
    found = {}
    for path in glob.glob(osp.join(grib_dir, 'hrrr.*', '*', '*.grib2')):
        vt = _grib_valid_time(path)
        if vt is not None and start_utc <= vt <= end_utc:
            #one file per valid time; a later lead wins nothing, so keep the first
            found.setdefault(vt, path)
    return [(vt, found[vt]) for vt in sorted(found)]


def grib_timing_table(ign_utc, gribs):
    """Timing table with make_timing_table's columns, from GRIB valid times.

    `ign_seconds` follows the same convention: -9999 before the ignition, then
    seconds measured from the step that contains it, which is ForeFire's t=0.
    """
    if len(gribs) < 2:
        raise ForeFireGribError(
            'need at least 2 GRIBs to define a time step, found {}'.format(len(gribs)))
    ig_time = as_naive_utc(ign_utc)
    times = [vt for vt, _ in gribs]
    t_step = int((times[1] - times[0]).total_seconds())
    ignition_seconds = int((ig_time - times[0]).total_seconds())

    rows, ws_ignition = [], 0.0
    for (vt, path) in gribs:
        ws = int((vt - times[0]).total_seconds())
        if ws + t_step < ignition_seconds:
            ig_s = -9999
        elif ws < ignition_seconds and ws + t_step > ignition_seconds:
            ws_ignition = ws
            ig_s = ignition_seconds - ws
        else:
            ig_s = ws - ws_ignition
        rows.append({
            'grib': path,
            'UTC_str': vt.strftime('%Y-%m-%dT%H:%M:%SZ'),
            'wrf_seconds': ws,
            'ign_seconds': int(ig_s),
        })
    return pd.DataFrame(rows)


# ------------------------------------------------------------- orchestration

def build_step_ncs(lat, lon, ign_utc, out_dir, grid_code, steps=8, domain_km=30.0,
                   fire_dx=30.0, grib_dir=None, time_zone='America/Denver',
                   wn_mesh=None, windrf=None, work_dir=None, overwrite=False):
    """Build the per-step ForeFire netcdfs for one fire.

    Returns the timing table.  Writes `FF_<step>_<grid_code>_<UTC>.nc` into
    out_dir, which is what make_script_set looks for and skips rebuilding.
    """
    ig_time = as_naive_utc(ign_utc)
    #one extra valid time: the sequence needs a step beyond the last one it runs,
    #because each step's END_STEP is the next row's ign_seconds.
    gribs = find_gribs(ig_time, ig_time + pd.Timedelta(hours=steps + 1), grib_dir)
    if len(gribs) < steps + 1:
        logging.warning('asked for %d steps but the cache covers %d valid times '
                        'from %s -- the run will be shorter',
                        steps, len(gribs), ig_time)
    tt = grib_timing_table(ign_utc, gribs)

    if not osp.isdir(out_dir):
        os.makedirs(out_dir)
    work_dir = work_dir or osp.join(out_dir, '_grib_work')
    if not osp.isdir(work_dir):
        os.makedirs(work_dir)

    domain = build_domain(lat, lon, domain_km=domain_km, fire_dx=fire_dx)
    logging.info('domain %d x %d at %.0f m, UTM EPSG:%d',
                 domain['nx'], domain['ny'], domain['dx'], domain['epsg'])

    #static fields are cut once and reused by every step
    dem_tif = _warp(ELEV_TIF, osp.join(work_dir, 'dem.tif'), domain, 'bilinear', fmt='GTiff')
    dem_asc = _warp(dem_tif, osp.join(work_dir, 'dem.asc'), domain, 'near')
    altitude = read_asc(dem_asc)
    fuel_asc = _warp(FUEL_TIF, osp.join(work_dir, 'fuel.asc'), domain, 'near')
    fuel = remap_fuel(read_asc(fuel_asc))
    logging.info('fuels: %d of %d cells burnable; terrain %.0f-%.0f m',
                 int(((fuel >= 1) & (fuel <= 13)).sum()), fuel.size,
                 float(np.nanmin(altitude)), float(np.nanmax(altitude)))

    built = []
    for i, row in tt.iterrows():
        #mirror make_script_set: the last row only supplies END_STEP, and spinup
        #rows before the ignition are not run
        if i >= len(tt) - 1 or row['ign_seconds'] <= 0:
            continue
        name = 'FF_{}_{}_{}.nc'.format(str(i).zfill(2), grid_code, row['UTC_str'])
        target = osp.join(out_dir, name)
        if osp.exists(target) and not overwrite:
            logging.info('%s exists, skipping', name)
            built.append(target)
            continue
        logging.info('step %d: %s <- %s', i, row['UTC_str'], osp.basename(row['grib']))
        u, v = windninja_uv(row['grib'], dem_tif, domain, work_dir, time_zone,
                            mesh_m=wn_mesh)
        u, v = apply_windrf(u, v, fuel, windrf=windrf)
        write_ff_nc(target, domain, fuel, altitude, u, v)
        built.append(target)

    logging.info('built %d netcdfs in %s', len(built), out_dir)
    tt.to_csv(osp.join(out_dir, 'timing_table_grib.csv'), index=False)
    return tt


def main():
    p = argparse.ArgumentParser(
        description='Build ForeFire input netcdfs from an ignition and HRRR GRIBs.')
    p.add_argument('lat', type=float, help='ignition latitude')
    p.add_argument('lon', type=float, help='ignition longitude')
    p.add_argument('ign_utc', help='ignition time, e.g. 2026-09-15T17:00:00Z')
    p.add_argument('out_dir', help='ForeFire run directory to write the netcdfs into')
    p.add_argument('--grid-code', default=None, help='defaults to a name built from the ignition')
    p.add_argument('--steps', type=int, default=8,
                   help='forecasts in the time-lagged ensemble (default 8)')
    p.add_argument('--domain-km', type=float, default=30.0)
    p.add_argument('--fire-dx', type=float, default=30.0,
                   help='fire grid spacing in m (default 30, LANDFIRE native)')
    p.add_argument('--wn-mesh', type=float, default=None,
                   help='WindNinja solve mesh in m; defaults to --fire-dx. Coarsening '
                        'here is interpolation, not resolved detail')
    p.add_argument('--grib-dir', default=None, help='default: the FMDA HRRR cache')
    p.add_argument('--time-zone', default='America/Denver',
                   help='required by WindNinja; only affects its output filenames')
    p.add_argument('--windrf', default=None,
                   help='14 comma-separated per-category factors, overriding the '
                        'namelist defaults; set this per fire')
    p.add_argument('--overwrite', action='store_true')
    p.add_argument('--verbose', action='store_true')
    args = p.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format='%(levelname)s %(message)s')

    windrf = None
    if args.windrf:
        windrf = [float(x) for x in args.windrf.split(',')]
    grid_code = args.grid_code or 'FIRE_{}_{}'.format(
        args.ign_utc.replace(':', '').replace('-', ''), abs(int(args.lat * 100)))

    build_step_ncs(args.lat, args.lon, args.ign_utc, args.out_dir, grid_code,
                   steps=args.steps, domain_km=args.domain_km, fire_dx=args.fire_dx,
                   grib_dir=args.grib_dir, time_zone=args.time_zone,
                   wn_mesh=args.wn_mesh, windrf=windrf, overwrite=args.overwrite)


if __name__ == '__main__':
    main()
