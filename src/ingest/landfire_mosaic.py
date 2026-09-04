"""
Build a LANDFIRE fuels mosaic from the staged product directories.

LANDFIRE refreshes its fuels products in a rolling fashion: a given release
covers only the areas updated that cycle, and the rest of the raster is fill.
LF2025_FBFM13_CONUS, for instance, carries valid fuel codes over roughly the
western third of CONUS and nothing elsewhere, despite its name and a
CONUS-wide bounding box. No single release is therefore usable on its own as a
CONUS fuels source, and choosing between releases by state does not work
either, because the update regions cut across state lines.

This module discovers the staged releases, stacks them newest-first, and fills
each pixel from the newest release that has data there. The result is one
raster carrying the newest available fuel category everywhere any release has
one.

Two fill values appear in these products and both are treated as gaps:

    32767   declared nodata; in a rolling release this marks
            'inside the footprint but not updated this cycle'
    -9999   NOT declared in the GeoTIFF metadata; marks
            'outside the product footprint' (ocean, Canada, Mexico)

Because -9999 is undeclared, any gap test based on the raster's own nodata
value alone silently accepts it as a fuel category. The output collapses both
to a single declared nodata so downstream consumers cannot make that mistake.

Overviews are built with nearest-neighbour resampling. Fuel model codes are
categorical, so averaging them would invent categories that do not exist.

Usage:

    python -m ingest.landfire_mosaic --list
    python -m ingest.landfire_mosaic --product FBFM13 --region CONUS
    python -m ingest.landfire_mosaic --product FBFM13 --region CONUS --build

--list reports what is staged. Without --build the run is a dry run: it
validates the stack and estimates the mosaic from raster overviews without
writing anything. --build writes the mosaic and, with --geo-vars, a geo_vars
JSON pointing at it.
"""
from __future__ import absolute_import
from __future__ import print_function

import argparse
import glob
import json
import os
import os.path as osp
import re
import sys
import time
from collections import namedtuple

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.windows import Window

# Default location of the staged LANDFIRE releases.
STAGING_DIR = '/data/jhaley/wrfxpy/landfire'

# Values treated as absent data. See the module docstring: only 32767 is
# declared in the GeoTIFF metadata, so both must be listed explicitly.
FILL_VALUES = (32767, -9999)

# Single declared nodata for the mosaic we write.
OUT_NODATA = 32767

# Sidecar recording which releases went into a mosaic. Without it, a forecast's
# fuels are traceable only to 'some mosaic', because the raster itself carries
# no record of its inputs.
MANIFEST_SUFFIX = '.manifest.json'

# LF<year>_<product>[_<version>]_<region>[_<suffix>]
PRODUCT_RE = re.compile(
    r'^LF(?P<year>\d{4})'
    r'_(?P<product>[A-Za-z0-9]+?)'
    r'(?:_(?P<version>\d+))?'
    r'_(?P<region>CONUS|AK|HI|PRVI)'
    r'(?:_(?P<suffix>.+))?$'
)

Product = namedtuple('Product',
                     'name year product version region suffix path tif archive')


def parse_name(name):
    """Parse a LANDFIRE release directory or archive name, or return None."""
    m = PRODUCT_RE.match(name)
    if not m:
        return None
    d = m.groupdict()
    return {
        'year': int(d['year']),
        'product': d['product'],
        'version': int(d['version']) if d['version'] else 0,
        'region': d['region'],
        'suffix': d['suffix'] or '',
    }


def find_tif(product_dir):
    """Return the single GeoTIFF inside a release's Tif/ subdirectory."""
    tifs = sorted(glob.glob(osp.join(product_dir, 'Tif', '*.tif')))
    if len(tifs) == 1:
        return tifs[0]
    return None


def discover(staging_dir=STAGING_DIR):
    """
    Find every LANDFIRE release staged under staging_dir.

    Returns (products, problems). Extracted directories with a readable
    GeoTIFF become Product entries; archives with no extracted counterpart,
    and names that do not parse, are reported in problems so nothing is
    silently ignored.
    """
    products, problems = [], []
    entries = sorted(os.listdir(staging_dir))
    extracted = set(e for e in entries if osp.isdir(osp.join(staging_dir, e)))

    for entry in entries:
        path = osp.join(staging_dir, entry)
        is_archive = entry.endswith('.zip')
        base = entry[:-4] if is_archive else entry

        fields = parse_name(base)
        if fields is None:
            if is_archive or osp.isdir(path):
                problems.append(('unparsed name', entry))
            continue

        if is_archive:
            if base not in extracted:
                problems.append(('archive not extracted', entry))
            continue

        tif = find_tif(path)
        if tif is None:
            problems.append(('no single Tif/*.tif', entry))
            continue

        products.append(Product(name=base, path=path, tif=tif,
                                archive=osp.exists(path + '.zip'), **fields))
    return products, problems


def select_stack(products, product, region):
    """
    Return matching releases, newest first.

    Ordered by (year, version, is a regional refresh). A same-year suffixed
    release such as LF2023_FBFM13_240_CONUS_sw_update is a later refresh of
    its base, so it is preferred over that base.
    """
    matches = [p for p in products
               if p.product.upper() == product.upper()
               and p.region.upper() == region.upper()]
    return sorted(matches,
                  key=lambda p: (p.year, p.version, 1 if p.suffix else 0),
                  reverse=True)


def describe_grid(path):
    """Read the grid description needed to check that layers are compatible."""
    with rasterio.open(path) as src:
        return {
            'crs': src.crs.to_string() if src.crs else None,
            'res': (src.res[0], src.res[1]),
            'left': src.bounds.left,
            'top': src.bounds.top,
            'right': src.bounds.right,
            'bottom': src.bounds.bottom,
            'width': src.width,
            'height': src.height,
            'dtype': src.dtypes[0],
            'nodata': src.nodata,
        }


def validate_stack(stack):
    """
    Check the stack can be mosaicked without resampling.

    Requires one CRS, one pixel size, and origins that differ by a whole
    number of pixels. Returns (grids, errors); a non-empty errors list means
    the caller must stop, because pasting misaligned layers would put fuel
    categories on the wrong ground.
    """
    grids = [describe_grid(p.tif) for p in stack]
    errors = []

    crs_set = set(g['crs'] for g in grids)
    if len(crs_set) > 1:
        errors.append('layers disagree on CRS: %s' % sorted(crs_set))

    res_set = set(g['res'] for g in grids)
    if len(res_set) > 1:
        errors.append('layers disagree on pixel size: %s' % sorted(res_set))

    if not errors:
        xres, yres = grids[0]['res']
        x0 = min(g['left'] for g in grids)
        y0 = max(g['top'] for g in grids)
        for p, g in zip(stack, grids):
            dx = (g['left'] - x0) / xres
            dy = (y0 - g['top']) / yres
            if abs(dx - round(dx)) > 1e-6 or abs(dy - round(dy)) > 1e-6:
                errors.append('%s is off the shared lattice by (%.4f, %.4f) px'
                              % (p.name, dx, dy))

    dtypes = set(g['dtype'] for g in grids)
    if len(dtypes) > 1:
        errors.append('layers disagree on dtype: %s' % sorted(dtypes))

    return grids, errors


def union_grid(grids):
    """Return the output grid covering every layer, on the shared lattice."""
    xres, yres = grids[0]['res']
    left = min(g['left'] for g in grids)
    top = max(g['top'] for g in grids)
    right = max(g['right'] for g in grids)
    bottom = min(g['bottom'] for g in grids)
    width = int(round((right - left) / xres))
    height = int(round((top - bottom) / yres))
    transform = rasterio.transform.from_origin(left, top, xres, yres)
    return {'left': left, 'top': top, 'right': right, 'bottom': bottom,
            'width': width, 'height': height, 'transform': transform,
            'res': (xres, yres)}


def gap_mask(block):
    """True where a block holds no fuel category. See FILL_VALUES."""
    return np.isin(block, FILL_VALUES)


def estimate(stack, grids, shape=(800, 1200)):
    """
    Estimate layer contributions from raster overviews, without writing.

    Decimated reads use the pyramid levels, so this is fast but approximate:
    it answers 'roughly how much does each release contribute' rather than
    giving exact pixel counts.
    """
    out = union_grid(grids)
    filled = np.zeros(shape, dtype=bool)
    report = []
    for p, g in zip(stack, grids):
        with rasterio.open(p.tif) as src:
            block = src.read(1, out_shape=shape, boundless=False)
        # Place the layer in the union frame by scaling its offset. Coarse,
        # which is why this is an estimate and not the built product.
        frac_x = (g['left'] - out['left']) / (out['right'] - out['left'])
        frac_y = (out['top'] - g['top']) / (out['top'] - out['bottom'])
        col0 = int(round(frac_x * shape[1]))
        row0 = int(round(frac_y * shape[0]))
        h = int(round(shape[0] * (g['top'] - g['bottom'])
                      / (out['top'] - out['bottom'])))
        w = int(round(shape[1] * (g['right'] - g['left'])
                      / (out['right'] - out['left'])))
        h, w = min(h, shape[0] - row0), min(w, shape[1] - col0)
        if h <= 0 or w <= 0:
            report.append((p.name, 0.0))
            continue
        sub = block[:h, :w]
        target = (slice(row0, row0 + h), slice(col0, col0 + w))
        usable = ~gap_mask(sub) & ~filled[target]
        report.append((p.name, 100.0 * usable.sum() / filled.size))
        filled[target] |= usable
    return report, 100.0 * filled.mean()


def build(stack, grids, out_path, strip_rows=256):
    """
    Write the mosaic, newest release first, filling gaps from older ones.

    Processes the output in horizontal strips so the full raster (tens of GB
    uncompressed) is never held in memory. Writes to a temporary name and
    renames into place, so an interrupted run cannot leave a truncated raster
    where something later reads fuels.
    """
    out = union_grid(grids)
    profile = {
        'driver': 'GTiff', 'dtype': 'int16', 'count': 1,
        'width': out['width'], 'height': out['height'],
        'crs': grids[0]['crs'], 'transform': out['transform'],
        'nodata': OUT_NODATA,
        'compress': 'lzw', 'BIGTIFF': 'IF_SAFER',
    }
    # Tiling needs a raster at least one block across; the small regional
    # products (and tests) can be smaller than the block size.
    if out['width'] >= 256 and out['height'] >= 256:
        profile.update({'tiled': True, 'blockxsize': 256, 'blockysize': 256})
    contributed = dict((p.name, 0) for p in stack)
    gaps = 0
    tmp_path = '%s.%d.tmp' % (out_path, os.getpid())

    srcs = [rasterio.open(p.tif) for p in stack]
    try:
        with rasterio.open(tmp_path, 'w', **profile) as dst:
            for row0 in range(0, out['height'], strip_rows):
                rows = min(strip_rows, out['height'] - row0)
                strip = np.full((rows, out['width']), OUT_NODATA, dtype='int16')
                filled = np.zeros((rows, out['width']), dtype=bool)

                for p, g, src in zip(stack, grids, srcs):
                    xres, yres = out['res']
                    off_col = int(round((g['left'] - out['left']) / xres))
                    off_row = int(round((out['top'] - g['top']) / yres))
                    # rows of this strip that fall inside this layer
                    src_row0 = row0 - off_row
                    r_lo = max(0, -src_row0)
                    r_hi = min(rows, g['height'] - src_row0)
                    if r_hi <= r_lo:
                        continue
                    window = Window(0, src_row0 + r_lo, g['width'], r_hi - r_lo)
                    block = src.read(1, window=window)
                    target = (slice(r_lo, r_hi),
                              slice(off_col, off_col + g['width']))
                    usable = ~gap_mask(block) & ~filled[target]
                    strip[target] = np.where(usable, block, strip[target])
                    filled[target] |= usable
                    contributed[p.name] += int(usable.sum())

                gaps += int((~filled).sum())
                dst.write(strip, 1, window=Window(0, row0, out['width'], rows))

            # Only levels that leave a usefully sized image. Nearest
            # neighbour because fuel model codes are categorical: averaging
            # them would invent categories that do not exist.
            levels = [f for f in (2, 4, 8, 16, 32, 64, 128)
                      if min(out['width'], out['height']) // f >= 128]
            if levels:
                print('   building overviews %s (nearest, categorical data)'
                      % levels)
                dst.build_overviews(levels, Resampling.nearest)
    finally:
        for src in srcs:
            src.close()

    os.replace(tmp_path, out_path)
    total = out['width'] * out['height']
    return contributed, gaps, total


def write_manifest(out_path, product, region, stack, grids, contributed,
                   gaps, total, note=None):
    """
    Record which releases went into a mosaic, next to the mosaic itself.

    Input files are identified by size and mtime rather than a checksum:
    these rasters are gigabytes each and hashing them would add minutes to
    every build for a guarantee that size plus mtime already gives in
    practice. If you need certainty that an input has not changed, hash it
    yourself and compare.
    """
    out = union_grid(grids)
    layers = []
    for i, (p, g) in enumerate(zip(stack, grids)):
        st = os.stat(p.tif)
        layers.append({
            'rank': i + 1,
            'name': p.name,
            'year': p.year,
            'version': p.version,
            'suffix': p.suffix,
            'tif': p.tif,
            'tif_bytes': st.st_size,
            'tif_mtime': time.strftime('%Y-%m-%dT%H:%M:%S',
                                       time.localtime(st.st_mtime)),
            'pixels_contributed': contributed[p.name],
            'percent_of_output': round(100.0 * contributed[p.name] / total, 4),
        })
    manifest = {
        'built_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'built_by': 'ingest/landfire_mosaic.py',
        'product': product,
        'region': region,
        'output': out_path,
        'output_grid': dict(crs=grids[0]['crs'],
                            pixel_size=out['res'][0],
                            width=out['width'], height=out['height'],
                            left=out['left'], top=out['top'],
                            right=out['right'], bottom=out['bottom'],
                            nodata=OUT_NODATA),
        'gaps_treated_as_missing': list(FILL_VALUES),
        'layers': layers,
        'pixels_total': total,
        'pixels_no_data_in_any_layer': gaps,
        'percent_with_fuel': round(100.0 * (total - gaps) / total, 4),
    }
    if note:
        manifest['note'] = note
    path = out_path + MANIFEST_SUFFIX
    with open(path, 'w') as fh:
        json.dump(manifest, fh, indent=2, separators=(',', ': '))
        fh.write('\n')
    return path


def write_geo_vars(out_json, nfuel_cat, zsf):
    """Write a geo_vars JSON pointing NFUEL_CAT at the mosaic."""
    with open(out_json, 'w') as fh:
        json.dump({'NFUEL_CAT': nfuel_cat, 'ZSF': zsf}, fh,
                  indent=2, separators=(',', ': '))
        fh.write('\n')


def report_products(products, problems):
    print('Staged LANDFIRE releases:')
    fmt = '   {:38s} {:>6s} {:>8s} {:>6s} {:>10s}'
    print(fmt.format('name', 'year', 'product', 'region', 'refresh'))
    for p in sorted(products, key=lambda p: (p.product, p.region, -p.year)):
        print(fmt.format(p.name, str(p.year), p.product, p.region,
                         p.suffix or '-'))
    if problems:
        print('\nNot usable as a layer:')
        for why, entry in problems:
            print('   {:26s} {}'.format(why, entry))


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Build a LANDFIRE fuels mosaic, newest release first, '
                    'filling gaps from older releases.')
    parser.add_argument('--staging', default=STAGING_DIR,
                        help='directory of staged LANDFIRE releases')
    parser.add_argument('--product', default='FBFM13',
                        help='product to mosaic, e.g. FBFM13 (default) or FBFM40')
    parser.add_argument('--region', default='CONUS',
                        help='region to mosaic: CONUS (default), AK, HI, PRVI')
    parser.add_argument('--list', action='store_true',
                        help='report what is staged and exit')
    parser.add_argument('--build', action='store_true',
                        help='write the mosaic; without this the run is a dry run')
    parser.add_argument('--out', default=None,
                        help='output GeoTIFF (default: <staging>/LFmosaic_<product>_<region>/Tif/...)')
    parser.add_argument('--geo-vars', default=None,
                        help='also write this geo_vars JSON pointing at the mosaic')
    parser.add_argument('--zsf', default=None,
                        help='ZSF elevation GeoTIFF to record in --geo-vars')
    args = parser.parse_args(argv)

    products, problems = discover(args.staging)
    if args.list:
        report_products(products, problems)
        return 0

    stack = select_stack(products, args.product, args.region)
    if not stack:
        print('No staged releases match product=%s region=%s'
              % (args.product, args.region))
        report_products(products, problems)
        return 1

    print('Layer stack for %s %s, newest first:' % (args.product, args.region))
    for i, p in enumerate(stack):
        print('   %d. %s' % (i + 1, p.name))

    grids, errors = validate_stack(stack)
    if errors:
        print('\nStack is not mosaickable without resampling:')
        for e in errors:
            print('   ' + e)
        print('Refusing to proceed: pasting misaligned layers would place '
              'fuel categories on the wrong ground.')
        return 1
    print('\nGrid check: one CRS (%s), pixel size %g, all layers on the '
          'shared lattice.' % (grids[0]['crs'], grids[0]['res'][0]))

    out = union_grid(grids)
    print('Output grid: %d x %d px, %.0f..%.0f x, %.0f..%.0f y'
          % (out['width'], out['height'], out['left'], out['right'],
             out['bottom'], out['top']))

    if not args.build:
        print('\nDry run. Estimated contribution by layer (from overviews, '
              'approximate):')
        report, covered = estimate(stack, grids)
        for name, pct in report:
            print('   %-38s %5.1f%% of the output grid' % (name, pct))
        print('   %-38s %5.1f%%' % ('total with fuel data', covered))
        print('   %-38s %5.1f%%' % ('remaining gaps', 100.0 - covered))
        print('\nRe-run with --build to write the mosaic.')
        return 0

    out_path = args.out
    if out_path is None:
        out_dir = osp.join(args.staging,
                           'LFmosaic_%s_%s' % (args.product, args.region), 'Tif')
        os.makedirs(out_dir, exist_ok=True)
        out_path = osp.join(out_dir,
                            'LFmosaic_%s_%s.tif' % (args.product, args.region))

    print('\nWriting %s' % out_path)
    contributed, gaps, total = build(stack, grids, out_path)
    print('\nPixel contribution by layer:')
    for p in stack:
        n = contributed[p.name]
        print('   %-38s %14d px  %5.1f%%' % (p.name, n, 100.0 * n / total))
    print('   %-38s %14d px  %5.1f%%' % ('no data in any layer', gaps,
                                         100.0 * gaps / total))
    print('\nWrote %s' % out_path)
    manifest = write_manifest(out_path, args.product, args.region, stack,
                              grids, contributed, gaps, total)
    print('Wrote %s' % manifest)

    if args.geo_vars:
        if not args.zsf:
            print('--geo-vars needs --zsf, the elevation GeoTIFF to record')
            return 1
        write_geo_vars(args.geo_vars, out_path, args.zsf)
        print('Wrote %s' % args.geo_vars)
    return 0


if __name__ == '__main__':
    sys.exit(main())
