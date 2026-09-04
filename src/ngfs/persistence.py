"""
Storing and retrieving the system state, and writing each run's output files.

State is one gzip-compressed pickle of the whole ngfs_day object per run, in the
configured ngfs_directory. get_old_incidents rebuilds a run's starting point
from the last seven days of those files: the union of every started_inc_ids
list, the incident objects, and the newest run's detection frame.

Things that are easy to get wrong here:

  * mtime is load-bearing. Both the seven-day filter and the choice of "newest
    state" sort by mtime, so any tool that rewrites a state file must preserve
    its timestamp, or it promotes stale state to newest and corrupts the ledger.
  * Pickles record different class paths depending on their writer. The monolith
    runs as a script, so its pickles name __main__.ngfs_day, while
    package-written ones name ngfs.ngfs_day.ngfs_day; reading the former from
    the package raises AttributeError. pickle_writer determines ownership from
    the file's first 128 bytes rather than from its name -- do not rely on the
    '_testing' filename suffix.
  * ngfs_start.py globs only *.pkl, so compressing state in the monolith's own
    directory would hide it and cause fires to be forecast a second time.
    compress_state_pickles refuses such directories.

Also holds the run's output writers -- detection_summary, print_base_map and
save_incident_text -- which ngfs_day delegates to.

Run as a module to compress a backlog of uncompressed pickles. It verifies each
compressed copy by SHA-256 of the decompressed bytes and deletes the original
only on a match; without --delete it reports and changes nothing, and files
newer than --min-age-hours are never touched, so it is safe to run while the
loop is live:

    python -m ngfs.persistence <directory> [--delete] [--limit N]
"""
import glob
import gzip
import hashlib
import time
import os
import pandas as pd
import numpy as np
from mpl_toolkits.basemap import Basemap
#from ngfs import config_manager
import matplotlib.pyplot as plt
from PIL import Image
from datetime import timedelta, datetime
from ngfs import constants as cons

####### Functions  #######
def state_pickle_files(ngfs_directory):
    """
    Returns the saved state pickles in `ngfs_directory`, oldest first.

    Accepts compressed and uncompressed pickles alike; pandas infers the codec
    from the file extension when reading them back.
    """
    paths = []
    for pattern in cons.PICKLE_PATTERNS:
        paths.extend(glob.glob(os.path.join(ngfs_directory, pattern)))
    return sorted(paths, key=os.path.getmtime)


#Pickles record the module path of the class they hold. Files written by this
#package name ngfs.ngfs_day; files written by the precursor src/ngfs_start.py
#name __main__, because that script defines its classes in the script namespace.
#That distinction matters for compression: ngfs_start.py globs only '*.pkl', so
#compressing a file it owns would hide that state from it entirely.
PACKAGE_PICKLE_MARKER = b'ngfs.ngfs_day'
MONOLITH_PICKLE_MARKER = b'__main__'


def pickle_writer(path):
    """
    Identifies which program wrote a state pickle, by reading its class path.

    Returns 'package' for pickles this package can read back, 'monolith' for
    ones belonging to src/ngfs_start.py, or 'unknown' when neither marker is
    found in the header.
    """
    with open(path, 'rb') as handle:
        header = handle.read(128)
    if PACKAGE_PICKLE_MARKER in header:
        return 'package'
    if MONOLITH_PICKLE_MARKER in header:
        return 'monolith'
    return 'unknown'


def file_digest(path, opener=open):
    """Returns the SHA-256 of a file's contents, read through `opener`."""
    digest = hashlib.sha256()
    with opener(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def compress_state_pickles(ngfs_directory, min_age_hours=6, delete=False,
                           limit=None):
    """
    Compresses previously saved uncompressed state pickles in place.

    Each file is gzipped at the byte level rather than being unpickled and
    re-saved: that preserves the contents exactly, needs no class to resolve,
    and produces the same gzip-wrapped pickle that pandas reads back. The
    compressed copy is verified by comparing the SHA-256 of its decompressed
    bytes against the original, and the original is removed only when that
    matches and `delete` is set.

    Refuses any directory holding pickles written by src/ngfs_start.py, whose
    reader globs only '*.pkl' and would lose sight of a compressed file.

    `min_age_hours` skips recently written state so a run in progress is never
    touched; `limit` caps how many files are processed, for a cautious first
    pass. Files that already have a compressed counterpart are skipped, so an
    interrupted sweep can simply be run again.
    """
    uncompressed = sorted(glob.glob(os.path.join(ngfs_directory, '*.pkl')),
                          key=os.path.getmtime)
    if not uncompressed:
        print(f'No uncompressed state pickles found in {ngfs_directory}')
        return {}

    #Guard on the newest pickle, since that is the one a reader picks up: if it
    #belongs to the monolith then this is a directory the monolith serves, and
    #compressing anything here risks hiding state from a reader that globs only
    #'*.pkl'. Individual stray files are skipped in the loop instead.
    newest_owner = pickle_writer(uncompressed[-1])
    if newest_owner != 'package':
        raise ValueError(
            f'Refusing to compress {ngfs_directory}: its most recent pickle was '
            f'written by {newest_owner} (src/ngfs_start.py names its classes '
            f'__main__). That reader globs only "*.pkl", so compressing state '
            f'it owns would hide it and cause re-forecasting. Point this at the '
            f'directory used by ngfs_start_2.py instead.')

    cutoff = time.time() - min_age_hours * 3600
    summary = {'compressed': 0, 'skipped': 0, 'failed': [], 'deleted': 0,
               'bytes_before': 0, 'bytes_after': 0, 'bytes_freed': 0}

    for path in uncompressed:
        target = path + '.gz'
        if os.path.getmtime(path) > cutoff:
            print(f'\tSkipping {os.path.basename(path)}: newer than '
                  f'{min_age_hours} h')
            summary['skipped'] += 1
            continue
        owner = pickle_writer(path)
        if owner != 'package':
            print(f'\tSkipping {os.path.basename(path)}: written by {owner}')
            summary['skipped'] += 1
            continue
        if os.path.exists(target):
            #compressed by an earlier pass, so this run can still finish the job:
            #re-verify that copy against the original before discarding it, which
            #also lets an interrupted sweep be cleaned up by running again
            if delete:
                original_size = os.path.getsize(path)
                if file_digest(target, opener=gzip.open) == file_digest(path):
                    os.remove(path)
                    summary['deleted'] += 1
                    summary['bytes_freed'] += original_size
                    print(f'\tRemoved {os.path.basename(path)}: verified '
                          f'against existing compressed copy')
                else:
                    print(f'\tFAILED {os.path.basename(path)}: existing '
                          f'{os.path.basename(target)} does not match; '
                          f'original kept')
                    summary['failed'].append(path)
            summary['skipped'] += 1
            continue
        if limit is not None and summary['compressed'] >= limit:
            break

        original_size = os.path.getsize(path)
        stat = os.stat(path)
        tmp_target = f'{target}.{os.getpid()}.tmp'
        try:
            with open(path, 'rb') as source:
                with gzip.open(tmp_target, 'wb',
                               compresslevel=cons.PICKLE_COMPRESSION['compresslevel']) as sink:
                    for chunk in iter(lambda: source.read(1024 * 1024), b''):
                        sink.write(chunk)
            #keep the original timestamps: both readers sort candidates by mtime
            #and apply an age filter, so a fresh mtime would promote old state
            os.utime(tmp_target, (stat.st_atime, stat.st_mtime))
            if file_digest(tmp_target, opener=gzip.open) != file_digest(path):
                raise ValueError('decompressed contents differ from original')
            os.replace(tmp_target, target)
        except Exception as exc:
            if os.path.exists(tmp_target):
                os.remove(tmp_target)
            print(f'\tFAILED {os.path.basename(path)}: {exc!r}')
            summary['failed'].append(path)
            continue

        compressed_size = os.path.getsize(target)
        summary['compressed'] += 1
        summary['bytes_before'] += original_size
        summary['bytes_after'] += compressed_size
        print(f'\t{os.path.basename(path)}: '
              f'{original_size / 1048576.0:.1f} MB -> '
              f'{compressed_size / 1048576.0:.1f} MB '
              f'({original_size / float(compressed_size):.1f}x)')

        if delete:
            os.remove(path)
            summary['deleted'] += 1
            summary['bytes_freed'] += original_size

    print(f'Compressed {summary["compressed"]} pickle(s), '
          f'skipped {summary["skipped"]}, failed {len(summary["failed"])}')
    if summary['compressed']:
        print(f'\tcompressed {summary["bytes_before"] / 1073741824.0:.2f} GB '
              f'down to {summary["bytes_after"] / 1073741824.0:.2f} GB')
    if delete:
        print(f'\tremoved {summary["deleted"]} verified original(s), freeing '
              f'{summary["bytes_freed"] / 1073741824.0:.2f} GB')
    else:
        pending = summary['bytes_before'] / 1073741824.0
        print(f'\toriginals kept; re-running with delete=True would verify and '
              f'remove them, freeing {pending:.2f} GB')
    return summary


def get_old_incidents(ngfs_directory):
    """
    Looks at older, saved ngfs_days objects and finds incidents that have already been processed or forecasted.
    """
    from ngfs.ngfs_day import ngfs_day
    print('Reading previous pickle file(s)')
    
    # Get and sort pickle files by modification time
    full_pick_list = state_pickle_files(ngfs_directory)
    current_time = time.time()

    # Filter out old pickle files not generated by specific CSV files
    pick_list = [
        i for i in full_pick_list
        if 'GOES' not in i and (current_time - os.path.getmtime(i)) / 3600 < 24 * 7
    ]

    #take most recent pickle file, even if it is older than 7 days
    if not len(pick_list):
      pick_list = [i for i in full_pick_list[-1:] if 'GOES' not in i]

    #no previous state at all: a first run, or a fresh ngfs_directory
    if not len(pick_list):
        print(f'No previous state pickles found in {ngfs_directory}; '
              f'starting with no history')
        return [], [], None

    print(f'Number of pickle files to possibly look at: {len(pick_list)}')

    # Reverse the list to get the newest files first
    pick_list.reverse()

    old_ngfs_incidents = []
    started_inc_ids = []
    latest_day = None

    for i in pick_list:
        print(f'\tReading {i}')
        try:
            df = pd.read_pickle(i)
        except Exception as exc:
            #a pickle can be unreadable because it is truncated, or because it
            #belongs to src/ngfs_start.py, whose classes are named __main__ and
            #cannot be resolved here. Either way, try the next candidate rather
            #than losing the whole ledger.
            print(f'\tCould not read {i}: {exc!r}')
            continue
        latest_day = df
        old_ngfs_incidents.extend(df.incidents)

        if hasattr(df, 'started_inc_ids'):
            started_inc_ids.extend(df.started_inc_ids)
            print(f'\tTracking {len(started_inc_ids)} old incident ID strings')
            break
        else:
            print('\tPickle file does not have started_inc_ids')
        print(len(started_inc_ids))

    if latest_day is None:
        print(f'None of the {len(pick_list)} candidate pickle(s) could be read; '
              f'starting with no history')
        return [], [], None

    print(f'\tTimestamp of loaded ngfs_day object: {latest_day.timestamp} '
          f'with data size: {len(latest_day.data)}')
    # Process and print old incidents
    print('Found the previous incidents:')
    for inc in old_ngfs_incidents:
        if inc.incident_id_string in started_inc_ids or inc.started:
            inc.started = True
            if inc.incident_id_string not in started_inc_ids:
                started_inc_ids.append(inc.incident_id_string)
        print(f'\t{inc.incident_id_string} {inc.incident_name} Started = {inc.started}')


    return list(set(started_inc_ids)), old_ngfs_incidents, latest_day

def detection_summary(ngfs_day, hours=24):
    """
    Print summary information about detections from the various sources
    over the previous period of length `hours`.
    """
    print(f'Length of dat set: {len(ngfs_day.data)}')
    min_time = ngfs_day.timestamp - timedelta(hours=hours)
    # Determine time column
    for c in ngfs_day.data.columns:
        if 'time' in c:
            time_col = c
    # Filter once
    df = ngfs_day.data# [ngfs_day.data[time_col] > min_time].copy()
    sat_list = list(df.satellite.unique())
    det_summary = pd.DataFrame()
    for sat in sat_list:
        tmp = dict()
        print(sat)
        data = df[df.satellite == sat ].copy()
        print(f'Length of {sat} data set: {len(data)}')
        tmp['satellite_name'] = sat
        tmp['start_timestamp'] = data[time_col].min()
        tmp['end_timestamp'] = ngfs_day.timestamp
        tmp['known_incidents'] = len(data.known_incident_id.unique())
        tmp['total_detections'] = len(data)
        if "type_description" in data.columns:
            tmp['known_wildland_fire'] = len(data[data.type_description == 'Known Wildland Fire Incident'])
            tmp['possible_wildland_fire'] = len(data[data.type_description == 'Possible Wildland Fire'])
        else:
            tmp['known_wildland_fire'] = len(data[data.type == 1])
            tmp['possible_wildland_fire'] = len(data[data.type == 0])
        tmp['total_frp'] = np.sum(data.frp)
        for k in tmp.keys():
            print(f'\t{k}\t{tmp[k]}')
        #
        det_summary = det_summary.append(pd.DataFrame(tmp,index=[tmp['satellite_name']]))
    # Save
    det_summary.to_csv("ngfs/detection_summary_testing.csv", index=False)
    det_summary.to_csv(f"ngfs/detection_summary_{ngfs_day.date_str}_testing.csv", index=False)
    #return det_summary
def print_base_map(ngfs_day):
    #prints map of locations of incidents within ngfs_day object
    print('Starting map making')
    #find, new, started, and ongoing incidents
    ngfs_day.set_new()
    ngfs_day.set_started()
    ngfs_day.set_ongoing()

    #projects coordinates onto the map space
    def get_coords(latlons, mask, map):
        #map is a Basemap object, this returns projected coordinates
        #mask is boolean array
        filtered = latlons[mask]
        return map(filtered[:, 1], filtered[:, 0])
    #scatters locations onto map
    def scatter_if_any(x, y, size, label, map):
        #scatters projected coordinates and adds legend text to plot
        if len(x) > 0:
            map.scatter(x, y, s=size, label=label, edgecolors='black')
    #CONUS map
    m = Basemap(llcrnrlon=-119, llcrnrlat=22, urcrnrlon=-64, urcrnrlat=49, projection='lcc', lat_1=33, lat_2=45, lon_0=-95)
    m.latlon = True
    m.shadedrelief()
    m.drawstates()
    m.drawcountries()
    m.drawcoastlines()

    #nx2 matrix with lat,lon as columns
    latlons = ngfs_day.incident_ign_latlons()
    #print('latlons',latlons)

    #tuple with incident status, icon size, legend text label for plotting on map
    incident_data = [
        (ngfs_day.ongoing, 15, 'Ongoing Incident'),
        (np.logical_and(~ngfs_day.started,ngfs_day.new), 30, 'New Incident'),
        (np.logical_and(ngfs_day.started,ngfs_day.new), 30, 'New, Forecast Started')
    ]

    #scatter the CONUS incidents
    for data, size, label in incident_data:
        #data is a boolean array
        x, y = get_coords(latlons, data, m)
        scatter_if_any(x, y, size, label, m)

    plt.legend(loc="lower left")
    plt.title(f'NIFC Fire Incidents \n{ngfs_day.date_str}\n {ngfs_day.sat_name}')

    print('Saving ', ngfs_day.map_save_str)
    plt.savefig(ngfs_day.map_save_str, bbox_inches='tight')  #always saved
    plt.cla()

    #if there are fires in Alaska, make a map and join it to the existing CONUS map
    if max(latlons[:, 0]) > 54.0:
        n = Basemap(llcrnrlon=-164, llcrnrlat=54, urcrnrlon=-130, urcrnrlat=73, projection='lcc', lat_1=63, lat_2=68, lon_0=-151)
        n.latlon = True
        n.shadedrelief()
        n.drawstates()
        n.drawcountries()
        n.drawcoastlines()
        #scatter the Alaska incidents on Alaskan map
        for data, size, label in incident_data:
            x, y = get_coords(latlons, data, n)
            scatter_if_any(x, y, size, label, n)

        plt.legend(loc="lower right")
        plt.title(f'NIFC Alaska Fire Incidents \n{ngfs_day.date_str}\n GOES-18')

        sv_str = ngfs_day.map_save_str.replace('ngfs/', 'ngfs/Alaska_')
        print('Saving ', sv_str)
        plt.savefig(sv_str, bbox_inches='tight')

        conus_str = ngfs_day.map_save_str.replace('ngfs/', 'ngfs/CONUS_')
        os.system(f'cp {ngfs_day.map_save_str} {conus_str}')

        #join alaska and CONUS image after waiting for file to write
        time.sleep(20)
        img0 = Image.open(ngfs_day.map_save_str)
        img1 = Image.open(sv_str)

        height = img0.size[1]
        width = int(np.round(height * img1.size[0] / img1.size[1]))
        img1 = img1.resize((width, height), Image.LANCZOS)
        img1.save(sv_str)

        image_new = Image.new("RGB", (img0.size[0] + img1.size[0], height), "white")
        image_new.paste(img1, (0, 0))
        image_new.paste(img0, (width, 0))
        image_new.save(ngfs_day.map_save_str)

def save_incident_text(ngfs_day):
      """
      Saves the estimated ignition point information for all new incidents as a CSV file.
      """
      #Needs to add
      # additional fields to the csv file such as viirs pixel location and assumed ignition time
      # use pd.to_csv() function
      ign_pix = pd.DataFrame()
      new_ign_lat = []
      new_ign_lon = []
      new_ign_time = []
      viirs_pixel = []

      for inc in ngfs_day.incidents:
         #print(inc.name,inc.new)
         if (inc.new & inc.started):
            #print(inc.ignition_pixel)
            ign_pix = ign_pix.append(inc.ignition_pixel)
            #for adding columns to the csv
            try:
               new_ign_lat = np.append(new_ign_lat,inc.new_ign_latlon[0])
               new_ign_lon = np.append(new_ign_lon,inc.new_ign_latlon[1])
               new_ign_time.append(inc.ign_utc)
               viirs_pix = True
               #viirs_pixel.append(True)
            except AttributeError:
               new_ign_lat = np.append(new_ign_lat,inc.ign_latlon[0])
               new_ign_lon = np.append(new_ign_lon,inc.ign_latlon[1])
               new_ign_time.append(inc.ign_utc)
               viirs_pix = False
            if hasattr(inc,'viirs_ignition_pixel'):
               viirs_pixel.append(inc.viirs_ignition_pixel)
            else:
               viirs_pixel.append(viirs_pix)
            
      #print('length of new vars',len(new_ign_lat),len(new_ign_lon),len(new_ign_time))
      #print('dataframe shape',ign_pix.shape)
      ign_pix['forecast_ign_lat'] = new_ign_lat
      ign_pix['forecast_ign_lon'] = new_ign_lon
      ign_pix['forecast_ign_UTC'] = new_ign_time
      ign_pix['viirs_pixel_ign'] = viirs_pixel
      print(ign_pix)
      #only save non-empty dataframe
      if not ign_pix.empty:
         time_str = time.strftime("%H_%M", time.localtime())
         csv_save_str = f'ngfs/forecast_ignition_pixels_{ngfs_day.date_str}_{time_str}.csv'
         ign_pix.to_csv(csv_save_str, index=False)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Compress previously saved uncompressed state pickles. '
                    'Reports what it would do unless --delete is given.')
    parser.add_argument('directory',
                        help='directory of state pickles, e.g. ngfs')
    parser.add_argument('--min-age-hours', type=float, default=6,
                        help='leave pickles newer than this alone (default 6)')
    parser.add_argument('--limit', type=int, default=None,
                        help='stop after this many files, for a cautious pass')
    parser.add_argument('--delete', action='store_true',
                        help='remove each original once its compressed copy '
                             'has been verified byte-for-byte')
    args = parser.parse_args()

    compress_state_pickles(args.directory, min_age_hours=args.min_age_hours,
                           delete=args.delete, limit=args.limit)
