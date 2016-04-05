# Copyright (C) 2013-2016 Martin Vejmelka, UC Denver
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


from datetime import datetime, timedelta
import pytz
import os
import os.path as osp
import glob
import numpy as np
import math


def ensure_dir(path):
    """
    Ensure all directories in path exist, for convenience return path itself.

    :param path: the path whose directories should exist
    :return: the path back for convenience
    """
    path_dir = osp.dirname(path)
    if not osp.exists(path_dir):
        os.makedirs(path_dir)
    return path


def make_dir(dir):
    """
    Create a directory if it does not exist. Creates any intermediate dirs as necessary.

    :param dir: the directory to be created
    """
    if not osp.exists(dir):
        os.makedirs(dir)


def symlink_unless_exists(link_tgt, link_loc):
    """
    Create a symlink at link_loc pointing to link_tgt unless file already exists.
    
    :param link_tgt: link target
    :param link_loc: link location
    """
    if not osp.lexists(link_loc):
        os.symlink(link_tgt, link_loc)


def esmf_to_utc(esmf):
    """
    Converts an ESMF datetime into a UTC datetime.

    :param esmf: date & time from YYYY-MM-DD_hh:mm:ss (ESMF) format
    :return: a python datetime with UTC timezone
    """
    # ESMF date: YYYY-MM-DD_hh:mm:ss
    year, mon, day = int(esmf[0:4]), int(esmf[5:7]), int(esmf[8:10])
    hour, min, sec = int(esmf[11:13]), int(esmf[14:16]), int(esmf[17:19])
    return datetime(year, mon, day, hour, min, sec, tzinfo=pytz.utc)


def utc_to_esmf(utc):
    """
    Converts a UTC datetime into ESMF format.

    :param utc: python UTC datetime
    :return: a string in ESMF format
    """
    return "%04d-%02d-%02d_%02d:%02d:%02d" % (utc.year, utc.month, utc.day, utc.hour, utc.minute, utc.second)


def compute_fc_hours(start_utc, end_utc):
    """
    Compute the number of forecast hours.  Rounds up in case of incomplete hours.

    :param start_utc: start time of forecast
    :param end_utc: end time of forecast
    :return: number of forecast hours (rounded up to nearest hour)
    """
    delta = end_utc - start_utc
    return delta.days * 24 + int(delta.seconds + 3600 - 1) / 3600


def symlink_matching_files(tgt_dir, src_dir, glob_pattern):
    """
    Retrieves all files in directory src_dir matching glob_pattern and creates symlinks in tgt_dir.

    :param tgt_dir: target directory
    :param src_dir: source directory, must be absolute path 
    :param glob_pattern: the shell glob pattern (ls style) to match against files
    """
    files = glob.glob(osp.join(src_dir, glob_pattern))
    map(lambda f: symlink_unless_exists(f, osp.join(tgt_dir, osp.basename(f))), files)


def update_time_keys(time_utc, which, num_domains):
    """
    Returns a dict that can be use to set the start_xxxx or end_xxxx keys in time_control.

    The update is the same for all domains.

    :param time_utc: the UTC time to place into the dictionary
    :param which: must be either 'start' or 'end', controls which time is set
    :param num_domains: number of domains for which we should setup time_control
    """
    res = {}
    year, mon, day = time_utc.year, time_utc.month, time_utc.day
    hour, minute, sec = time_utc.hour, time_utc.minute, time_utc.second
    key_seq = map(lambda x: which + x, ['_year', '_month', '_day', '_hour', '_minute', '_second'])
    return {key: [value] * num_domains for (key, value) in zip(key_seq, [year, mon, day, hour, minute, sec])}


def update_time_control(start_utc, end_utc, num_domains):
    """
    Generate dictionary keys that can be directly update()d into time_control section of input namelist.

    This fill in start_yyyy, end_yyyy, and run_yyy keys in time_control section.

    :param start_utc: simulation start time
    :param end_utc: simulation end time
    :param num_domains: number of domains for which we should setup time_control
    """
    tc_dict = update_time_keys(start_utc, 'start', num_domains)
    tc_dict.update(update_time_keys(end_utc, 'end', num_domains))

    total_seconds = (end_utc - start_utc).total_seconds()
    days, rest = divmod(int(total_seconds), 86400)
    hours, rest = divmod(rest, 3600)
    minutes, seconds = divmod(rest, 60)

    # note, the multiplicity of run_xxxx parameters is 1 irrespective of domain count
    tc_dict['run_days'] = days
    tc_dict['run_hours'] = hours
    tc_dict['run_minutes'] = minutes
    tc_dict['run_seconds'] = seconds

    return tc_dict


def update_namelist(nml, with_keys):
    """
    Update namelist with keys from the nested dictionary with_keys.  Each key is overwritten.

    The nested dictionary is first keyed by section, then by key & value.

    :param nml: the namelist to update
    :param with_keys: the nested dictionary with update instructions
    """
    for section, section_dict in with_keys.iteritems():
        nml[section].update(section_dict)


def update_ignitions(ign_specs, max_dom):
    """
    Build a dictionary with which we can update the input namelist with passed ignition specifications.

    :param ign_specs: the ignition specifications
    :param max_dom: the maximum number of domains
    """
    def set_ignition_val(dom_id, v):
        dvals = [0] * max_dom
        dvals[dom_id-1] = v
        return dvals

    keys = [ "fire_ignition_start_lat", "fire_ignition_end_lat",
             "fire_ignition_start_lon", "fire_ignition_end_lon", 
             "fire_ignition_start_time", "fire_ignition_end_time", 
             "fire_ignition_radius", "fire_ignition_ros" ]

    nml_igns = { 'ifire' : [0] * max_dom, 'fire_num_ignitions' : [0] * max_dom }
    for dom_str, dom_igns in ign_specs.iteritems():
        # ensure fire model is switched on in every domain with ignitions
        dom_id = int(dom_str)
        nml_igns['ifire'][dom_id-1] = 2
        nml_igns['fire_num_ignitions'][dom_id-1] = len(dom_igns)
        nml_igns['fire_fuel_read'][dom_id-1] = -1 # real fuel data from WPS
        nml_igns['fire_fuel_cat'][dom_id-1] = 1 # arbitrary, won't be used


        # for each ignition 
        for ndx,ign in enumerate(dom_igns):
            start, dur = ign['start_delay_s'], ign['duration_s']
            lat, lon = ign['lat'], ign['long']
            vals = [ lat, lat, lon, lon, start, start+dur, 200, 1 ]
            kv = dict(zip([x + str(ndx+1) for x in keys], [set_ignition_val(dom_id, v) for v in vals]))
            nml_igns.update(kv)

    return { 'fire' : nml_igns }


def timespec_to_utc_hour(ts_str, from_time = None):
    """
    Converts relative time specifications into a UTC datetime which is rounded down to the nearest hour.

    ts_str can be an ESMF time, in which case this function returns esmf_to_utc(ts_str) rounded
    down, or it can be in the format 'T+YY' or 'T-YY' and then it's a relative specification of minutes
    from from_time.

    :param ts_str: the time specification string
    :param from_time: an optional origin time w.r.t. which the timespec is resolved
    """
    if ts_str[0] == 'T':
        if from_time is None:
            from_time = datetime.utcnow()
        min_shift = int(ts_str[1:])
        dt = from_time + timedelta(minutes = min_shift)
        return dt.replace(minute = 0, second = 0)
    else:
        dt = esmf_to_utc(ts_str)
        return dt.replace(minute = 0, second = 0)



#
#  Geospatial utilities
#

def great_circle_distance(lon1, lat1, lon2, lat2):
    """
    Computes the great circle distance between two points given as (lon1,lat1), (lon2,lat2)
    in kilometers.
    
        d = great_circle_distance(lon1, lat1, lon2, lat2)
    """
    rlat1, rlat2 = np.pi * lat1 / 180.0, np.pi * lat2 / 180.0
    rlon1, rlon2 = np.pi * lon1 / 180.0, np.pi * lon2 / 180.0
    
    a = math.sin(0.5*(rlat1 - rlat2))**2 + math.cos(rlat1)*math.cos(rlat2)*math.sin(0.5*(rlon1 - rlon2))**2
    c = 2 * math.atan2(a**0.5, (1-a)**0.5)
    return 6371.0 * c



def find_closest_grid_point(slon, slat, glon, glat):
    """
    Finds the closest grid point to the given station longitude/lattitude.
    """
    closest = np.argmin((slon - glon)**2 + (slat - glat)**2)
    return np.unravel_index(closest, glon.shape)


# FIXME: coordinate lookups and transforms go here

