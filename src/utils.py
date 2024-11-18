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
import pprint
import logging
import dill
import json
import pickle
import inspect
import shutil
import psutil
import requests
import socket
import collections
from clamp2mesh import nearest_idx, get_subgrid_coordinates, clamp2mesh

try:
    import fire_init
    _fire_init_plugin = True
except:
    _fire_init_plugin = False

class Dict(dict):
    """
    A dictionary that allows member access to its keys.
    A convenience class.
    """

    def __init__(self, d):
        """
        Updates itself with d.
        """
        self.update(d)

    def __getattr__(self, item):
        return self[item]

    def __setattr__(self, item, value):
        self[item] = value

    def __getitem__(self, item):
        if item in self:
            return super().__getitem__(item)
        else:
            for key in self:
                if isinstance(key,(range,tuple)) and item in key:
                    return super().__getitem__(key)
            raise KeyError(item)

    def keys(self):
        if any([isinstance(key,(range,tuple)) for key in self]):
            keys = []
            for key in self:
                if isinstance(key,(range,tuple)):
                    for k in key:
                        keys.append(k)
                else:
                    keys.append(key)
            return keys
        else:
            return super().keys()

def save(obj, file):
    with open(file,'wb') as output:
        dill.dump(obj, output )

def load(file):
    with open(file,'rb') as input:
        returnitem = dill.load(input)
        return returnitem

def traceargs():
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    for i in args:
        print(("    %s:\n%s" % (i, pprint.pformat(values[i]))))

def dump(obj,title):
    frame = inspect.currentframe()
    outframe = inspect.getouterframes(frame, 2)
    logging.info(outframe[1][1] + ':' + outframe[1][3] + ":" + title + ':\n' + pprint.pformat(obj,width=-1))

def check_obj(obj, title):
    if pprint.isreadable(obj):
        logging.info(title + " is readable")
    else:
        logging.info(title + " is NOT readable")
    try:
        json.dumps(obj)
        logging.info(title + " is JSON serializable")
    except TypeError as err:
        logging.error(title + " is not JSON serializable")
        logging.error(err)
    try:
        s=pickle.dumps(obj,pickle.HIGHEST_PROTOCOL)
        pickle.loads(s)
        logging.info(title + " can be pickled")
    except:
        logging.error(title + " could not be picked and unpickled")

def kill_process(pid):
    if pid is not None:
        logging.info('Killing process %s and children' % pid)
        try:
            parent = psutil.Process(pid)
            try:
                for child in parent.children(recursive=True):
                    child.kill()
                parent.kill()
            except:
                logging.error('Could not get them all, check for runaways')
        except:
            logging.error('Process %s does not exist' % pid)

def process_create_time(pid):
    """
    Check if process exists and return its creation time.

    :param pid: process number
    :return: the process creation time, None if the process not exist, -1 if pid is None.
    """
    if pid is None:
        return -1
    else:
        try:
            return psutil.Process(pid).create_time()
        except:
            return None

def file_exists(path):
    """
    Ensure a file in a path exists and is readable.

    :param path: path to an existent file
    :return: boolean if the file exists and is readable
    """
    return (os.path.exists(path) and os.access(path,os.R_OK))

def ensure_dir(path):
    """
    Ensure all directories in path if a file exist, for convenience return path itself.

    :param path: the path whose directories should exist
    :return: the path back for convenience
    """
    path_dir = osp.dirname(path)
    if not osp.exists(path_dir):
        os.makedirs(path_dir)
    return path

def delete(dir):
    if osp.exists(dir):
        try:
            shutil.rmtree(dir)
        except Exception as e:
            # This is often caused by hidden files on NSF mounted volumes, which is harmless.
            logging.warning(str(e))

def make_clean_dir(dir):
    """
    Create a clean directory; delete first if it exists.
    """
    logging.info('Deleting existing directory %s to make a clean one' % dir)
    delete(dir)
    if not osp.exists(dir):
        os.makedirs(dir)

def make_dir(dir):
    """
    Create a directory if it does not exist. Creates any intermediate dirs as necessary.
    For convenience return the director path itself

    :param dir: the directory to be created
    :retun: same as input
    """
    if not osp.exists(dir):
        os.makedirs(dir)
    return dir


def symlink_unless_exists(link_tgt, link_loc):
    """
    Create a symlink at link_loc pointing to link_tgt unless file already exists.

    :param link_tgt: link target
    :param link_loc: link location
    """
    logging.info('Linking %s -> %s' % (link_loc, link_tgt))
    if osp.isfile(link_tgt) or osp.isdir(link_tgt):
        if not osp.lexists(link_loc):
            os.symlink(link_tgt, link_loc)
        else:
            logging.warning('Link %s already exists' % link_loc)
    else:
        logging.error('Link target %s does not exist' % link_tgt)

def remove(tgt):
    """
    os.remove wrapper
    """
    if osp.isfile(tgt):
        logging.info('remove: file %s exists, removing' % tgt)
        os.remove(tgt)

def force_copy(src,tgt):
    """
    remove target if exists and copy there, making directories as needed
    """
    remove(tgt)
    shutil.copy(src,ensure_dir(tgt))

def append2file(addl,base):
    """
    Append a file to another
    """
    logging.info("appending file %s to %s" % (addl,base)) 
    with open(base,'a') as base_file:
        with open(addl,'r') as addl_file:
            base_file.write(addl_file.read())

def link2copy(src):
    """
    replace link by a copy of the target file
    """
    try:
        link_target = os.readlink(src)
    except OSError as e:
        return
    logging.info("replacing soft link %s -> %s by a copy" % (src,link_target))
    os.remove(src)
    shutil.copy(link_target,src)

def move(src,tgt):
    """
    shutil.move wrapper
    """
    logging.info('moving %s to %s' % (src, tgt))
    remove(tgt)
    shutil.move(src,tgt)

def cache_file(path, cache_dir):
    """
    Move file at the path to cache directory and replace by symlink
    except when it is symlink to the cache direcory already
    except if the file is symlink to elsewhere, copy the file to cache directory and replace by symlink

    :path: file name
    :param cache_dir: source directory, must be absolute path
    """
    if not osp.isdir(cache_dir):
        logging.error('%s is not directory' % str(cache_dir))
        raise Exception('%s is not directory' % str(cache_dir))
    if not osp.isfile(path):
        logging.error('%s is not file' % str(path))
        raise Exception('%s is not directory' % str(cache_dir))
    dst = osp.join(cache_dir,osp.basename(path))
    if osp.islink(path):
        if osp.dirname(osp.realpath(path)) is osp.realpath(cache_dir):
            logging.debug('%s is link to %s already' % (path, cache_dir))
            if osp.basename(osp.realpath(path)) is not osp.basename(path):
                logging.error('link %s -> %s does not match' % (path,osp.realpath(path)))
                raise Exception('link %s -> %s does not match' % (path,osp.realpath(path)))
        else:
            src = osp.realpath(path)
            logging.info('Copying %s to %s' % (src, dst))
            shutil.copyfile(src,dst)
    else:
        logging.info('Moving %s to %s' % (path, dst))
        shutil.move(path,dst)
    symlink_unless_exists(dst, path)

def esmf_to_utc(esmf):
    """
    Converts an ESMF datetime into a UTC datetime.

    :param esmf: date & time from YYYY-MM-DD_hh:mm:ss (ESMF) format
    :return: a python datetime with UTC timezone
    """
    if esmf is None:
        return None
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

def utc_to_utcf(utc):
    """
    Converts a UTC datetime into UTC string format %Y-%m-%dT%H:%M:%SZ.

    :param utc: python UTC datetime
    :return: a string in UTC string format %Y-%m-%dT%H:%M:%SZ
    """
    return "%04d-%02d-%02dT%02d:%02d:%02dZ" % (utc.year, utc.month, utc.day, utc.hour, utc.minute, utc.second)


def round_time_to_hour(utc, up=False, period_hours=1):
    """
    Round the utc time to the nearest hour up or down.

    :param utc: the datetime
    :param up: round up
    :param hours: round to multiple of this from the start of day
    :return: a new datetime rounded as specified
    """
    if utc is None:
        return None
    tm = utc + timedelta(hours=1, seconds=-1) if up else utc
    tm = tm.replace(minute=0, second=0)
    h = period_hours * ((tm.hour + period_hours - 1 if up else tm.hour) // period_hours)
    tm = tm + timedelta(hours=h-tm.hour) if h > tm.hour else tm - timedelta(hours=tm.hour - h)
    return tm


def timedelta_hours(timedelta_utc, up = True):
    """
    Compute difference of times in hours, rounds up in case of incomplete hours.

    :param delta_utc: difference of two datetime objects
    :param up: round up if True, down in False
    :return: number of hours rounded to whole hour
    """
    return int(timedelta_utc.total_seconds() + ((3600 - 0.001)if up else 0 )) // 3600


def symlink_matching_files(tgt_dir, src_dir, glob_pattern):
    """
    Retrieves all files in directory src_dir matching glob_pattern and creates symlinks in tgt_dir.

    :param tgt_dir: target directory
    :param src_dir: source directory, must be absolute path
    :param glob_pattern: the shell glob pattern (ls style) to match against files
    """
    files = glob.glob(osp.join(src_dir, glob_pattern))
    list(map(lambda f: symlink_unless_exists(f, osp.join(tgt_dir, osp.basename(f))), files))


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
    key_seq = [which + x for x in ['_year', '_month', '_day', '_hour', '_minute', '_second']]
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
    for section, section_dict in with_keys.items():
        nml[section].update(section_dict)


def render_ignitions(js, max_dom):
    """
    Build a dictionary with which we can update the input namelist with passed ignition specifications.

    :param js: the job state
    :param max_dom: the maximum number of domains
    """
    duration_default = 60.
    radius_default = 60.
    ros_default = 1.
    
    def set_ignition_val(dom_id, v):
        dvals = [0] * max_dom
        dvals[dom_id-1] = v
        return dvals

    # extract the original ignition specs and the original start time
    dom_specs = js.domains
    ign_specs = js.ignitions
    orig_start_time = js.orig_start_utc

    keys = [ "fire_ignition_start_lat", "fire_ignition_end_lat",
             "fire_ignition_start_lon", "fire_ignition_end_lon",
             "fire_ignition_start_time", "fire_ignition_end_time",
             "fire_ignition_radius", "fire_ignition_ros" ]

    nml_fire = { 'ifire' : [0] * max_dom, 'fire_num_ignitions' : [0] * max_dom,
                 'fire_fuel_read' : [0] * max_dom, 'fire_fuel_cat' : [1] * max_dom,
                 'fmoist_run' : [False] * max_dom, 'fmoist_interp' : [False] * max_dom,
                 'fire_fmc_read' : [0] * max_dom, 'fmoist_dt' : [600] * max_dom,
                 'fire_viscosity' : [0] * max_dom, 'fire_wind_log_interp': [0] * max_dom,
                 'fire_use_windrf': [0] * max_dom }
    
    if js.use_realtime:
        fire_perimeter_time = js.get('fire_perimeter_time', 7200.)
        nml_fire.update({'fire_perimeter_time': [0] * max_dom})
        nml_fire['fire_perimeter_time'][max_dom-1] = fire_perimeter_time
        nml_fire.update({'fire_boundary_guard': [0] * max_dom})
        nml_fire['fire_boundary_guard'][max_dom-1] = -1
        ign_specs = {str(max_dom): []}
    elif js.use_tign_ignition:
        ign_times = [
            (
                timespec_to_utc(ign['time_utc'], orig_start_time),
                ign.get('duration_s',duration_default)
            ) 
            for dom_igns in ign_specs.values() for ign in dom_igns
        ]
        if len(ign_times):
            fire_tign_in_time = max([int((ign_time - js.start_utc).total_seconds())+ign_dur for ign_time,ign_dur in ign_times])
            nml_fire.update({'fire_tign_in_time': [0] * max_dom})
            nml_fire['fire_tign_in_time'][max_dom-1] = fire_tign_in_time
        ign_specs = {str(max_dom): []}
    
    for dom_str, dom_igns in ign_specs.items():
        dom_id = int(dom_str)
        # ensure fire model is switched on in every domain with ignitions
        nml_fire['ifire'][dom_id-1] = js.get('ifire',1)
        nml_fire['fire_wind_log_interp'][dom_id-1] = 1
        nml_fire['fire_use_windrf'][dom_id-1] = 2
        if not (js.use_tign_ignition or js.use_realtime):
            nml_fire['fire_num_ignitions'][dom_id-1] = len(dom_igns)
        nml_fire['fire_fuel_read'][dom_id-1] = -1 # real fuel data from WPS
        nml_fire['fire_fuel_cat'][dom_id-1] = 1 # arbitrary, won't be used
        nml_fire['fmoist_run'][dom_id-1] = True # use the fuel moisture model
        nml_fire['fmoist_interp'][dom_id-1] = True # interpolate fm onto fire mesh
        nml_fire['fire_fmc_read'][dom_id-1] = 0 # use wrfinput and/or running moisture model

        # for each ignition
        for ndx,ign in enumerate(dom_igns):
            start_time = timespec_to_utc(ign['time_utc'], orig_start_time)
            start = int((start_time - js.start_utc).total_seconds())
            dur = ign.get('duration_s', duration_default)
            radius = ign.get('radius', radius_default)
            ros = ign.get('ros', ros_default)
            if 'latlon' in ign.keys():
                lat,lon = ign['latlon']
                nlon,nlat = clamp2mesh(glob.glob(osp.join(js.wps_dir,'geo_em.d{:02d}*'.format(max_dom)))[0], float(lon), float(lat))
                vals = [ nlat, nlat, nlon, nlon, start, start+dur, radius, ros]
            elif 'start_latlon' in ign.keys() and 'end_latlon' in ign.keys():
                lat,lon = ign['start_latlon']
                slon,slat = clamp2mesh(glob.glob(osp.join(js.wps_dir,'geo_em.d{:02d}*'.format(max_dom)))[0], float(lon), float(lat))
                lat,lon = ign['end_latlon']
                elon,elat = clamp2mesh(glob.glob(osp.join(js.wps_dir,'geo_em.d{:02d}*'.format(max_dom)))[0], float(lon), float(lat))
                vals = [ slat, elat, slon, elon, start, start+dur, radius, ros]
            kv = dict(list(zip([x + str(ndx+1) for x in keys], [set_ignition_val(dom_id, v) for v in vals])))
            nml_fire.update(kv)

    # activate fire at domains with subgrid ratio specified
    for dom_str, dom_spec in dom_specs.items():
        dom_id = int(dom_str)
        if any([sr > 0 for sr in dom_spec['subgrid_ratio']]):
            nml_fire['ifire'][dom_id-1] = 1 # switch on fire model
            nml_fire['fmoist_run'][dom_id-1] = True # use the fuel moisture model
            nml_fire['fmoist_interp'][dom_id-1] = True # interpolate fm onto fire mesh
            nml_fire['fire_fmc_read'][dom_id-1] = 0 # use wrfinput and/or running moisture model

    return { 'fire' : nml_fire }

def process_ignitions(js):
    """
    Use fire_init package to process ignitions from user defined fire ignition points and lines.

    :param js: the job state
    """
    if not _fire_init_plugin: 
        logging.error('fire_init plugin not installed')
        raise Exception('fire_init plugin is required for selected ignitions')
    duration_default = 60.
    radius_default = 60.
    ros_default = 1.
    orig_start_time = js.orig_start_utc
    for dom,igns in js.ignitions.items():
        wrf_path = osp.join(js.wrf_dir, 'wrfinput_d%02d' % int(dom))
        fxlon,fxlat = get_subgrid_coordinates(wrf_path, strip=False)
        fire_data = {'points': [], 'lines': []}
        points = []
        lines = []
        for ign in igns:
            if 'line_id' in ign.keys():
                lines.append(ign)
            else:
                points.append(ign)
        line_ids = [ign['line_id'] for ign in lines]
        if len(line_ids):
            for line_id in np.unique(line_ids):
                line_points = [ign for ign in lines if ign['line_id'] == line_id]
                line_data = {
                    'time': [], 'ix': [], 'iy': [],
                    'lon': [], 'lat': [],  'duration': [], 
                    'radius': [], 'ros': []
                }
                for line_point in line_points:
                    time_utc = timespec_to_utc(line_point['time_utc'], orig_start_time)
                    time = int((time_utc - js.start_utc).total_seconds())
                    lat,lon = line_point['latlon']
                    iy,ix = nearest_idx(fxlon, fxlat, float(lon), float(lat))
                    lon = fxlon[iy, ix]
                    lat = fxlat[iy, ix]
                    duration = line_point.get('duration_s', duration_default)
                    radius = line_point.get('radius', radius_default)
                    ros = line_point.get('ros', ros_default)
                    line_data['time'].append(time)
                    line_data['ix'].append(ix)
                    line_data['iy'].append(iy)
                    line_data['lon'].append(lon)
                    line_data['lat'].append(lat)
                    line_data['duration'].append(duration)
                    line_data['radius'].append(radius)
                    line_data['ros'].append(ros)
                fire_data['lines'].append(line_data)
        if len(points):
            for point in points:
                time_utc = timespec_to_utc(point['time_utc'], orig_start_time)
                time = int((time_utc - js.start_utc).total_seconds())
                lat,lon = point['latlon']
                ix,iy = nearest_idx(fxlon, fxlat, float(lon), float(lat))
                duration = point.get('duration_s', duration_default)
                radius = point.get('radius', radius_default)
                ros = point.get('ros', ros_default)
                point_data = {
                    'time': time, 'ix': ix, 'iy': iy, 
                    'lon': lon, 'lat': lat, 'duration': duration, 
                    'radius': radius, 'ros': ros
                }
                fire_data['points'].append(point_data)
        TIGN_G = fire_init.process_ignitions_to_tign(wrf_path, fire_data)   
    if js.burn_plot_boundary is not None:
        coords = [c[::-1] for c in js.burn_plot_boundary]
        if coords[0] != coords[-1]:
            coords.append(coords[0])
        wrf_path = osp.join(js.wrf_dir, 'wrfinput_d%02d' % js.max_dom)
        FUEL_MASK = fire_init.process_burn_plot_boundary(wrf_path, [[coords]])  
    fire_init.integrate_init(wrf_path, TIGN_G, FUEL_MASK, no_fuel_cat=js.fire_nml['fuel_scalars']['no_fuel_cat'])

def timespec_to_utc(ts_str, from_time = None):
    """
    Converts relative time specifications into a UTC datetime.

    ts_str can be an ESMF time, in which case this function returns esmf_to_utc(ts_str) rounded
    down, or it can be in the format 'T+YY' or 'T-YY' and then it's a relative specification of minutes
    from from_time.

    :param ts_str: the time specification string
    :param from_time: an optional origin time w.r.t. which the timespec is resolved
    """
    if ts_str is None:
        return None
    if ts_str[0] == 'T':
        if from_time is None:
            from_time = datetime.utcnow().replace(tzinfo=pytz.UTC)
        min_shift = int(ts_str[1:])
        return from_time + timedelta(minutes = min_shift)
    else:
        return esmf_to_utc(ts_str)

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


def load_sys_cfg():
    # load the system configuration
    sys_cfg = None
    try:
        sys_cfg = Dict(json.load(open('etc/conf.json')))
    except IOError:
        import sys
        logging.critical('Cannot find system configuration, have you created etc/conf.json?')
        sys.exit(2)
    
    # set defaults
    sys_cfg.sys_install_path = sys_cfg.get('sys_install_path',os.getcwd())
    # configuration defaults + make directories if they do not exist
    sys_cfg.workspace_path = make_dir(osp.abspath(sys_cfg.get('workspace_path','wksp')))
    sys_cfg.ingest_path = make_dir(osp.abspath(sys_cfg.get('ingest_path','ingest')))
    sys_cfg.cache_path = make_dir(osp.abspath(sys_cfg.get('cache_path','cache')))
    sys_cfg.ref_utc = sys_cfg.get('ref_utc',None)
    sys_cfg.ungrib_only = sys_cfg.get('ungrib_only',False)
    sys_cfg.use_wgrib2 = sys_cfg.get('use_wgrib2',False)
    sys_cfg.sat_only = sys_cfg.get('sat_only',False)
    sys_cfg.use_realtime = sys_cfg.get('use_realtime',False)
    sys_cfg.use_tign_ignition = sys_cfg.get('use_tign_ignition',False)
    return sys_cfg

class response_object(object):
    status_code = 0
    def __init__(self,status_code):
        self.status_code = status_code

def readhead(url,msg_level=1,retries=10):
    if msg_level > 0:
        logging.info('reading http head of %s ' % url)
    for n in range(retries):
        try:
            ret=requests.head(url)
            ret.raise_for_status()
        except Exception as e:
            if msg_level > 0:
                logging.error(e)
            ret = response_object(-1)
        if ret.status_code == 200:
            break
    return ret

def json2xml(d):
    if isinstance(d,dict):
        s = ''
        for key in d:
            s += '<%s>%s</%s>' % (key, json2xml(d[key]),key)
        return s
    else:
        return str(d)

def checkip():
    try:
        r = requests.get('http://checkip.dyndns.org')
    except Exception as e:
        logging.error(e)
        logging.error('Cannot get the public IP address, are you connected to the internet?')
        return None
    if r.status_code != requests.codes.ok:
        logging.error('checkip: Bad request')
        return None
    s = 'Current IP Address: '
    i = r.text.find(s) + len(s)
    j = r.text.find('</body>')
    ip = r.text[i:j]
    logging.info('Your public IP address is %s' % ip)
    return ip

def get_ip_address():
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    s.connect(("8.8.8.8", 80))
    return s.getsockname()[0]

def json_join(path,json_list):
    """
    Join local jsons in a singular json and remove the previous jsons

    :param path: local path to the jsons
    :param json_list: list of json names to join
    """
    manifest = Dict({})
    for jj in json_list:
        json_path = osp.join(path,str(jj)+'.json')
        try:
            f = json.load(open(json_path))
            manifest.update({jj: f})
        except:
            logging.warning('no satellite data for source %s in manifest json file %s' % (jj,json_path))
            manifest.update({jj: {}})
            pass
        remove(json_path)
    return manifest

def duplicates(replist):
    """
    Give dictionary of repeated elements (keys) and their indexes (array in values)

    :param replist: list to look for repetitions
    """
    counter=collections.Counter(replist)
    dups=[i for i in counter if counter[i]!=1]
    result={}
    for item in dups:
        result[item]=[i for i,j in enumerate(replist) if j==item]
    return result

def number_minutes(t_int,t_fin,dt):
    """
    Total number of minutes between two datetimes using a specific time step in minutes

    :param t_int: initial datetime
    :param t_fin: final datetime
    :param dt: time step in minutes
    """
    return int(np.floor((t_fin-t_int).total_seconds()/60./int(dt)))

def serial_json(obj):
    """
    JSON serializer for objects not serializable by default json code

    :param obj: object
    """
    if isinstance(obj, datetime):
        return utc_to_esmf(obj)
    raise TypeError("Type %s not serializable" % type(obj))

def inq(x):
    """
    Inquire numpy array for shape min max
    :param x: object
    :return string:
    """
    try:
        s= 'size %s min %g max %g' % (str(x.shape), np.min(x), np.max(x))
    except:
        s=str(x)
    return s

def addquotes(s):
    """
    add quotes to string
    """
    return '"'+s+'"'

def split_path(path):
    """
    split path into a list of components
    """
    path_list = []
    path = osp.dirname(path)
    while osp.basename(path):
        path_list.append(osp.basename(path))
        path = osp.dirname(path)
    path_list.reverse()
    return path_list
