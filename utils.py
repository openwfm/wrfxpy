from datetime import datetime
import pytz
import os
import os.path as osp
import glob


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
    key_seq = map(lambda x: which + x, [ '_year', '_month', '_day', '_hour', '_minute', '_second' ])
    return { key : [value] * num_domains for (key, value) in zip(key_seq, [year, mon, day, hour, minute, sec])}


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

