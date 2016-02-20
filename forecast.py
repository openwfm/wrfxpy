#!/usr/bin/env python
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


from wrf_cloner import WRFCloner
from wrf_exec import Geogrid, Ungrib, Metgrid, Real, WRF
from utils import utc_to_esmf, symlink_matching_files, symlink_unless_exists, update_time_control, \
                  update_namelist, compute_fc_hours, esmf_to_utc, update_ignitions, make_dir, \
                  timespec_to_utc_hour
from postproc import Postprocessor
from grib_source import HRRR
from var_wisdom import get_wisdom_variables

import f90nml
from datetime import datetime, timedelta
import time, re, json, sys, logging
import os.path as osp


def make_job_id(grid_code, start_utc, fc_hrs):
    """
    Computes the job id from grid code, UTC start time and number of forecast hours.
    """
    return 'wfc-' + grid_code + '-' + utc_to_esmf(start_utc) + '-{0:02d}'.format(fc_hrs)


def execute(args):
    """
    Executes a weather/fire simulation.

    The args dictionary contains

    :param args: a dictionary with the following keys
        grid_code: the (unique) code of the grid that is used
        sys_install_dir: system installation directory
        start_utc: start time of simulation in UTC
        end_utc: end time of simulation in UTC
        workspace_dir: workspace directory
        wps_install_dir: installation directory of WPS that will be used
        wrf_install_dir: installation directory of WRF that will be used
        grib_source: a GribSource object that will be used to obtain GRIB files
        wps_namelist_path: the path to the namelist.wps file that will be used as template
        wrf_namelist_path: the path to the namelist.input file that will be used as template
        fire_namelist_path: the path to the namelist.fire file that will be used as template
        geogrid_path: the path to the geogrid data directory providing terrain/fuel data
    :return:
    """
    logging.basicConfig(level=logging.INFO)
    wksp_dir, grib_source = args['workspace_dir'], args['grib_source']

    # compute the job id
    grid_code, start_utc, end_utc = args['grid_code'], args['start_utc'], args['end_utc']
    fc_hrs = compute_fc_hours(start_utc, end_utc)
    job_id = make_job_id(grid_code, start_utc, fc_hrs)

    logging.info("job %s starting [%d hours to forecast]." % (job_id, fc_hrs))

    # read in all namelists
    wps_nml = f90nml.read(args['wps_namelist_path'])
    wrf_nml = f90nml.read(args['wrf_namelist_path'])
    fire_nml = f90nml.read(args['fire_namelist_path'])

    num_doms = int(wps_nml['share']['max_dom'])

    logging.info("number of domains is %d." % num_doms)

    # build directories in workspace
    wps_dir = osp.abspath(osp.join(wksp_dir, job_id, 'wps'))
    wrf_dir = osp.abspath(osp.join(wksp_dir, job_id, 'wrf'))

    logging.info("cloning WPS into %s" % wps_dir)

    # step 1: clone WPS and WRF directories
    cln = WRFCloner(args)
    cln.clone_wps(wps_dir, grib_source.vtables(), [])

    # step 2: patch namelist for geogrid and execute geogrid
    wps_nml['geogrid']['geog_data_path'] = args['geogrid_path']
    f90nml.write(wps_nml, osp.join(wps_dir, 'namelist.wps'), force=True)
    if 'precomputed' in args:
        logging.info('precomputed grids found, linking in ...')
        for grid, path in args['precomputed'].iteritems():
           symlink_unless_exists(osp.abspath(path), osp.join(wps_dir, grid))
    else:
        logging.info("running GEOGRID")
        Geogrid(wps_dir).execute().check_output()

    logging.info("retrieving GRIB files.")

    # step 3: retrieve required GRIB files from the grib_source, symlink into GRIBFILE.XYZ links into wps
    manifest = grib_source.retrieve_gribs(start_utc, end_utc)
    grib_source.symlink_gribs(manifest, wps_dir)

    logging.info("running UNGRIB")

    # step 4: patch namelist for ungrib end execute ungrib
    wps_nml['share']['start_date'] = [utc_to_esmf(start_utc)] * num_doms
    wps_nml['share']['end_date'] = [utc_to_esmf(end_utc)] * num_doms
    wps_nml['share']['interval_seconds'] = 3600
    f90nml.write(wps_nml, osp.join(wps_dir, 'namelist.wps'), force=True)

    Ungrib(wps_dir).execute().check_output()

    logging.info("running METGRID")

    # step 5: execute metgrid
    Metgrid(wps_dir).execute().check_output()

    logging.info("cloning WRF into %s" % wrf_dir)

    # step 6: clone wrf directory, symlink all met_em* files
    cln.clone_wrf(wrf_dir, [])
    symlink_matching_files(wrf_dir, wps_dir, "met_em*")

    logging.info("running REAL")

    # step 7: patch input namelist, fire namelist and execute real.exe
    time_ctrl = update_time_control(start_utc, end_utc, num_doms)
    wrf_nml['time_control'].update(time_ctrl)
    update_namelist(wrf_nml, grib_source.namelist_keys())
    f90nml.write(fire_nml, osp.join(wrf_dir, 'namelist.fire'), force=True)

    # render ignition specification into the wrf namelist
    if 'ignitions' in args:
        update_namelist(wrf_nml, update_ignitions(args['ignitions'], num_doms))
   
    f90nml.write(wrf_nml, osp.join(wrf_dir, 'namelist.input'), force=True)
    Real(wrf_dir).execute().check_output()

    logging.info("submitting WRF job")

    # step 8: execute wrf.exe on parallel backend
    qman, nnodes, ppn, wall_time_hrs = args['qman'], args['num_nodes'], args['ppn'], args['wall_time_hrs']
    task_id = "sim-" + grid_code + "-" + utc_to_esmf(start_utc)[:10]
    WRF(wrf_dir, qman).submit(task_id, nnodes, ppn, wall_time_hrs)

    logging.info("WRF job submitted with id %s, waiting for rsl.error.0000" % task_id)

    # step 9: wait for appearance of rsl.error.0000 and open it
    wrf_out = None
    while wrf_out is None:
        try:
            wrf_out = open(osp.join(wrf_dir, 'rsl.error.0000'))
            break
        except IOError:
            logging.info('forecast: waiting 10 seconds for rsl.error.0000 file')
        
        time.sleep(5)
    
    logging.info('Detected rsl.error.0000')

    # step 10: track log output and check for history writes fro WRF
    pp_instr, pp = {}, None
    if 'postproc' in args:
        pp_instr = args['postproc']
        pp_dir = osp.join(wksp_dir, job_id, "products")
        make_dir(pp_dir)
        pp = Postprocessor(pp_dir, 'wfc-' + grid_code )

    while True:
        line = wrf_out.readline().strip()
        if not line:
            time.sleep(0.2)
            continue

        if "SUCCESS COMPLETE WRF" in line:
            logging.info("WRF completion detected.")
            break

        if "Timing for Writing wrfout" in line:
            esmf_time,domain_str = re.match(r'.*wrfout_d.._([0-9_\-:]{19}) for domain\ +(\d+):' ,line).groups()
            dom_id = int(domain_str)
            logging.info("Detected history write in for domain %d for time %s." % (dom_id, esmf_time))
            if str(dom_id) in pp_instr:
                var_list = [str(x) for x in pp_instr[str(dom_id)]]
                logging.info("Executing postproc instructions for vars %s for domain %d." % (str(var_list), dom_id))
                wrfout_path = osp.join(wrf_dir,"wrfout_d%02d_%s" % (dom_id, utc_to_esmf(start_utc))) 
                pp.vars2kmz(wrfout_path, dom_id, esmf_time, var_list)
                pp.vars2png(wrfout_path, dom_id, esmf_time, var_list)


def verify_inputs(args):
    """
    Check if arguments (eventually) supplied to execute(...) are valid - if not exception is thrown.

    Arguments:
      args -- dictionary of arguments
    """
    # we don't check if job_id is a valid path
    required_files = [('sys_install_dir', 'Non-existent system installation directory %s'),
                      ('workspace_dir', 'Non-existent workspace directory %s'),
                      ('wps_install_dir', 'Non-existent WPS installation directory %s'),
                      ('wrf_install_dir', 'Non-existent WRF installation directory %s'),
                      ('wps_namelist_path', 'Non-existent WPS namelist template %s'),
                      ('wrf_namelist_path', 'Non-existent WRF namelist template %s'),
                      ('fire_namelist_path', 'Non-existent fire namelist template %s'),
                      ('geogrid_path', 'Non-existent geogrid data (WPS-GEOG) path %s')]

    # check each path that should exist
    for key, err in required_files:
        if not osp.exists(args[key]):
            raise OSError(err % args[key])

    # check if the postprocessor knows how to handle all variables
    wvs = get_wisdom_variables()
    failing = False
    if 'postproc' in args:
        for dom in args['postproc'].keys():
            for vname in args['postproc'][dom]:
                if vname not in wvs:
                    logger.error('unrecognized variable %s in postproc key for domain %s.' % (vname, dom))
                    failing = True
    if failing:
        raise ValueError('One or more unrecognized variables in postproc.')


def test():

    args = {'grid_code': 'colo2dv1',
            'workspace_dir': 'wksp',
            'wps_install_dir': '/share_home/mvejmelka/Packages/wrf-fire.openwfm.clamping2/WPS',
            'wrf_install_dir': '/share_home/mvejmelka/Packages/wrf-fire.openwfm.clamping2/WRFV3',
            'grib_source': HRRR('ingest'),
            'wps_namelist_path': 'etc/nlists/colorado-3k.wps',
            'wrf_namelist_path': 'etc/nlists/colorado-3k.input',
            'fire_namelist_path': 'etc/nlists/colorado-3k.fire',
            'precomputed' : { 'geo_em.d01.nc' : 'precomputed/colorado/geo_em.d01',
                              'geo_em.d02.nc' : 'precomputed/colorado/geo_em.d02' },
            'geogrid_path': '/share_home/mvejmelka/Packages/WPS-GEOG',
            'sys_install_dir': '/share_home/mvejmelka/Projects/wrfxpy',
            'num_nodes': 10,
            'ppn': 12,
            'wall_time_hrs': 3,
            'qman': 'sge',
            'start_utc': datetime(2016, 1, 17, 16, 0, 0),
            'end_utc': datetime(2016, 1, 17, 22, 0, 0)}

    verify_inputs(args)

    execute(args)


def process_arguments(args):
    """
    Convert arguments passed into program via the JSON configuration file.

    Transforms unicode strings into standard strings.

    :param args: the input arguments
    """

    # process the arguments
    args['grib_source'] = HRRR('ingest')

    # resolve possible relative time specifications
    args['start_utc'] = timespec_to_utc_hour(args['start_utc'])
    args['end_utc'] = timespec_to_utc_hour(args['end_utc'], args['start_utc'])

    for k, v in args.iteritems():
        if type(v) == unicode:
            args[k] = v.encode('ascii')


if __name__ == '__main__':

    # load configuration JSON
    cfg_str = open(sys.argv[1]).read()
    args = json.loads(cfg_str, 'ascii')

    # HRRR source is compulsory
    assert args['grib_source'] == 'HRRR'

    process_arguments(args)

    # sanity check
    verify_inputs(args)

    # execute the job
    execute(args)

