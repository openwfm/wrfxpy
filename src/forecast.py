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


from wrf.wrf_cloner import WRFCloner
from wrf.wrf_exec import Geogrid, Ungrib, Metgrid, Real, WRF
from wrf.wps_domains import WPSDomainLCC, WPSDomainConf

from utils import utc_to_esmf, symlink_matching_files, symlink_unless_exists, update_time_control, \
                  update_namelist, compute_fc_hours, esmf_to_utc, update_ignitions, make_dir, \
                  timespec_to_utc_hour
from vis.postprocessor import Postprocessor
from vis.var_wisdom import get_wisdom_variables

from ingest.grib_source import HRRR, NAM218, NARR
from fmda.fuel_moisture_da import assimilate_fm10_observations

import f90nml
from datetime import datetime, timedelta
import time, re, json, sys, logging
import os.path as osp
from multiprocessing import Process, Queue

import smtplib
from email.mime.text import MIMEText


class Dict(dict):
    """
    A dictionary that allows member access to its keys.
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


class JobState(Dict):
    """
    A coherent structure that holds information about the job.
    """

    def __init__(self, args):
        """
        Initialize the job state from the arguments dictionary.

        :param args: the forecast job arguments
        """
        super(JobState, self).__init__(args)
        self.fc_hrs = compute_fc_hours(self.start_utc, self.end_utc)
        self.grib_source = self.resolve_grib_source(self.grib_source)
        self.job_id = 'wfc-' + self.grid_code + '-' + utc_to_esmf(self.start_utc) + '-{0:02d}'.format(self.fc_hrs)
        self.emails = self.parse_emails(args)
        self.domains = args['domains']
        self.ignitions = args.get('ignitions', None)
        self.fmda = self.parse_fmda(args)
        self.postproc = args['postproc']
        self.wrfxpy_dir = args['sys_install_dir']

    
    def resolve_grib_source(self, gs_name):
        """
        Creates the right GribSource object from the name.
        
        :param gs_name: the name of the grib source
        """
        if gs_name == 'HRRR':
            return HRRR('ingest')
        elif gs_name == 'NAM':
            return NAM218('ingest')
        elif gs_name == 'NARR':
            return NARR('ingest')
        else:
            raise ValueError('Unrecognized grib_source %s' % gs_name)


    def parse_fmda(self, args):
        """
        Parse information inside the FMDA blob, if any.

        :param args: the forecast job argument dictionary
        """
        if 'fuel_moisture_da' in args:
            fmda = args['fuel_moisture_da']
            self.fmda = Dict({'token' : fmda['mesowest_token'], 'domains' : fmda['domains']})
        else:
            self.fmda = None


    def parse_emails(self, args):
        """
        Parse the definition of e-mail notifications

        :param args: the forecast job argument dictionary
        """
        if 'email_notifications' in args:
            emails = args['email_notifications']
            self.emails = Dict({'to' : emails['to'], 'events' : emails['events'],
                                'server' : emails.get('smtp_server', 'localhost'),
                                'origin' : emails.get('from', 'wrfxpy@gross.ucdenver.edu')})
        else:
            self.emails = None


def send_email(js, event, body):
    """
    Sends an e-mail with body <body> according to the e-mail parameters (constructed in execute) if the stated <event>
    is contained in the appropriate array.

    :param js: the JobState structure containing confiuration info
    :param event: name of the event firing the e-mail, the e-mail will not be sent unless <event> appears in the events array
    :param body: the body that will be placed into the e-mail
    """
    if js.emails is not None:
        if event in js.emails.events:
            mail_serv = smtplib.SMTP(js.emails.server)
            msg = MIMEText(body)
            msg['Subject'] = 'Job %s event %s notification' % (js.job_id, event)
            msg['From'] = js.emails.origin
            msg['To'] = js.emails.to
            mail_serv.sendmail(js.emails.origin, [js.emails.to], msg.as_string())
            mail_serv.quit()


def retrieve_gribs_and_run_ungrib(js, q):
    """
    This function retrieves required GRIB files and runs ungrib.

    It returns either 'SUCCESS' or 'FAILURE' on completion.

    :param js: the JobState object containing the forecast configuration
    :param q: the multiprocessing Queue into which we will send either 'SUCCESS' or 'FAILURE'
    """
    try:
        logging.info("retrieving GRIB files.")

        # step 3: retrieve required GRIB files from the grib_source, symlink into GRIBFILE.XYZ links into wps
        manifest = js.grib_source.retrieve_gribs(js.start_utc, js.end_utc)
        js.grib_source.symlink_gribs(manifest, js.wps_dir)

        send_email(js, 'grib2', 'Job %s - %d GRIB2 files downloaded.' % (js.job_id, len(manifest)))
        logging.info("running UNGRIB")

        # step 4: patch namelist for ungrib end execute ungrib
        f90nml.write(js.wps_nml, osp.join(js.wps_dir, 'namelist.wps'), force=True)

        Ungrib(js.wps_dir).execute().check_output()

        send_email(js, 'ungrib', 'Job %s - ungrib complete.' % js.job_id)
        logging.info('UNGRIB complete')
        q.put('SUCCESS')

    except Exception as e:
        logging.error('GRIB2/UNGRIB step failed with exception %s' % repr(e))
        q.put('FAILURE')


def run_geogrid(js, q):
    """
    This function runs geogrid or links in precomputed grid files as required.

    :param js: the JobState object containing the forecast configuration
    :param q: the multiprocessing Queue into which we will send either 'SUCCESS' or 'FAILURE'
    """
    try:
        logging.info("running GEOGRID")
        Geogrid(js.wps_dir).execute().check_output()
        logging.info('GEOGRID complete')
        
        send_email(js, 'geogrid', 'GEOGRID complete.')
        q.put('SUCCESS')

    except Exception as e:
        logging.error('GEOGRID step failed with exception %s' % repr(e))
        q.put('FAILURE')


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
        grib_source: a string identifying a valid GRIB2 source
        wps_namelist_path: the path to the namelist.wps file that will be used as template
        wrf_namelist_path: the path to the namelist.input file that will be used as template
        fire_namelist_path: the path to the namelist.fire file that will be used as template
        geogrid_path: the path to the geogrid data directory providing terrain/fuel data
        email_notification: dictionary containing keys address and events indicating when a mail should be fired off
    :return:
    """
    logging.basicConfig(level=logging.INFO)

    # initialize the job state from the arguments
    js = JobState(args)

    logging.info("job %s starting [%d hours to forecast]." % (js.job_id, js.fc_hrs))
    send_email(js, 'start', 'Job %s started.' % js.job_id)

    # read in all namelists
    js.wps_nml = f90nml.read(args['wps_namelist_path'])
    js.wrf_nml = f90nml.read(args['wrf_namelist_path'])
    js.fire_nml = f90nml.read(args['fire_namelist_path'])
    js.ems_nml = None
    if 'emissions_namelist_path' in args:
        js.ems_nml = f90nml.read(args['emissions_namelist_path'])
    
    # Parse and setup the domain configuration
    js.domain_conf = WPSDomainConf(js.domains)

    num_doms = len(js.domain_conf)
    js.wps_nml['share']['start_date'] = [utc_to_esmf(js.start_utc)] * num_doms
    js.wps_nml['share']['end_date'] = [utc_to_esmf(js.end_utc)] * num_doms
    js.wps_nml['share']['interval_seconds'] = 3600

    logging.info("number of domains defined is %d." % num_doms)

    # build directories in workspace
    js.wps_dir = osp.abspath(osp.join(js.workspace_dir, js.job_id, 'wps'))
    js.wrf_dir = osp.abspath(osp.join(js.workspace_dir, js.job_id, 'wrf'))

    logging.info("cloning WPS into %s" % js.wps_dir)

    # step 1: clone WPS and WRF directories
    cln = WRFCloner(args)
    cln.clone_wps(js.wps_dir, js.grib_source.vtables(), [])

    # step 2: process domain information and patch namelist for geogrid
    js.wps_nml['geogrid']['geog_data_path'] = args['geogrid_path']
    js.domain_conf.prepare_for_geogrid(js.wps_nml, js.wrf_nml, js.wrfxpy_dir, js.wps_dir)
    f90nml.write(js.wps_nml, osp.join(js.wps_dir, 'namelist.wps'), force=True)

    # do steps 2 & 3 & 4 in parallel (two execution streams)
    #  -> GEOGRID ->
    #  -> GRIB2 download ->  UNGRIB ->

    proc_q = Queue()
    geogrid_proc = Process(target=run_geogrid, args=(js, proc_q))
    grib_proc = Process(target=retrieve_gribs_and_run_ungrib, args=(js, proc_q))

    geogrid_proc.start()
    grib_proc.start()

    # wait until both tasks are done
    geogrid_proc.join()
    grib_proc.join()

    if proc_q.get() != 'SUCCESS':
        return

    if proc_q.get() != 'SUCCESS':
        return

    proc_q.close()

    # step 5: execute metgrid after ensuring all grids will be processed
    js.domain_conf.prepare_for_metgrid(js.wps_nml)
    f90nml.write(js.wps_nml, osp.join(js.wps_dir, 'namelist.wps'), force=True)

    logging.info("running METGRID")
    Metgrid(js.wps_dir).execute().check_output()

    send_email(js, 'metgrid', 'Job %s - metgrid complete.' % js.job_id)
    logging.info("cloning WRF into %s" % js.wrf_dir)

    # step 6: clone wrf directory, symlink all met_em* files
    cln.clone_wrf(js.wrf_dir, [])
    symlink_matching_files(js.wrf_dir, js.wps_dir, "met_em*")

    logging.info("running REAL")

    # step 7: patch input namelist, fire namelist, emissions namelist (if required)
    #         and execute real.exe
    time_ctrl = update_time_control(js.start_utc, js.end_utc, num_doms)
    js.wrf_nml['time_control'].update(time_ctrl)
    update_namelist(js.wrf_nml, js.grib_source.namelist_keys())
    if 'ignitions' in args:
        update_namelist(js.wrf_nml, update_ignitions(js.ignitions, num_doms))

    # if we have an emissions namelist, automatically turn on the tracers
    if js.ems_nml is not None:
        f90nml.write(js.ems_nml, osp.join(js.wrf_dir, 'namelist.fire_emissions'), force=True)
        js.wrf_nml['dynamics']['tracer_opt'] = [2] * num_doms

    f90nml.write(js.wrf_nml, osp.join(js.wrf_dir, 'namelist.input'), force=True)

    f90nml.write(js.fire_nml, osp.join(js.wrf_dir, 'namelist.fire'), force=True)
    
    Real(js.wrf_dir).execute().check_output()

    # step 8: if requested, do fuel moisture DA
    if js.fmda is not None:
        logging.info('running fuel moisture data assimilation')
        for dom in js.fmda.domains:
            assimilate_fm10_observations(osp.join(wrf_dir, 'wrfinput_d%02d' % dom), None, js.fmda.token)

    logging.info('submitting WRF job')
    send_email(js, 'wrf_submit', 'Job %s - wrf job submitted.' % js.job_id)

    # step 8: execute wrf.exe on parallel backend
    js.task_id = "sim-" + js.grid_code + "-" + utc_to_esmf(js.start_utc)[:10]
    WRF(js.wrf_dir, js.qman).submit(js.task_id, js.num_nodes, js.ppn, js.wall_time_hrs)

    send_email(js, 'wrf_exec', 'Job %s - wrf job starting now with id %s.' % (js.job_id, js.task_id))
    logging.info("WRF job submitted with id %s, waiting for rsl.error.0000" % js.task_id)

    # step 9: wait for appearance of rsl.error.0000 and open it
    wrf_out = None
    while wrf_out is None:
        try:
            wrf_out = open(osp.join(js.wrf_dir, 'rsl.error.0000'))
            break
        except IOError:
            logging.info('forecast: waiting 10 seconds for rsl.error.0000 file')
        
        time.sleep(5)
    
    logging.info('Detected rsl.error.0000')

    # step 10: track log output and check for history writes fro WRF
    pp_instr, pp = {}, None
    if js.postproc is not None:
        js.pp_dir = osp.join(js.workspace_dir, js.job_id, "products")
        make_dir(js.pp_dir)
        pp = Postprocessor(js.pp_dir, 'wfc-' + js.grid_code)

    while True:
        line = wrf_out.readline().strip()
        if not line:
            time.sleep(0.2)
            continue

        if "SUCCESS COMPLETE WRF" in line:
            send_email(js, 'complete', 'Job %s - wrf job complete SUCCESS.' % js.job_id)
            logging.info("WRF completion detected.")
            break

        if "Timing for Writing wrfout" in line:
            esmf_time,domain_str = re.match(r'.*wrfout_d.._([0-9_\-:]{19}) for domain\ +(\d+):' ,line).groups()
            dom_id = int(domain_str)
            logging.info("Detected history write for domain %d for time %s." % (dom_id, esmf_time))
            if str(dom_id) in js.postproc:
                var_list = [str(x) for x in js.postproc[str(dom_id)]]
                logging.info("Executing postproc instructions for vars %s for domain %d." % (str(var_list), dom_id))
                wrfout_path = osp.join(js.wrf_dir,"wrfout_d%02d_%s" % (dom_id, utc_to_esmf(js.start_utc))) 
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

    optional_files = [('emissions_namelist_path', 'Non-existent namelist template %s')]

    # check each path that should exist
    for key, err in required_files:
        if not osp.exists(args[key]):
            raise OSError(err % args[key])

    # check each path that should exist
    for key, err in optional_files:
        if key in args:
            if not osp.exists(args[key]):
                raise OSError(err % args[key])

    # check for valid grib source
    if args['grib_source'] not in ['HRRR', 'NAM', 'NARR']:
        raise ValueError('Invalid grib source, must be one of HRRR, NAM, NARR')

    # if precomputed key is present, check files linked in
    if 'precomputed' in args:
      for key,path in args['precomputed'].iteritems():
          if not osp.exists(path):
              raise OSError('Precomputed entry %s points to non-existent file %s' % (key,path))

    # check if the postprocessor knows how to handle all variables
    wvs = get_wisdom_variables()
    failing = False
    if 'postproc' in args:
        for dom in args['postproc'].keys():
            for vname in args['postproc'][dom]:
                if vname not in wvs:
                    logging.error('unrecognized variable %s in postproc key for domain %s.' % (vname, dom))
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
            'wrf_namelist_path': 'etc/nlists/colorado-3k.input.tracers',
            'fire_namelist_path': 'etc/nlists/colorado-3k.fire',
            'emissions_namelist_path' : 'etc/nlists/colorado-3k.fire_emissions',
            'precomputed' : { 'geo_em.d01.nc' : 'precomputed/colorado/geo_em.d01.nc',
                              'geo_em.d02.nc' : 'precomputed/colorado/geo_em.d02.nc' },
            'geogrid_path': '/share_home/mvejmelka/Packages/WPS-GEOG',
            'sys_install_dir': '/share_home/mvejmelka/Projects/wrfxpy-dev',
            'postproc' : {
                "2" : [ "T2", "PSFC", "WINDSPD", "WINDVEC", "FIRE_AREA", "FIRE_HFX", "SMOKE_INT", "F_ROS" ]
            },
            "ignitions" : {
                "2" : [ {
                    "start_delay_s" : 600,
                    "duration_s" : 240,
                    "lat" : 39.894264,
                    "long" : -103.903222
                    } ]
            },
            'num_nodes': 10,
            'ppn': 12,
            'wall_time_hrs': 3,
            'qman': 'sge',
            'start_utc': datetime(2016, 1, 17, 16, 0, 0),
            'end_utc': datetime(2016, 1, 17, 22, 0, 0)
           }

    verify_inputs(args)

    execute(args)


def process_arguments(args):
    """
    Convert arguments passed into program via the JSON configuration file.

    Transforms unicode strings into standard strings.

    :param args: the input arguments
    """
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

    # configure the basic logger
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    process_arguments(args)

    # sanity check
    verify_inputs(args)

    # execute the job
    execute(args)

