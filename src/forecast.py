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
                  update_namelist, compute_fc_hours, esmf_to_utc, render_ignitions, make_dir, \
                  timespec_to_utc, round_time_to_hour, Dict, dump, save, load, check_obj
from vis.postprocessor import Postprocessor
from vis.var_wisdom import get_wisdom_variables

from ingest.grib_source import HRRR, NAM218, NAM227, NARR
from fmda.fuel_moisture_da import assimilate_fm10_observations

from ssh_shuttle import send_product_to_server

import f90nml
from datetime import datetime, timedelta
import time, re, json, sys, logging
import os.path as osp
from multiprocessing import Process, Queue
import glob

import smtplib
from email.mime.text import MIMEText

import traceback
import pprint



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
        self.grib_source = self.resolve_grib_source(self.grib_source)
        self.start_utc = round_time_to_hour(self.start_utc, up=False, period_hours=self.grib_source.period_hours);
        self.end_utc = round_time_to_hour(self.end_utc, up=True, period_hours=self.grib_source.period_hours);
        self.fc_hrs = compute_fc_hours(self.start_utc, self.end_utc)
        self.job_id = 'wfc-' + self.grid_code + '-' + utc_to_esmf(self.start_utc) + '-{0:02d}'.format(self.fc_hrs)
        self.emails = self.parse_emails(args)
        self.domains = args['domains']
        self.ignitions = args.get('ignitions', None)
        self.fmda = self.parse_fmda(args)
        self.postproc = args['postproc']
        self.wrfxpy_dir = args['sys_install_path']
        self.args = args

    def resolve_grib_source(self, gs_name):
        """
        Creates the right GribSource object from the name.
        
        :param gs_name: the name of the grib source
        """
        if gs_name == 'HRRR':
            return HRRR('ingest')
        elif gs_name == 'NAM':
            return NAM218('ingest')
        elif gs_name == 'NAM227':
            return NAM227('ingest')
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
        traceback.print_exc()
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


def find_fresh_wrfout(path, dom_id):
    """
    Find the latest wrfout for postprocessing.

    :param path: the wrf path directory
    :param dom_id: the domain for which we search wrfouts
    :return: the path to the fresh (latest) wrfout
    """
    wrfouts = sorted(glob.glob(osp.join(path, 'wrfout_d%02d*' % dom_id)))
    return osp.join(path, wrfouts[-1])


def execute(args):
    """
    Executes a weather/fire simulation.

    The args dictionary contains

    :param args: a dictionary with the following keys
    :param grid_code: the (unique) code of the grid that is used
    :param sys_install_path: system installation directory
    :param start_utc: start time of simulation in UTC
    :param end_utc: end time of simulation in UTC
    :param workspace_path: workspace directory
    :param wps_install_path: installation directory of WPS that will be used
    :param wrf_install_path: installation directory of WRF that will be used
    :param grib_source: a string identifying a valid GRIB2 source
    :param wps_namelist_path: the path to the namelist.wps file that will be used as template
    :param wrf_namelist_path: the path to the namelist.input file that will be used as template
    :param fire_namelist_path: the path to the namelist.fire file that will be used as template
    :param wps_geog_path: the path to the geogrid data directory providing terrain/fuel data
    :param email_notification: dictionary containing keys address and events indicating when a mail should be fired off
    """

    # step 0 initialize the job state from the arguments
    js = JobState(args)

    logging.info("job %s starting [%d hours to forecast]." % (js.job_id, js.fc_hrs))
    send_email(js, 'start', 'Job %s started.' % js.job_id)

    # read in all namelists
    js.wps_nml = f90nml.read(js.args['wps_namelist_path'])
    js.wrf_nml = f90nml.read(js.args['wrf_namelist_path'])
    js.fire_nml = f90nml.read(js.args['fire_namelist_path'])
    js.ems_nml = None
    if 'emissions_namelist_path' in js.args:
        js.ems_nml = f90nml.read(js.args['emissions_namelist_path'])
    
    # Parse and setup the domain configuration
    js.domain_conf = WPSDomainConf(js.domains)

    num_doms = len(js.domain_conf)
    js.wps_nml['share']['start_date'] = [utc_to_esmf(js.start_utc)] * num_doms
    js.wps_nml['share']['end_date'] = [utc_to_esmf(js.end_utc)] * num_doms
    js.wps_nml['share']['interval_seconds'] = 3600

    logging.info("number of domains defined is %d." % num_doms)

    # build directories in workspace
    js.wps_dir = osp.abspath(osp.join(js.workspace_path, js.job_id, 'wps'))
    js.wrf_dir = osp.abspath(osp.join(js.workspace_path, js.job_id, 'wrf'))

    check_obj(args,'args')
    check_obj(js,'Initial job state')

    # step 1: clone WPS and WRF directories
    logging.info("cloning WPS into %s" % js.wps_dir)
    cln = WRFCloner(js.args)
    cln.clone_wps(js.wps_dir, js.grib_source.vtables(), [])

    # step 2: process domain information and patch namelist for geogrid
    js.wps_nml['geogrid']['geog_data_path'] = js.args['wps_geog_path']
    js.domain_conf.prepare_for_geogrid(js.wps_nml, js.wrf_nml, js.wrfxpy_dir, js.wps_dir)
    f90nml.write(js.wps_nml, osp.join(js.wps_dir, 'namelist.wps'), force=True)

    # do steps 2 & 3 & 4 in parallel (two execution streams)
    #  -> GEOGRID ->
    #  -> GRIB2 download ->  UNGRIB ->

    proc_q = Queue()
    geogrid_proc = Process(target=run_geogrid, args=(js, proc_q))
    grib_proc = Process(target=retrieve_gribs_and_run_ungrib, args=(js, proc_q))


    logging.info('starting GEOGRID and GRIB2/UNGRIB')
    geogrid_proc.start()
    grib_proc.start()

    # wait until both tasks are done
    logging.info('waiting until both tasks are done')
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
    if 'ignitions' in js.args:
        update_namelist(js.wrf_nml, render_ignitions(js, num_doms))

    # if we have an emissions namelist, automatically turn on the tracers
    if js.ems_nml is not None:
        logging.debug('namelist.fire_emissions given, turning on tracers')
        f90nml.write(js.ems_nml, osp.join(js.wrf_dir, 'namelist.fire_emissions'), force=True)
        js.wrf_nml['dynamics']['tracer_opt'] = [2] * num_doms

    f90nml.write(js.wrf_nml, osp.join(js.wrf_dir, 'namelist.input'), force=True)

    f90nml.write(js.fire_nml, osp.join(js.wrf_dir, 'namelist.fire'), force=True)

    # try to run Real twice as it sometimes fails the first time
    # it's not clear why this error happens 
    try:
        Real(js.wrf_dir).execute().check_output()
    except Exception as e:
        logging.error('Real step failed with exception %s, retrying ...' % str(e))
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
    WRF(js.wrf_dir, js.qsys).submit(js.task_id, js.num_nodes, js.ppn, js.wall_time_hrs)

    send_email(js, 'wrf_exec', 'Job %s - wrf job starting now with id %s.' % (js.job_id, js.task_id))
    logging.info("WRF job submitted with id %s, waiting for rsl.error.0000" % js.task_id)

    postprocess_output(js,args)

def postprocess_output(js,args):

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
    pp = None
    already_sent_files, max_pp_dom = [], -1
    if js.postproc is not None:
        js.pp_dir = osp.join(js.workspace_path, js.job_id, "products")
        make_dir(js.pp_dir)
        pp = Postprocessor(js.pp_dir, 'wfc-' + js.grid_code)
	max_pp_dom = max([int(x) for x in filter(lambda x: len(x) == 1, js.postproc)])

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
            if js.postproc is not None and str(dom_id) in js.postproc:
                var_list = [str(x) for x in js.postproc[str(dom_id)]]
                logging.info("Executing postproc instructions for vars %s for domain %d." % (str(var_list), dom_id))
                wrfout_path = find_fresh_wrfout(js.wrf_dir, dom_id)
            try:
                pp.process_vars(wrfout_path, dom_id, esmf_time, var_list)
            except Exception as e:
                logging.warning('Failed to postprocess for time %s with error %s.' % (esmf_time, str(e)))

            # if this is the last processed domain for this timestamp in incremental mode, upload to server
            if dom_id == max_pp_dom and js.postproc.get('shuttle', None) == 'incremental':
                desc = js.postproc['description'] if 'description' in js.postproc else js.job_id
                sent_files_1 = send_product_to_server(args, js.pp_dir, js.job_id, js.job_id, desc, already_sent_files)
                logging.info('sent %d files to visualization server.'  % len(sent_files_1))
                already_sent_files = filter(lambda x: not x.endswith('json'), already_sent_files + sent_files_1)

    # if we are to send out the postprocessed files after completion, this is the time
    if js.postproc.get('shuttle', None) == 'on_completion':
        desc = js.postproc['description'] if 'description' in js.postproc else js.job_id
        send_product_to_server(args, js.pp_dir, js.job_id, js.job_id, desc)


def verify_inputs(args,sys_cfg):
    """
    Check if arguments (eventually) supplied to execute(...) are valid - if not exception is thrown.

    Arguments:
      args -- dictionary of arguments
    """
    # dump(sys_cfg,'sys_cfg')
    # dump(args,'args')

    for key in sys_cfg:
        if key in args:
            if  sys_cfg[key] != args[key]:
               logging_error('system configuration %s=%s attempted change to %s' 
                   % (key, sys_cfg[key], args[key]))
               raise ValueError('System configuration values may not be overwritten.') 


    # we don't check if job_id is a valid path
    required_files = [('sys_install_path', 'Non-existent system installation directory %s'),
                      ('workspace_path', 'Non-existent workspace directory %s'),
                      ('wps_install_path', 'Non-existent WPS installation directory %s'),
                      ('wrf_install_path', 'Non-existent WRF installation directory %s'),
                      ('wps_namelist_path', 'Non-existent WPS namelist template %s'),
                      ('wrf_namelist_path', 'Non-existent WRF namelist template %s'),
                      ('fire_namelist_path', 'Non-existent fire namelist template %s'),
                      ('wps_geog_path', 'Non-existent geogrid data (WPS-GEOG) path %s')]

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
    if args['grib_source'] not in ['HRRR', 'NAM','NAM227', 'NARR']:
        raise ValueError('Invalid grib source, must be one of HRRR, NAM, NAM227, NARR')

    # if precomputed key is present, check files linked in
    if 'precomputed' in args:
      for key,path in args['precomputed'].iteritems():
          if not osp.exists(path):
              raise OSError('Precomputed entry %s points to non-existent file %s' % (key,path))

    # check if the postprocessor knows how to handle all variables
    wvs = get_wisdom_variables()
    failing = False
    if 'postproc' in args:
        for dom in filter(lambda x: len(x) == 1, args['postproc'].keys()):
            for vname in args['postproc'][dom]:
                if vname not in wvs:
                    logging.error('unrecognized variable %s in postproc key for domain %s.' % (vname, dom))
                    failing = True
    if failing:
        raise ValueError('One or more unrecognized variables in postproc.')


def process_arguments(args):
    """
    Convert arguments passed into program via the JSON configuration file.

    Transforms unicode strings into standard strings.

    :param args: the input arguments
    """
    # resolve possible relative time specifications
    start_utc = timespec_to_utc(args['start_utc'])
    args['orig_start_utc'] = start_utc
    args['start_utc'] = round_time_to_hour(start_utc)
    args['end_utc'] = round_time_to_hour(timespec_to_utc(args['end_utc'], args['start_utc']), True)

    for k, v in args.iteritems():
        if type(v) == unicode:
            args[k] = v.encode('ascii')


if __name__ == '__main__':

    # configure the basic logger
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    # logging.basicConfig(level=logging.DEBUG)

    # load the system configuration
    sys_cfg = None
    try:
        sys_cfg = json.load(open('etc/conf.json'))
    except IOError:
        logging.critical('Cannot find system configuration, have you created etc/conf.json?')
        sys.exit(2)

    # load configuration JSON
    # note: the execution flow allows us to override anything in the etc/conf.json file
    # dump(sys_cfg,'sys_cfg')
    job_args = json.load(open(sys.argv[1]), 'ascii')
    # dump(job_args,'job_args')
    args = sys_cfg
    args.update(job_args)

    process_arguments(args)

    # sanity check, also that nothing in etc/conf got overrident
    verify_inputs(args,sys_cfg)

    # execute the job
    execute(args)
    
    logging.info('forecast.py done')

