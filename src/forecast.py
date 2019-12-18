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
                  update_namelist, timedelta_hours, esmf_to_utc, render_ignitions, make_dir, \
                  timespec_to_utc, round_time_to_hour, Dict, dump, save, load, check_obj, \
                  make_clean_dir, process_create_time, load_sys_cfg, ensure_dir, move, \
                  json_join, number_minutes, serial_json, link2copy, append2file
from write_geogrid import write_table
from vis.postprocessor import Postprocessor
from vis.var_wisdom import get_wisdom_variables

from ingest.NAM218 import NAM218
from ingest.HRRR import HRRR
from ingest.NAM227 import NAM227
from ingest.CFSR import CFSR_P, CFSR_S
from ingest.NARR import NARR

from ingest.MODIS import Terra, Aqua
from ingest.VIIRS import SNPP

from fmda.fuel_moisture_da import assimilate_fm10_observations

from ssh_shuttle import send_product_to_server, ssh_command

import f90nml
from datetime import datetime, timedelta
import time, re, json, sys, logging
import os.path as osp
import os
import stat
from multiprocessing import Process, Queue
import glob
import netCDF4 as nc4
import shutil
import numpy as np

import smtplib
from email.mime.text import MIMEText

import traceback
import pprint
from cleanup import parallel_job_running, delete_visualization



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
        self.grib_source = self.resolve_grib_source(self.grib_source,args)
        self.satellite_source = self.resolve_satellite_source(args)
        logging.info('Simulation requested from %s to %s' % (str(self.start_utc), str(self.end_utc)))
        self.start_utc = round_time_to_hour(self.start_utc, up=False, period_hours=self.grib_source[0].period_hours);
        self.end_utc = round_time_to_hour(self.end_utc, up=True, period_hours=self.grib_source[0].period_hours);
        self.cycle_start_utc = round_time_to_hour(self.get('cycle_start_utc',None), period_hours=self.grib_source[0].cycle_hours);
        logging.info('Simulation times rounded  %s to %s' % (str(self.start_utc), str(self.end_utc)))
        self.fc_hrs = timedelta_hours(self.end_utc - self.start_utc)
        if 'job_id' in args:
            logging.info('job_id %s given in the job description' % args['job_id'])
            self.job_id = args['job_id']
        else:
            logging.warning('job_id not given, creating.')
            self.job_id = 'wfc-' + self.grid_code + '-' + utc_to_esmf(self.start_utc) + '-{0:02d}'.format(self.fc_hrs)
        if 'restart' in args:
            logging.info('restart %s given in the job description' % args['restart'])
            self.restart = args['restart']
        else:
            self.restart = False
            logging.info('restart not in arguments, default restart option %s' % self.restart)
        self.emails = self.parse_emails(args)
        self.domains = args['domains']
        self.ignitions = args.get('ignitions', None)
        self.fmda = self.parse_fmda(args)
        self.postproc = args['postproc']
        self.wrfxpy_dir = args['sys_install_path']
        if 'clean_dir' in args:
            self.clean_dir = args['clean_dir']
        else:
            self.clean_dir = True
        self.args = args
        logging.debug('JobState initialized: ' + str(self))



    def resolve_grib_source(self, gs_name, js):
        """
        Creates the right GribSource object from the name.

        :param gs_name: the name of the grib source
        :param js: configuration json
        """
        if gs_name == 'HRRR':
            return [HRRR(js)]
        elif gs_name == 'NAM' or gs_name == 'NAM218' :
            return [NAM218(js)]
        elif gs_name == 'NAM227':
            return [NAM227(js)]
        elif gs_name == 'NARR':
            return [NARR(js)]
        elif gs_name == 'CFSR':
            return [CFSR_P(js),CFSR_S(js)]
        else:
            raise ValueError('Unrecognized grib_source %s' % gs_name)

    def resolve_satellite_source(self, js):
        """
        Creates all the JPSSSource objects from the list of names.

        :param sat_list: the list of JPSS sources
        :param js: configuration json
        """
        sat_list = self.parse_satellite_source(js)
        sat = []
        if 'Terra' in sat_list:
            terra=Terra(js)
            sat.append(terra)
        if 'Aqua' in sat_list:
            aqua=Aqua(js)
            sat.append(aqua)
        if 'SNPP' in sat_list:
            snpp=SNPP(js)
            sat.append(snpp)
        return sat

    def parse_satellite_source(self, args):
        """
        Parse information inside the satellite source, if any.

        :param args: the forecast job argument dictionary
        """
        if 'satellite_source' in args:
            sats = args['satellite_source']
            return sats
        else:
            return []

    def parse_fmda(self, args):
        """
        Parse information inside the FMDA blob, if any.

        :param args: the forecast job argument dictionary
        """
        if 'fuel_moisture_da' in args:
            fmda = args['fuel_moisture_da']
            return Dict({'token' : fmda['mesowest_token'], 'domains' : fmda['domains']})
        else:
            return None


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

def retrieve_satellite(js, sat_source, q):
    """
    This function retrieves required Satellite files.

    It returns either 'SUCCESS' or 'FAILURE' on completion.

    :param js: the JobState object containing the forecast configuration
    :param sat_source: the SatSource object
    :param q: the multiprocessing Queue into which we will send either 'SUCCESS' or 'FAILURE'
    """
    try:
        logging.info("retrieving satellite files from %s" % sat_source.id)
        # retrieve satellite granules intersecting the last domain
        max_dom = max([int(x) for x in filter(lambda x: len(x) == 1, js.bounds.keys())])
        manifest = sat_source.retrieve_data_sat(js.bounds[str(max_dom)], js.start_utc, js.end_utc)
        # write a json file with satellite information
        sat_file = sat_source.id+'.json'
        json.dump(manifest, open(osp.join(js.jobdir,sat_file),'w'), indent=4, separators=(',', ': '))

        send_email(js, 'satellite', 'Job %s - satellite retrieving complete.' % js.job_id)
        logging.info('satellite retrieval complete for %s' % sat_source.id)
        q.put('SUCCESS')

    except Exception as e:
        logging.error('satellite retrieving step failed with exception %s' % repr(e))
        traceback.print_exc()
        q.put('FAILURE')

def retrieve_gribs_and_run_ungrib(js, grib_source, q):
    """
    This function retrieves required GRIB files and runs ungrib.

    It returns either 'SUCCESS' or 'FAILURE' on completion.

    :param js: the JobState object containing the forecast configuration
    :param grib_source: the GribSource object containing ungrib configuration
    :param q: the multiprocessing Queue into which we will send either 'SUCCESS' or 'FAILURE'
    """
    wps_dir = osp.abspath(js.wps_dir)
    grib_dir = osp.join(wps_dir,grib_source.id)
    make_clean_dir(grib_dir)
    wps_nml = js.wps_nml
    try:
        logging.info("retrieving GRIB files from %s" % grib_source.id)

        download_whole_cycle = js.get('download_whole_cycle',False)
        manifest = grib_source.retrieve_gribs(js.start_utc, js.end_utc, js.ref_utc, js.cycle_start_utc, download_whole_cycle)
        logging.info('manifest: ' + str(manifest))
        grib_file = grib_source.id+'.json'
        json.dump(manifest, open(osp.join(js.jobdir,grib_file),'w'), indent=4, separators=(',', ': '), default=serial_json)

        cache_colmet = len(manifest) > 1
        have_all_colmet = False
        if cache_colmet:
            have_all_colmet = len(manifest.colmet_missing) == 0
            colmet_dir = osp.join(grib_source.cache_dir, manifest.colmet_prefix)

        logging.info('cache colmet %s, have all colmet %s' % (cache_colmet, have_all_colmet))

        if not have_all_colmet:
            # this is also if we do not cache
            grib_source.symlink_gribs(manifest.grib_files, grib_dir)

            send_email(js, 'grib2', 'Job %s - %d GRIB2 files downloaded.' % (js.job_id, len(manifest)))
            logging.info("running UNGRIB for %s" % grib_source.id)

            logging.info("step 4: patch namelist for ungrib end execute ungrib on %s files" % grib_source.id)

            update_namelist(wps_nml, grib_source.namelist_wps_keys())
            if cache_colmet:
                wps_nml['share']['start_date'] = [utc_to_esmf(manifest.colmet_files_utc[0])] * js.num_doms
                wps_nml['share']['end_date'] = [utc_to_esmf(manifest.colmet_files_utc[-1])] * js.num_doms

            # logging.info("namelist.wps for UNGRIB: %s" % json.dumps(wps_nml, indent=4, separators=(',', ': ')))
            f90nml.write(wps_nml, osp.join(grib_dir, 'namelist.wps'), force=True)
            grib_source.clone_vtables(grib_dir)
            symlink_unless_exists(osp.join(wps_dir,'ungrib.exe'),osp.join(grib_dir,'ungrib.exe'))

            print(grib_dir + ':')
            os.system('ls -l %s' % grib_dir)

            Ungrib(grib_dir).execute().check_output()

            print(grib_dir + ':')
            os.system('ls -l %s' % grib_dir)

            if cache_colmet:
                # move output to cache directory
                make_dir(colmet_dir)
                for f in manifest.colmet_files:
                    move(osp.join(grib_dir,f),osp.join(colmet_dir,f))
                # now all colmet files should be in the cache

        if cache_colmet:
            for f in manifest.colmet_files:
                symlink_unless_exists(osp.join(colmet_dir,f),osp.join(wps_dir,f))
        else:
            # move output
            for f in glob.glob(osp.join(grib_dir,grib_source.prefix() + '*')):
                move(f,wps_dir)


        send_email(js, 'ungrib', 'Job %s - ungrib complete.' % js.job_id)
        logging.info('UNGRIB complete for %s' % grib_source.id)
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


def find_wrfout(path, dom_id, esmf_time):
    """
    Find wrfout for postprocessing.

    :param path: the wrf path directory
    :param dom_id: the domain for which we search wrfouts
    :esmf_time: time string to match variable Times
    :return: the path to the fresh (latest) wrfout
    """
    logging.info('find_wrfout: looking for the first wrfout for domain %s time %s' % (dom_id,esmf_time))
    wrfouts = sorted(glob.glob(osp.join(path, 'wrfout_d%02d*' % dom_id)),None,None,True) # reverse order
    wrfouts
    for wrfout in wrfouts:
        wrfout_time = re.match(r'.*wrfout_d.._([0-9_\-:]{19})' ,wrfout).groups()[0]
        if esmf_time >= wrfout_time:
            logging.info('find_wrfout: found %s' % wrfout)
            return wrfout
    logging.warning('wrfout for time %s domain %s not found' % (esmf_time, dom_id))
    logging.warning('Available wrfouts are: %s' % wrfouts)
    return None

def make_job_file(js):
    """
    Create minimal dictionary for the job state
    :param js: job state from JobState(args)
    :return: the dictionary
    """
    jsub=Dict({})
    jsub.job_id = js.job_id
    jsub.pid = os.getpid()
    jsub.process_create_time = process_create_time(jsub.pid)
    jsub.job_num = None
    jsub.old_job_num = None
    jsub.state = 'Preparing'
    jsub.qsys = js.qsys
    jsub.postproc = js.postproc
    jsub.grid_code = js.grid_code
    jsub.jobfile = osp.abspath(osp.join(js.workspace_path, js.job_id,'job.json'))
    return jsub

def make_kmz(args):
    ssh_command('wrfxweb/make_kmz.sh ' + args)

def read_namelist(path):
    logging.info('Reading namelist %s' % path)
    return f90nml.read(path)

def fmda_add_to_geogrid(js):
    """
    Add fmda datasets to geogrid if specified
    """
    if 'fmda_geogrid_path' in js:
        fmda_geogrid_path = osp.abspath(js['fmda_geogrid_path'])
        logging.info('fmda_geogrid_path is %s' % fmda_geogrid_path)
    else:
        return
        logging.info('fmda_geogrid_path not given')
    try:
        index_path = osp.join(fmda_geogrid_path,'index.json')
        index = json.load(open(index_path,'r'))
        logging.info('Loaded fmda geogrid index at %s' % index_path)
    except:
        logging.error('Cannot open %s' % index_path)
        sys.exit(1)
    geo_path = osp.dirname(osp.dirname(fmda_geogrid_path))+'-geo.nc'
    logging.info('fmda_add_to_geogrid reading longitudes and latitudes from NetCDF file %s' % geo_path )
    with nc4.Dataset(geo_path) as d:
        lats = d.variables['XLAT'][:,:]
        lons = d.variables['XLONG'][:,:]
    ndomains = len(js['domains'])
    lat,lon = js['domains'][str(ndomains)]['center_latlon']
    bbox = (np.min(lats), np.min(lons), np.max(lats), np.max(lons))
    logging.info('fmda_add_to_geogrid: fmda bounding box is %s %s %s %s' % bbox)
    i, j = np.unravel_index((np.abs(lats-lat)+np.abs(lons-lon)).argmin(),lats.shape)
    if i<=1 or j<=1 or i >= lats.shape[0]-2 or j >= lats.shape[1]-2:
        logging.error('fmda_add_to_geogrid: WRF domain center %s %s at %i %i is outside or near FMDA boundary' % (lat,lon,i,j) )
        sys.exit(1)
    """
    for varname,varindex in index.iteritems():
        # update index
        varindex['known_y']=float(i)
        varindex['known_x']=float(j)
        varindex['known_lat']=lats[i-1,j-1]
        varindex['known_lon']=lons[i-1,j-1]
        logging.info('fmda_add_to_geogrid: updating index known_x=%s known_y=%s known_lat=%s known_lon=%s' %
           (varindex['known_x'],varindex['known_y'],varindex['known_lat'],varindex['known_lon']))
        varindex_path=osp.join(fmda_geogrid_path,varname,'index')
        write_table(varindex_path,varindex)
    """
    # update geogrid table
    geogrid_tbl_path = osp.join(js.wps_dir, 'geogrid/GEOGRID.TBL')
    link2copy(geogrid_tbl_path)
    geogrid_tbl_json_path = osp.join(fmda_geogrid_path,'geogrid_tbl.json')
    logging.info('fmda_add_to_geogrid: updating GEOGRID.TBL at %s from %s' %
        (geogrid_tbl_path,geogrid_tbl_json_path))
    geogrid_tbl_json = json.load(open(geogrid_tbl_json_path,'r'))
    for varname,vartable in geogrid_tbl_json.iteritems():
        write_table(geogrid_tbl_path,vartable,mode='a',divider_after=True)

def execute(args,job_args):
    """
    Executes a weather/fire simulation.

    :param args: a dictionary with all to start the simulationfollowing keys
    :param job_args: a the original json given the forecast

    Keys in args:
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

    logging.info('step 0 initialize the job state from the arguments')
    js = JobState(args)

    jobdir = osp.abspath(osp.join(js.workspace_path, js.job_id))
    js.jobdir = jobdir
    if js.clean_dir or not osp.exists(osp.join(js.jobdir,'input.json')):
        make_clean_dir(js.jobdir)

    json.dump(job_args, open(osp.join(js.jobdir,'input.json'),'w'), indent=4, separators=(',', ': '))
    jsub = make_job_file(js)
    json.dump(jsub, open(jsub.jobfile,'w'), indent=4, separators=(',', ': '))

    logging.info("job %s starting [%d hours to forecast]." % (js.job_id, js.fc_hrs))
    sys.stdout.flush()
    send_email(js, 'start', 'Job %s started.' % js.job_id)

    # read in all namelists
    js.wps_nml = read_namelist(js.args['wps_namelist_path'])
    js.wrf_nml = read_namelist(js.args['wrf_namelist_path'])
    js.fire_nml = read_namelist(js.args['fire_namelist_path'])
    js.ems_nml = None
    if 'emissions_namelist_path' in js.args:
        js.ems_nml = read_namelist(js.args['emissions_namelist_path'])

    # Parse and setup the domain configuration
    js.domain_conf = WPSDomainConf(js.domains)

    js.num_doms = len(js.domain_conf)
    js.wps_nml['share']['interval_seconds'] = js.grib_source[0].interval_seconds

    logging.info("number of domains defined is %d." % js.num_doms)

    # build directories in workspace
    js.wps_dir = osp.abspath(osp.join(js.jobdir, 'wps'))
    js.wrf_dir = osp.abspath(osp.join(js.jobdir, 'wrf'))

    #check_obj(args,'args')
    #check_obj(js,'Initial job state')

    logging.info("step 1: clone WPS and WRF directories")
    logging.info("cloning WPS into %s" % js.wps_dir)
    cln = WRFCloner(js.args)
    cln.clone_wps(js.wps_dir, [])
    js.grib_source[0].clone_vtables(js.wps_dir)

    logging.info("step 2: process domain information and patch namelist for geogrid")
    js.wps_nml['share']['start_date'] = [utc_to_esmf(js.start_utc)] * js.num_doms
    js.wps_nml['share']['end_date'] = [utc_to_esmf(js.end_utc)] * js.num_doms
    js.wps_nml['geogrid']['geog_data_path'] = js.args['wps_geog_path']
    js.domain_conf.prepare_for_geogrid(js.wps_nml, js.wrf_nml, js.wrfxpy_dir, js.wps_dir)
    f90nml.write(js.wps_nml, osp.join(js.wps_dir, 'namelist.wps'), force=True)

    # do steps 2 & 3 & 4 in parallel (three execution streams)
    #  -> Satellite retrieval ->
    #  -> GEOGRID ->
    #  -> GRIB2 download ->  UNGRIB ->

    proc_q = Queue()

    if js.satellite_source:
        js.bounds = Dict({})
        for k,domain in enumerate(js.domain_conf.domains):
            latloni = domain.ij_to_latlon(0,0)
            latlonf = domain.ij_to_latlon(domain.domain_size[0],domain.domain_size[1])
            bounds = (latloni[1],latlonf[1],latloni[0],latlonf[0])
            js.bounds[str(k+1)] = bounds
        sat_proc = {}
    for satellite_source in js.satellite_source:
        sat_proc[satellite_source.id] = Process(target=retrieve_satellite, args=(js, satellite_source, proc_q))

    geogrid_proc = Process(target=run_geogrid, args=(js, proc_q))

    grib_proc = {}
    for grib_source in js.grib_source:
        grib_proc[grib_source.id] = Process(target=retrieve_gribs_and_run_ungrib, args=(js, grib_source, proc_q))

    logging.info('starting GEOGRID and GRIB2/UNGRIB')
    logging.info('satellite sources %s' % [s.id for s in js.satellite_source])

    if js.ungrib_only:
        logging.info('ungrib_only set, skipping GEOGRID and SATELLITE, will exit after UNGRIB')
    else:
        geogrid_proc.start()
        for satellite_source in js.satellite_source:
            sat_proc[satellite_source.id].start()

    for grib_source in js.grib_source:
        grib_proc[grib_source.id].start()

    # wait until all tasks are done
    logging.info('waiting until all tasks are done')

    for grib_source in js.grib_source:
        grib_proc[grib_source.id].join()

    if js.ungrib_only:
        return
    else:
        geogrid_proc.join()
        for satellite_source in js.satellite_source:
            sat_proc[satellite_source.id].join()

    for satellite_source in js.satellite_source:
        if proc_q.get() != 'SUCCESS':
            return

    for grib_source in js.grib_source:
        if proc_q.get() != 'SUCCESS':
            return

    if proc_q.get() != 'SUCCESS':
        return

    proc_q.close()

    logging.info('execute: finished parallel GEOGRID, GRIB2/UNGRIB, and Satellite')

    if js.satellite_source:
        # joining all satellite manifests
        sat_manifest = Dict({})
        sat_manifest.granules = json_join(js.jobdir, args['satellite_source'])
        sat_manifest.bounds = js.bounds
        sat_manifest.time_interval = (utc_to_esmf(js.start_utc), utc_to_esmf(js.end_utc))
        sat_manifest.dt = Dict({})
        sat_manifest.sat_interval = Dict({})
        for k in args['domains'].keys():
            sat_manifest.dt[k] = args['domains'][k]['history_interval']
            if 'sat_interval' in args['domains'][k].keys():
                sat_manifest.sat_interval[k] = args['domains'][k]['sat_interval']
            else:
                sat_manifest.sat_interval[k] = (args['domains'][k]['history_interval'], args['domains'][k]['history_interval'])
        sat_manifest.satprod_satsource = js.satprod_satsource
        json.dump(sat_manifest, open(osp.join(js.jobdir, 'sat.json'),'w'), indent=4, separators=(',', ': '))

    logging.info("step 5: execute metgrid after ensuring all grids will be processed")
    update_namelist(js.wps_nml, js.grib_source[0].namelist_wps_keys())
    js.domain_conf.prepare_for_metgrid(js.wps_nml)
    logging.info("namelist.wps for METGRID: %s" % json.dumps(js.wps_nml, indent=4, separators=(',', ': ')))
    f90nml.write(js.wps_nml, osp.join(js.wps_dir, 'namelist.wps'), force=True)

    logging.info("running METGRID")
    Metgrid(js.wps_dir).execute().check_output()

    send_email(js, 'metgrid', 'Job %s - metgrid complete.' % js.job_id)
    logging.info("METGRID complete")

    logging.info("cloning WRF into %s" % js.wrf_dir)

    logging.info("step 6: clone wrf directory, symlink all met_em* files, make namelists")
    cln.clone_wrf(js.wrf_dir, [])
    symlink_matching_files(js.wrf_dir, js.wps_dir, "met_em*")
    time_ctrl = update_time_control(js.start_utc, js.end_utc, js.num_doms)
    js.wrf_nml['time_control'].update(time_ctrl)
    js.wrf_nml['time_control']['interval_seconds'] = js.grib_source[0].interval_seconds
    update_namelist(js.wrf_nml, js.grib_source[0].namelist_keys())
    if 'ignitions' in js.args:
        update_namelist(js.wrf_nml, render_ignitions(js, js.num_doms))

    # if we have an emissions namelist, automatically turn on the tracers
    if js.ems_nml is not None:
        logging.debug('namelist.fire_emissions given, turning on tracers')
        f90nml.write(js.ems_nml, osp.join(js.wrf_dir, 'namelist.fire_emissions'), force=True)
        js.wrf_nml['dynamics']['tracer_opt'] = [2] * js.num_doms

    f90nml.write(js.wrf_nml, osp.join(js.wrf_dir, 'namelist.input'), force=True)

    f90nml.write(js.fire_nml, osp.join(js.wrf_dir, 'namelist.fire'), force=True)

    # step 7: execute real.exe

    logging.info("running REAL")
    # try to run Real twice as it sometimes fails the first time
    # it's not clear why this error happens
    try:
        Real(js.wrf_dir).execute().check_output()
    except Exception as e:
        logging.error('Real step failed with exception %s, retrying ...' % str(e))
        Real(js.wrf_dir).execute().check_output()


    logging.info('step 7b: if requested, do fuel moisture DA')
    logging.info('fmda = %s' % js.fmda)
    if js.fmda is not None:
        logging.info('running fuel moisture data assimilation')
        for dom in js.fmda.domains:
            logging.info('assimilate_fm10_observations for domain %s' % dom)
            assimilate_fm10_observations(osp.join(js.wrf_dir, 'wrfinput_d%02d' % int(dom)), None, js.fmda.token)

    # step 8: execute wrf.exe on parallel backend
    logging.info('submitting WRF job')
    send_email(js, 'wrf_submit', 'Job %s - wrf job submitted.' % js.job_id)

    js.task_id = "sim-" + js.grid_code + "-" + utc_to_esmf(js.start_utc)[:10]
    jsub.job_num=WRF(js.wrf_dir, js.qsys).submit(js.task_id, js.num_nodes, js.ppn, js.wall_time_hrs)

    send_email(js, 'wrf_exec', 'Job %s - wrf job starting now with id %s.' % (js.job_id, js.task_id))
    logging.info("WRF job %s submitted with id %s, waiting for rsl.error.0000" % (jsub.job_num, js.task_id))

    jobfile = osp.abspath(osp.join(js.jobdir, 'job.json'))
    json.dump(jsub, open(jobfile,'w'), indent=4, separators=(',', ': '))

    process_output(js.job_id)

    return jobfile

def process_output(job_id):
    args = load_sys_cfg()
    jobfile = osp.abspath(osp.join(args.workspace_path, job_id,'job.json'))
    satfile = osp.abspath(osp.join(args.workspace_path, job_id,'sat.json'))
    logging.info('process_output: loading job description from %s' % jobfile)
    try:
        js = Dict(json.load(open(jobfile,'r')))
    except Exception as e:
        logging.error('Cannot load the job description file %s' % jobfile)
        logging.error('%s' % e)
        sys.exit(1)
    logging.info('process_output: loading satellite description from %s' % satfile)
    try:
        jsat = Dict(json.load(open(satfile,'r')))
        available_sats = [sat.upper()+'_AF' for sat in jsat.granules.keys()]
        not_empty_sats = [sat.upper()+'_AF' for sat in jsat.granules.keys() if jsat.granules[sat]]
    except:
        logging.warning('Cannot load the satellite data in satellite description file %s' % satfile)
        available_sats = []
        not_empty_sats = []
        pass
    logging.info('process_output: available satellite data %s' % available_sats)
    logging.info('process_output: not empty satellite data %s' % not_empty_sats)
    js.old_pid = js.pid
    js.pid = os.getpid()
    js.state = 'Processing'
    json.dump(js, open(jobfile,'w'), indent=4, separators=(',', ': '))

    js.wrf_dir = osp.abspath(osp.join(args.workspace_path, js.job_id, 'wrf'))

    pp = None
    already_sent_files = []
    if js.postproc is None:
        logging.info('No postprocessing specified, exiting.')
        return

    # set up postprocessing
    delete_visualization(js.job_id)
    js.pp_dir = osp.join(args.workspace_path, js.job_id, "products")
    make_clean_dir(js.pp_dir)
    pp = Postprocessor(js.pp_dir, 'wfc-' + js.grid_code)
    js.manifest_filename= 'wfc-' + js.grid_code + '.json'
    logging.debug('Postprocessor created manifest %s',js.manifest_filename)

    if js.postproc.get('from', None) == 'wrfout':
        logging.info('Postprocessing all wrfout files.')
        failures = cases = 0
        # postprocess all wrfouts
        for wrfout_path in sorted(glob.glob(osp.join(js.wrf_dir,'wrfout_d??_????-??-??_??:??:??'))):
            logging.info("Found %s" % wrfout_path)
            domain_str,wrfout_esmf_time = re.match(r'.*wrfout_d(0[0-9])_([0-9_\-:]{19})',wrfout_path).groups()
            dom_id = int(domain_str)
            d = nc4.Dataset(wrfout_path)
            # extract ESMF string times
            times = [''.join(x) for x in d.variables['Times'][:]]
            d.close()
            for esmf_time in sorted(times):
                logging.info("Processing domain %d for time %s." % (dom_id, esmf_time))
                if js.postproc is not None and str(dom_id) in js.postproc:
                    cases += 1
                    if available_sats:
                        sat_list = [sat for sat in available_sats if sat in js.postproc[str(dom_id)]]
                        var_list = [str(x) for x in js.postproc[str(dom_id)] if not str(x) in sat_list]
                        sat_list = [sat for sat in sat_list if sat in not_empty_sats]
                        logging.info("Executing postproc instructions for sats %s for domain %d." % (str(sat_list), dom_id))
                    else:
                        sat_list = []
                        var_list = [str(x) for x in js.postproc[str(dom_id)]]
                    logging.info("Executing postproc instructions for vars %s for domain %d." % (str(var_list), dom_id))
                    try:
                        if sat_list:
                            pp.process_sats(jsat, dom_id, esmf_time, sat_list)
                        pp.process_vars(osp.join(js.wrf_dir,wrfout_path), dom_id, esmf_time, var_list)
                        # in incremental mode, upload to server
                        if js.postproc.get('shuttle', None) == 'incremental':
                            desc = js.postproc['description'] if 'description' in js.postproc else js.job_id
                            sent_files_1 = send_product_to_server(args, js.pp_dir, js.job_id, js.job_id, js.manifest_filename, desc, already_sent_files)
                            already_sent_files = filter(lambda x: not x.endswith('json'), already_sent_files + sent_files_1)
                    except Exception as e:
                        logging.warning('Failed to postprocess for time %s with error %s.' % (esmf_time, str(e)))
                        failures += 1

        if cases != failures:
            logging.info('number of postprocessing steps is %d and number of postprocessing failures is %d' % (cases,failures))
            # if we are to send out the postprocessed files after completion, this is the time
            if js.postproc.get('shuttle', None) == 'on_completion':
                desc = js.postproc['description'] if 'description' in js.postproc else js.job_id
                send_product_to_server(args, js.pp_dir, js.job_id, js.job_id, js.manifest_filename, desc)

        else:
            logging.error('All postprocessing steps failed')
            js.state = 'Postprocessing failed'

        json.dump(js, open(jobfile,'w'), indent=4, separators=(',', ': '))
        return

    # step 9: wait for appearance of rsl.error.0000 and open it
    wrf_out = None
    rsl_path = osp.join(js.wrf_dir, 'rsl.error.0000')
    while wrf_out is None:
        try:
            wrf_out = open(rsl_path)
            break
        except IOError:
            logging.info('process_output: waiting 5 seconds for rsl.error.0000 file')
        time.sleep(5)

    logging.info('process_output: Detected rsl.error.0000')
    js.run_utc = time.ctime(os.path.getmtime(rsl_path))
    js.processed_utc = time.asctime(time.gmtime())

    # step 10: track log output and check for history writes fro WRF
    wait_lines = 0
    wait_wrfout = 0
    failures = cases = 0
    while True:
        line = wrf_out.readline().strip()
        if not line:
            if not parallel_job_running(js):
                logging.warning('WRF did not run to completion.')
                break
            if not wait_lines:
                logging.info('Waiting for more output lines')
            wait_lines = wait_lines + 1
            time.sleep(5)
            continue
        wait_lines = 0

        if "SUCCESS COMPLETE WRF" in line:
            # send_email(js, 'complete', 'Job %s - wrf job complete SUCCESS.' % js.job_id)
            logging.info("WRF completion detected.")
            js.old_job_num = js.job_num
            js.job_num = None
            json.dump(js, open(jobfile,'w'), indent=4, separators=(',', ': '))
            break

        if "Timing for Writing wrfout" in line:
            wait_wrfout = 0
            esmf_time,domain_str = re.match(r'.*wrfout_d.._([0-9_\-:]{19}) for domain\ +(\d+):' ,line).groups()
            wrfout_path,domain_str = re.match(r'.*(wrfout_d.._[0-9_\-:]{19}) for domain\ +(\d+):' ,line).groups()
            dom_id = int(domain_str)
            logging.info("Detected history write for domain %d for time %s." % (dom_id, esmf_time))
            if js.postproc is not None and str(dom_id) in js.postproc:
                cases += 1
                if available_sats:
                    sat_list = [sat for sat in available_sats if sat in js.postproc[str(dom_id)]]
                    var_list = [str(x) for x in js.postproc[str(dom_id)] if not str(x) in sat_list]
                    sat_list = [sat for sat in sat_list if sat in not_empty_sats]
                    logging.info("Executing postproc instructions for sats %s for domain %d." % (str(sat_list), dom_id))
                else:
                    sat_list = []
                    var_list = [str(x) for x in js.postproc[str(dom_id)]]
                logging.info("Executing postproc instructions for vars %s for domain %d." % (str(var_list), dom_id))
                wrfout_path = find_wrfout(js.wrf_dir, dom_id, esmf_time)
                try:
                    if sat_list:
                        pp.process_sats(jsat, dom_id, esmf_time, sat_list)
                    pp.process_vars(osp.join(js.wrf_dir,wrfout_path), dom_id, esmf_time, var_list)
                except Exception as e:
                    logging.warning('Failed to postprocess for time %s with error %s.' % (esmf_time, str(e)))
                    failures += 1
                else:
                    try:
                        # in incremental mode, upload to server
                        if js.postproc.get('shuttle', None) == 'incremental':
                            desc = js.postproc['description'] if 'description' in js.postproc else js.job_id
                            sent_files_1 = send_product_to_server(args, js.pp_dir, js.job_id, js.job_id, js.manifest_filename, desc, already_sent_files)
                            already_sent_files = filter(lambda x: not x.endswith('json'), already_sent_files + sent_files_1)
                    except Exception as e:
                        logging.warning('Failed sending potprocess results to the server with error %s' % str(e))
        else:
            if not wait_wrfout:
                logging.info('Waiting for wrfout')
            wait_wrfout = wait_wrfout + 1

    # if we are to send out the postprocessed files after completion, this is the time
    if cases != failures:
        logging.info('number of postprocessing steps is %d and number of postprocessing failures is %d' % (cases,failures))
        if js.postproc.get('shuttle', None) == 'on_completion':
            desc = js.postproc['description'] if 'description' in js.postproc else js.job_id
            send_product_to_server(args, js.pp_dir, js.job_id, js.job_id, js.manifest_filename, desc)

        if js.postproc.get('shuttle', None) is not None:
            make_kmz(js.job_id)  # arguments can be added to the job id string

        js.state = 'Completed'

    else:
        logging.error('All postprocessing steps failed')
        js.state = 'Postprocessing failed'

    js.old_pid = js.pid
    js.pid = None
    json.dump(js, open(jobfile,'w'), indent=4, separators=(',', ': '))

def create_process_output_script(job_id):
    cfg = load_sys_cfg()
    script_path = osp.join(cfg.workspace_path, job_id,'job_process_output.sh')
    log_path = osp.join(cfg.workspace_path, job_id,'job_process_output.log')
    process_script = osp.join(cfg.sys_install_path,'process_output.sh')
    with open(script_path,'w') as f:
        f.write('#!/usr/bin/env bash\n')
        f.write('cd ' + cfg.sys_install_path + '\n')
        f.write('LOG=' + log_path + '\n')
        f.write(process_script + ' ' + job_id + ' &> $LOG \n')

    # make it executable
    st = os.stat(script_path)
    os.chmod(script_path, st.st_mode | stat.S_IEXEC)


def process_sat_output(job_id):
    args = load_sys_cfg()
    jobfile = osp.abspath(osp.join(args.workspace_path, job_id,'job.json'))
    satfile = osp.abspath(osp.join(args.workspace_path, job_id,'sat.json'))
    logging.info('process_output: loading job description from %s' % jobfile)
    try:
        js = Dict(json.load(open(jobfile,'r')))
    except Exception as e:
        logging.error('Cannot load the job description file %s' % jobfile)
        logging.error('%s' % e)
        sys.exit(1)
    logging.info('process_output: loading satellite description from %s' % satfile)
    try:
        jsat = Dict(json.load(open(satfile,'r')))
        available_sats = [sat.upper()+'_AF' for sat in jsat.granules.keys()]
        not_empty_sats = [sat.upper()+'_AF' for sat in jsat.granules.keys() if jsat.granules[sat]]
    except:
        logging.warning('Cannot load the satellite data in satellite description file %s' % satfile)
        available_sats = []
        not_empty_sats = []
        return
    logging.info('process_output: available satellite data %s' % available_sats)
    logging.info('process_output: not empty satellite data %s' % not_empty_sats)
    if not not_empty_sats:
        logging.warning('Do not have satellite data to postprocess')
        return
    js.old_pid = js.pid
    js.pid = os.getpid()
    js.state = 'Processing'
    json.dump(js, open(jobfile,'w'), indent=4, separators=(',', ': '))

    pp = None
    already_sent_files = []
    delete_visualization(js.job_id)
    js.pp_dir = osp.join(args.workspace_path, js.job_id, "products")
    make_clean_dir(js.pp_dir)
    pp = Postprocessor(js.pp_dir, 'wfc-' + js.grid_code)
    js.manifest_filename= 'wfc-' + js.grid_code + '.json'
    logging.debug('Postprocessor created manifest %s',js.manifest_filename)
    domains = sorted([int(x) for x in filter(lambda x: len(x) == 1, js.postproc)])
    for dom_id in domains:
        logging.info('Processing domain %s' % str(dom_id))
        dt = timedelta(minutes=jsat.dt[str(dom_id)])
        logging.info('dt for satellite postprocessing = %s' % dt)
        t_int = esmf_to_utc(jsat['time_interval'][0])
        t_fin = esmf_to_utc(jsat['time_interval'][1])
        ndts = number_minutes(t_int,t_fin,jsat.dt[str(dom_id)])
        times = [t_int + tt*dt for tt in range(ndts)]
        sat_list = [sat for sat in available_sats if sat in js.postproc[str(dom_id)]]
        if sat_list:
            for time in times:
                try:
                    esmf_time = utc_to_esmf(time)
                    logging.info('Posprocessing satellite data for time %s' % esmf_time)
                    pp.process_sats(jsat, dom_id, esmf_time, sat_list)
                except Exception as e:
                    logging.warning('Failed to postprocess for time %s with error %s.' % (esmf_time, str(e)))
                else:
                    try:
                        if js.postproc.get('shuttle', None) == 'incremental':
                            desc = js.postproc['description'] if 'description' in js.postproc else js.job_id
                            sent_files_1 = send_product_to_server(args, js.pp_dir, js.job_id, js.job_id, js.manifest_filename, desc, already_sent_files)
                            already_sent_files = filter(lambda x: not x.endswith('json'), already_sent_files + sent_files_1)
                    except Exception as e:
                        logging.warning('Failed sending potprocess results to the server with error %s' % str(e))

    # if we are to send out the postprocessed files after completion, this is the time
    if js.postproc.get('shuttle', None) == 'on_completion':
        desc = js.postproc['description'] if 'description' in js.postproc else js.job_id
        send_product_to_server(args, js.pp_dir, js.job_id, js.job_id, js.manifest_filename, desc)

    if js.postproc.get('shuttle', None) is not None:
        make_kmz(js.job_id)  # arguments can be added to the job id string

    js.old_pid = js.pid
    js.pid = None
    js.state = 'Completed'
    json.dump(js, open(jobfile,'w'), indent=4, separators=(',', ': '))


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
    if args['grib_source'] not in ['HRRR', 'NAM','NAM218', 'NAM227', 'NARR','CFSR']:
        raise ValueError('Invalid grib source %s, must be one of HRRR, NAM, NAM227, NARR, CFSR' % args['grib_source'])

    # check for valid satellite source
    if 'satellite_source' in args:
        for sat in args['satellite_source']:
            if sat not in ['Terra','Aqua','SNPP']:
                raise ValueError('Invalid satellite source %s, must be one of Terra, Aqua, SNPP' % sat)

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


def process_arguments(job_args,sys_cfg):
    """
    Convert arguments passed into program via the JSON configuration file and job json argument.
    This is processed after the configuration is updated by the job json file.

    Transforms unicode strings into standard strings.

    :param args: the input arguments
    """
    # note: the execution flow allows us to override anything in the etc/conf.json file
    # dump(sys_cfg,'sys_cfg')
    args = sys_cfg
    keys = job_args.keys()
    for key in keys:
        if job_args[key] is None:
            logging.warning('Job argument %s=None, ignoring' % key)
            del job_args[key]
    args.update(job_args)
    # logging.info('updated args = %s' % json.dumps(args, indent=4, separators=(',', ': ')))

    # resolve possible relative time specifications
    start_utc = timespec_to_utc(args['start_utc'])
    args['orig_start_utc'] = start_utc
    args['start_utc'] = round_time_to_hour(start_utc)
    args['end_utc'] = round_time_to_hour(timespec_to_utc(args['end_utc'], args['start_utc']), True)
    args['cycle_start_utc'] = timespec_to_utc(args.get('cycle_start_utc', None))
    args['satprod_satsource'] = Dict({})

    # add postprocess satellite data
    if 'satellite_source' in args:
        if 'postproc' in args:
            max_dom = max([int(x) for x in filter(lambda x: len(x) == 1, args['postproc'])])
            sats = args['satellite_source']
            for sat in sats:
                satprod = sat.upper()+'_AF'
                args['postproc'][str(max_dom)].append(satprod)
                args['satprod_satsource'].update({satprod: sat})

    # defaults
    if args['ref_utc'] is not None:
        args['ref_utc'] = timespec_to_utc(args['ref_utc'], args['start_utc'])

    for k, v in args.iteritems():
        if type(v) == unicode:
            args[k] = v.encode('ascii')

    # sanity check, also that nothing in etc/conf got overrident
    verify_inputs(args,sys_cfg)

    return args

if __name__ == '__main__':

    # configure the basic logger
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    # logging.basicConfig(level=logging.DEBUG)

    # load configuration JSON
    sys_cfg = load_sys_cfg()
    # logging.info('sys_cfg = %s' % json.dumps(sys_cfg, indent=4, separators=(',', ': ')))

    # load job JSON
    job_args = json.load(open(sys.argv[1]), 'ascii')
    # logging.info('job_args = %s' % json.dumps(job_args, indent=4, separators=(',', ': ')))

    # process arguments
    args = process_arguments(job_args,sys_cfg)
    # logging.info('processed args = %s' % str(args))

    # execute the job
    logging.info('calling execute')
    execute(args,job_args)

    logging.info('forecast.py done')

