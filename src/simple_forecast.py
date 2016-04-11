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

import json
import sys
import utils
from datetime import timedelta, datetime
import pytz
import socket


def print_header(header):
    print('\033[31m** %s **\033[0m' % header)

def print_question(question):
    print('\033[95m%s\033[0m' % question)

def print_answer(ans):
    print('\033[1;33m-> %s\033[0m' % ans)

def read_string(default):
    read_in = sys.stdin.readline().strip('\n')
    return read_in if read_in != '' else default
    
def read_integer(default):
    return int(read_string(default))


def read_location(default):
    return map(float, read_string(default).split(','))


def read_size(default):
    return map(int, read_string(default).split(','))


def select_grib_source(start_time):
    now = datetime.utcnow().replace(tzinfo=pytz.UTC)
    if now - start_time < timedelta(days=30):
        return 'NAM'
    else:
        return 'NARR'


def read_time_indicator(default):
    s = read_string(default)
    if len(s)==0:
        return utils.esmf_to_utc(default).replace(tzinfo=pytz.UTC)
    else:
        return utils.timespec_to_utc(s).replace(tzinfo=pytz.UTC)

def read_boolean(default):
    s = read_string(default)
    return s.lower() == 'yes'


def queuing_systems():
    return map(str, json.load(open('etc/clusters.json')).keys())


def questionnaire():
    """
    Give a questionnaire to the user (with sensible default) to create
    a simple domain configuration.
    
    :return: a dictionary with the configuration of a fire simulation
    """
    cfg = {}
    print_question('Enter a name for your job [default = experiment]:')
    cfg['grid_code'] = read_string('experiment')
    print_answer('Name is %s' % cfg['grid_code'])

    print_header('IGNITION section')
    newline()

    print_question('Enter the ignition point as lat, lon [default = 39.1, -104.3]:')
    ign_latlon = read_location('39.1, -104.3')

    print_question('Enter the ignition time in UTC timezone as an ESMF string or relative time')
    print('Examples: 2016-03-30_16:00:00 or T-60 or T+30 (minutes), [default = T-60]')
    ign_utc = read_time_indicator('T-60') 
    print('Ignition time is %s\n' % str(ign_utc))

    print_question('Enter the duration of the ignition process in seconds [default = 240]')
    ign_dur = read_integer('240')

    newline()
    print_header('SIMULATION section')

    start_utc = None
    while True:
        print_question('Enter the start time of the simulation in UTC timezone [default = 30 mins before ignition time]')
        start_utc = read_time_indicator(utils.utc_to_esmf(ign_utc - timedelta(minutes=30)))
        if start_utc < ign_utc:
            break
        print('Simulation start must be before ignition time %s' % utils.utc_to_esmf(ign_utc))
    cfg['start_utc'] = utils.utc_to_esmf(start_utc)
    
    end_utc = None 
    while True:
        print_question('Enter the end time of the simulation [default = start_time + 5 hours]')
        end_utc= read_time_indicator(utils.utc_to_esmf(start_utc + timedelta(hours=5)))
        if end_utc > ign_utc:
            break
        print('Simulation end must be after ignition time %s' % utils.utc_to_esmf(ign_utc))
    cfg['end_utc'] = utils.utc_to_esmf(end_utc)

    print_question('Please enter the cell size in meters for the atmospheric mesh [default 1000]')
    cell_size = read_integer('1000')

    print_question('Enter the number of grid cells in the longitudinal and latitudinal position [default 61, 61]')
    domain_size = read_size('61, 61')

    print_question('Enter the refinement ratio for the fire grid [default=40]')
    refinement = read_integer('40')

    print_question('Enter the interval between output frames in minutes [default=15]')
    history_interval = read_integer('15')

    cfg['grib_source'] = select_grib_source(start_utc)
    print_answer('Selected GRIB2 source %s' % cfg['grib_source'])

    
    def_geog_path = None
    try:
        def_geog_path = json.load(open('etc/conf.json'))['default_geogrid_path']
    except Exception as e:
        print(e)
        pass
    print_question('Enter the path to geogrid information (WPS-GEOG) [default=%s]' % def_geog_path)
    cfg['geogrid_path'] = read_string(def_geog_path)

    cfg['domains'] = { '1' : {
        'cell_size' : (cell_size,cell_size),
        'domain_size' : domain_size,
        'subgrid_ratio' : (refinement, refinement),
        'geog_res' : '.3s',
        'center_latlon' : ign_latlon,
        'truelats' : (ign_latlon[0], ign_latlon[0]),
        'stand_lon' : ign_latlon[1],
        'history_interval' : history_interval,
        'time_step' : max(1, 5 * cell_size / 1000)
        }
    }

    cfg['ignitions'] = { '1' : [ { 'start_delay_s' : (ign_utc - start_utc).total_seconds(),
                                   'duration_s' : ign_dur,
                                   'latlon' : ign_latlon } ] }


    print_header('NAMELISTS')

    print_question('WPS [default etc/nlists/default.wps]')
    cfg['wps_namelist_path'] = read_string('etc/nlists/default.wps')

    print_question('Input [default etc/nlists/default.input]')
    cfg['wrf_namelist_path'] = read_string('etc/nlists/default.input')

    print_question('Fire [default etc/nlists/default.fire]')
    cfg['fire_namelist_path'] = read_string('etc/nlists/default.fire')

    print_question('Emissions [default etc/nlists/default.fire_emissions]')
    cfg['emissions_namelist_path'] = read_string('etc/nlists/default.fire_emissions')

    print_header('PARALLEL JOB configuration')
    print_question('Enter number of parallel nodes [default=8]')
    cfg['num_nodes'] = read_integer('8')

    print_question('Enter number of cores per node [default=12]')
    cfg['ppn'] = read_integer('12')

    print_question('Enter the max walltime in hours [default=2]')
    cfg['wall_time_hrs'] = read_integer('2')

    qsys_opts = queuing_systems()
    while True:
        def_qsys = socket.gethostname().split('.')[0]
        print('Enter queuing system [choices are %s, default is %s]' % (qsys_opts, def_qsys))
        cfg['qsys'] = read_string(def_qsys)
        if cfg['qsys'] in qsys_opts:
            break
        print('Invalid queuing system selected, please try again')

    print_header('POSTPROCESSING')
    print_question('Which variables should wrfxpy postprocess? [default T2,PSFC,WINDSPD,WINDVEC,FIRE_AREA,FGRNHFX,FLINEINT,SMOKE_INT]')
    pp_vars = read_string('T2,PSFC,WINDSPD,WINDVEC,FIRE_AREA,FGRNHFX,FLINEINT,SMOKE_INT').split(',')

    print_question('Send variables to visualization server? [default=no]')
    shuttle = read_boolean('no')

    desc = ''
    if shuttle:
        print_question('Enter a short description of your job [default=experimental run]')
        desc = read_string('experimental run')
    
    cfg['postproc'] = { '1' : pp_vars }
    if shuttle:
        cfg['postproc']['shuttle'] = 'incremental'
        cfg['postproc']['description'] = desc

    return cfg


def newline():
    print('')


if __name__ == '__main__':
    print_header('INT wrfxpy interactive fire simulation setup')
    print('This script will build a single domain simulation configuration')
    print('centered around the ignition point.')
    newline()

    cfg = questionnaire()

    filename = cfg['grid_code'] + '.json'
    json.dump(cfg, open(filename, 'w'), indent=4, separators=(',', ': '))

    newline()

    print('INT to start the simulation, execute ./forecast.sh %s' % filename)

