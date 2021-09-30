#
# Angel Farguell, CU Denver
#

from __future__ import absolute_import
from __future__ import print_function
from ingest.MODIS import Terra, Aqua
from ingest.VIIRS import SNPP, SNPPHR, NOAA20
from ingest.VIIRSLL import SNPPLL,J01LL
from wrf.wps_domains import WPSDomainConf
from utils import load_sys_cfg, esmf_to_utc, Dict
import datetime as dt
import logging, sys, json, pytz

if __name__ == '__main__':
    if len(sys.argv) == 2:
        # load input JSON
        try:
            js_input = Dict(json.load(open(sys.argv[1])))
        except IOError:
            import sys
            logging.critical('Cannot find input json file.')
            sys.exit(2)

        # inputs
        sat_sources = list(js_input["satellite_source"])
        domain = WPSDomainConf(js_input["domains"]).domains[-1]
        latloni = domain.ij_to_latlon(0,0)
        latlonf = domain.ij_to_latlon(domain.domain_size[0],domain.domain_size[1])
        bounds = (latloni[1],latlonf[1],latloni[0],latlonf[0])
        from_utc = esmf_to_utc(js_input["start_utc"])
        to_utc = esmf_to_utc(js_input["end_utc"])
    elif len(sys.argv) >= 4:
        bounds = tuple([float(c) for c in sys.argv[1].split(',')])
        from_utc = dt.datetime.strptime(sys.argv[2],'%Y%m%d%H%M%S').replace(tzinfo=pytz.UTC)
        to_utc = dt.datetime.strptime(sys.argv[3],'%Y%m%d%H%M%S').replace(tzinfo=pytz.UTC)
        if len(sys.argv) == 4:
            sat_sources = ['Terra', 'Aqua', 'SNPP', 'SNPPLL', 'J01LL']
        else:
            sat_sources = [c for c in sys.argv[4].split(',')]
    else:
        print('Usage: ./retrieve_sat.sh input.json')
        print('   or: ./retrieve_sat.sh coord start_time end_time sat_sources')
        print('   notes:')
        print('     *) coord - min_lon,max_lon,min_lat,max_lat')
        print('     *) start_time - string, YYYYMMDDHHMMSS')
        print('     *) end_time - string, YYYYMMDDHHMMSS')
        print('     *) sat_spirces - string, sat1,sta2,sat3')
        print('     *) supported satellite sources - Terra, Aqua, SNPP, SNPPLL, and J01LL')
        sys.exit(-1)

    # configure the basic logger
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    # load configuration JSON
    sys_cfg = load_sys_cfg()

    # create satellite classes
    logging.info('Retrieving all the satellite data in:') 
    logging.info('* Bounding box (%s,%s,%s,%s), and' % bounds)
    logging.info('* Time interval (%s,%s)' % (from_utc, to_utc))
    manifest = {}
    if 'Terra' in sat_sources:
        logging.info('>> MODIS Terra <<')
        terra=Terra(sys_cfg)
        # retrieve granules
        m_terra=terra.retrieve_data_sat(bounds, from_utc, to_utc)
        manifest.update({'Terra': m_terra})
    if 'Aqua' in sat_sources:
        logging.info('>> MODIS Aqua <<')
        aqua=Aqua(sys_cfg)
        # retrieve granules
        m_aqua=aqua.retrieve_data_sat(bounds, from_utc, to_utc)
        manifest.update({'Aqua': m_aqua})
    if 'SNPP' in sat_sources:
        logging.info('>> S-NPP VIIRS <<')
        snpp=SNPP(sys_cfg)
        # retrieve granules
        m_snpp=snpp.retrieve_data_sat(bounds, from_utc, to_utc)
        manifest.update({'SNPP': m_snpp})
    if 'SNPPLL' in sat_sources:
        logging.info('>> S-NPP VIIRS LL <<')
        snppll=SNPPLL(sys_cfg)
        # retrieve granules
        m_snppll=snppll.retrieve_data_sat(bounds, from_utc, to_utc)
        manifest.update({'SNPPLL': m_snppll})
    if 'J01LL' in sat_sources:
        logging.info('>> J01 VIIRS LL <<')
        j01ll=J01LL(sys_cfg)
        # retrieve granules
        m_j01ll=j01ll.retrieve_data_sat(bounds, from_utc, to_utc)
        manifest.update({'J01LL': m_j01ll})
    if 'SNPP_HR' in sat_sources:
        logging.info('>> High resolution S-NPP VIIRS <<')
        snpphr=SNPPHR(sys_cfg)
        # retrieve granules
        m_snpphr=snpphr.retrieve_data_sat(bounds, from_utc, to_utc)
        manifest.update({'SNPPHR': m_snpphr})
    if 'NOAA-20' in sat_sources:
        logging.info('>> NOAA-20 VIIRS <<')
        noaa20=NOAA20(sys_cfg)
        # retrieve granules
        m_noaa20=noaa20.retrieve_data_sat(bounds, from_utc, to_utc)
        manifest.update({'NOAA20': m_noaa20})
