#
# Angel Farguell, CU Denver
#

from __future__ import absolute_import
from __future__ import print_function
from ingest.MODIS import Terra, Aqua
from ingest.VIIRS import SNPP, SNPPHR, NOAA20
from wrf.wps_domains import WPSDomainConf
from utils import load_sys_cfg, esmf_to_utc, Dict
import datetime as dt
import logging, sys, json

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
	elif len(sys.argv) == 4:
		bounds = tuple([float(c) for c in sys.argv[1].split(',')])
		st = dt.datetime.strptime(sys.argv[2],'%Y%m%d%H%M%S')
		et = dt.datetime.strptime(sys.argv[3],'%Y%m%d%H%M%S')
		st_esmf = '%d-%02d-%02d_%02d:%02d:%02d' % (st.year,st.month,st.day,st.hour,st.minute,st.second)
		et_esmf = '%d-%02d-%02d_%02d:%02d:%02d' % (st.year,st.month,st.day,st.hour,st.minute,st.second)
		from_utc = esmf_to_utc(st_esmf)
		to_utc = esmf_to_utc(et_esmf)
		sat_sources = ['Terra', 'Aqua', 'SNPP']
	else:
		print('Usage: ./retrieve_sat.sh input.json')
		print('   or: ./retrieve_sat.sh coord start_time end_time')
		print('	  notes:')
		print('	  	*) coord - min_lon,max_lon,min_lat,max_lat')
		print('	  	*) start_time - string, YYYYMMDDHHMMSS')
		print('	  	*) end_time - string, YYYYMMDDHHMMSS')
		print('		*) supported satellite sources - Terra, Aqua, and SNPP')
		sys.exit(-1)

	# configure the basic logger
	logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
	# load configuration JSON
	sys_cfg = load_sys_cfg()

	# create satellite classes
	logging.info('Retrieving all the satellite data in:') 
	logging.info('* Bounding box (%s,%s,%s,%s), and' % bounds)
	logging.info('* Time interval (%s,%s)' % (from_utc, to_utc))
	if 'Terra' in sat_sources:
		logging.info('>> MODIS Terra <<')
		terra=Terra(sys_cfg)
		# retrieve granules
		m_terra=terra.retrieve_data_sat(bounds, from_utc, to_utc)
	if 'Aqua' in sat_sources:
		logging.info('>> MODIS Aqua <<')
		aqua=Aqua(sys_cfg)
		# retrieve granules
		m_aqua=aqua.retrieve_data_sat(bounds, from_utc, to_utc)
	if 'SNPP' in sat_sources:
		logging.info('>> S-NPP VIIRS <<')
		snpp=SNPP(sys_cfg)
		# retrieve granules
		m_snpp=snpp.retrieve_data_sat(bounds, from_utc, to_utc)
	if 'SNPP_HR' in sat_sources:
		logging.info('>> High resolution S-NPP VIIRS <<')
		snpphr=SNPPHR(sys_cfg)
		# retrieve granules
		m_snpphr=snpphr.retrieve_data_sat(bounds, from_utc, to_utc)
		print(m_snpphr)
	if 'NOAA-20' in sat_sources:
		logging.info('>> NOAA-20 VIIRS <<')
		noaa20=NOAA20(sys_cfg)
		# retrieve granules
		m_noaa20=noaa20.retrieve_data_sat(bounds, from_utc, to_utc)

