#
# Angel Farguell, CU Denver
#

from ingest.MODIS import Terra, Aqua
from ingest.VIIRS import SNPP, NOAA20
from wrf.wps_domains import WPSDomainConf
from utils import load_sys_cfg, esmf_to_utc, Dict
import logging, sys, json

if __name__ == '__main__':
	if len(sys.argv) != 2:
		print('Usage: %s input.json' % sys.argv[0])
		print('		supported atellite sources: Terra, Aqua, S-NPP')
		sys.exit(-1)

	# load input JSON
	try:
		js_input = Dict(json.load(open(sys.argv[1])))
	except IOError:
		import sys
		logging.critical('Cannot find input json file.')
		sys.exit(2)

	# configure the basic logger
	logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
	# load configuration JSON
	sys_cfg = load_sys_cfg()

	# inputs
	sat_sources = list(js_input["satellite_source"])
	domain=WPSDomainConf(js_input["domains"]).domains[-1]
	latloni=domain.ij_to_latlon(0,0)
	latlonf=domain.ij_to_latlon(domain.domain_size[0],domain.domain_size[1])
	bounds = (latloni[1],latlonf[1],latloni[0],latlonf[0])
	from_utc = esmf_to_utc(js_input["start_utc"])
	to_utc = esmf_to_utc(js_input["end_utc"])

	# create JPSS classes
	logging.info('Retrieving all the satellite data for:')
	if 'Terra' in sat_sources:
		logging.info('> MODIS Terra')
		terra=Terra(sys_cfg)
		# retrieve granules
		m_terra=terra.retrieve_data_jpss(bounds, from_utc, to_utc)
	if 'Aqua' in sat_sources:
		logging.info('> MODIS Aqua')
		aqua=Aqua(sys_cfg)
		# retrieve granules
		m_aqua=aqua.retrieve_data_jpss(bounds, from_utc, to_utc)
	if 'S-NPP' in sat_sources:
		logging.info('> MODIS Aqua')
		snpp=SNPP(sys_cfg)
		# retrieve granules
		m_snpp=snpp.retrieve_data_jpss(bounds, from_utc, to_utc)

