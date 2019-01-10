from ingest.MODIS import Terra, Aqua
from ingest.VIIRS import SNPP
from utils import load_sys_cfg
import logging

if __name__ == '__main__':
	# configure the basic logger
	logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
	# load configuration JSON
	sys_cfg = load_sys_cfg()
	# create JPSS classes
	terra=Terra(sys_cfg)
	aqua=Aqua(sys_cfg)
	snpp=SNPP(sys_cfg)
	# example coordinates and time
	bounds=[-122.04251098632812, -120.97007751464844, 39.3486213684082, 40.169677734375]
	time=('2018-11-08T19:55:00Z', '2018-11-08T21:30:02Z')
	# search granules
	m_terra=terra.retrieve_data_jpss(bounds,time)
	m_aqua=aqua.retrieve_data_jpss(bounds,time)
	m_snpp=snpp.retrieve_data_jpss(bounds,time)
