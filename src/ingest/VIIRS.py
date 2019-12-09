#
# Angel Farguell, CU Denver
#

from ingest.sat_source import SatSource

class VIIRS(SatSource):
	"""
	VIIRS (Visible Infrared Imaging Radiometer Suite) satellite source.
	"""

	def __init__(self, arg):
		super(VIIRS, self).__init__(arg)

	# instance variables
	url='https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5000'

class SNPP(VIIRS):
	"""
	S-NPP VIIRS (Visible Infrared Imaging Radiometer Suite) satellite source.
	"""

	def __init__(self, arg):
		super(SNPP, self).__init__(arg)

	# instance variables
	id='SNPP'
	info_url='https://www.nasa.gov/mission_pages/NPP/mission_overview/index.html'
	info='S-NPP Visible Infrared Imaging Radiometer Suite (VIIRS)'
	prefix='VNP'
	geo_prefix='VNP03MODLL'
	fire_prefix='VNP14'
	ref_prefix='VNP09'
	platform='S-NPP'

class SNPPHR(VIIRS):
	"""
	High Resolution S-NPP VIIRS (Visible Infrared Imaging Radiometer Suite) satellite source.
	"""

	def __init__(self, arg):
		super(SNPPHR, self).__init__(arg)

	# instance variables
	id='SNPPHR'
	info_url='https://www.nasa.gov/mission_pages/NPP/mission_overview/index.html'
	info='High Resolution S-NPP Visible Infrared Imaging Radiometer Suite (VIIRS)'
	prefix='VNPHR'
	pre_geo_prefix='VNP03MODLL'
	geo_prefix='NPP_IMFTS_L1'
	fire_prefix='VNP14IMG'
	ref_prefix='VNP09'
	platform='S-NPP'

class NOAA20(VIIRS):
	"""
	NOAA-20 VIIRS (Visible Infrared Imaging Radiometer Suite) satellite source.
	"""

	def __init__(self, arg):
		super(NOAA20, self).__init__(arg)

	# instance variables
	id='NOAA20'
	info_url=''
	info='NOAA-20 Visible Infrared Imaging Radiometer Suite (VIIRS)'
	prefix=''
	geo_prefix=''
	fire_prefix=''
	ref_prefix=''
	platform='NOAA-20'

