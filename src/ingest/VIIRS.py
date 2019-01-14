#
# Angel Farguell, CU Denver
#

from ingest.jpss_source import JPSSSource, JPSSError

class VIIRS(JPSSSource):
	"""
	VIIRS (Visible Infrared Imaging Radiometer Suite) JPSS source.
	"""

	def __init__(self, arg):
		super(VIIRS, self).__init__(arg)

	# instance variables
	url='https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5000'

class SNPP(VIIRS):
	"""
	S-NPP VIIRS (Visible Infrared Imaging Radiometer Suite) JPSS source.
	"""

	def __init__(self, arg):
		super(SNPP, self).__init__(arg)

	# instance variables
	info_url=''
	info='S-NPP Visible Infrared Imaging Radiometer Suite (VIIRS)'
	prefix='VNP'
	geo_prefix='VNP03MODLL'
	fire_prefix='VNP14'
	ref_prefix='VNP09'
	platform=''

class NOAA20(VIIRS):
	"""
	NOAA-20 VIIRS (Visible Infrared Imaging Radiometer Suite) JPSS source.
	"""

	def __init__(self, arg):
		super(NOAA20, self).__init__(arg)

	# instance variables
	info_url=''
	info='NOAA-20 Visible Infrared Imaging Radiometer Suite (VIIRS)'
	prefix=''
	geo_prefix=''
	fire_prefix=''
	ref_prefix=''
	platform=''

