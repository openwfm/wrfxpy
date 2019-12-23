#
# Angel Farguell, CU Denver
#

from ingest.sat_source import SatSource

class MODIS(SatSource):
	"""
	MODIS (Moderate Resolution Imaging Spectroradiometer) satellite source.
	"""

	def __init__(self, arg):
		super(MODIS, self).__init__(arg)

	# instance variables
	url='https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6'

class Terra(MODIS):
	"""
	Terra MODIS (Moderate Resolution Imaging Spectroradiometer) satellite source.
	"""

	def __init__(self, arg):
		super(Terra, self).__init__(arg)

	# instance variables
	id='Terra'
	info_url='https://terra.nasa.gov/about/terra-instruments/modis'
	info='Terra Moderate Resolution Imaging Spectroradiometer (MODIS)'
	prefix='MOD'
	geo_prefix='MOD03'
	fire_prefix='MOD14'
	ref_prefix='MOD09'
	platform='Terra'

class Aqua(MODIS):
	"""
	Aqua MODIS (Moderate Resolution Imaging Spectroradiometer) satellite source.
	"""

	def __init__(self, arg):
		super(Aqua, self).__init__(arg)

	# instance variables
	id='Aqua'
	info_url='https://aqua.nasa.gov/modis'
	info='Aqua Moderate Resolution Imaging Spectroradiometer (MODIS)'
	prefix='MYD'
	geo_prefix='MYD03'
	fire_prefix='MYD14'
	ref_prefix='MYD09'
	platform='Aqua'

