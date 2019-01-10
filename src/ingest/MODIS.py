#
# Angel Farguell, CU Denver
#

from ingest.jpss_source import JPSSSource, JPSSError

class MODIS(JPSSSource):
	"""
	MODIS (Moderate Resolution Imaging Spectroradiometer) JPSS source.
	"""

	def __init__(self, arg):
		super(MODIS, self).__init__(arg)

	# instance variables
	url='https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6'
	version='6'

class Terra(MODIS):
	"""
	Terra MODIS (Moderate Resolution Imaging Spectroradiometer) JPSS source.
	"""

	def __init__(self, arg):
		super(Terra, self).__init__(arg)

	# instance variables
	info_url='https://terra.nasa.gov/about/terra-instruments/modis'
	info='Terra Moderate Resolution Imaging Spectroradiometer (MODIS)'
	prefix='MOD'
	platform='Terra'

class Aqua(MODIS):
	"""
	Aqua MODIS (Moderate Resolution Imaging Spectroradiometer) JPSS source.
	"""

	def __init__(self, arg):
		super(Aqua, self).__init__(arg)

	# instance variables
	info_url='https://aqua.nasa.gov/modis'
	info='Aqua Moderate Resolution Imaging Spectroradiometer (MODIS)'
	prefix='MYD'
	platform='Aqua'

