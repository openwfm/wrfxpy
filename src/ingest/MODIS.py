#
# Angel Farguell, CU Denver
#

from __future__ import absolute_import
from ingest.sat_source_cmr import SatSourceCMR

class MODIS(SatSourceCMR):
    """
    MODIS (Moderate Resolution Imaging Spectroradiometer) satellite source.
    """

    def __init__(self, arg):
        super(MODIS, self).__init__(arg)

    # instance variables
    geo_col='61'
    fire_col='61'
    geo_nrt_col='61'
    fire_nrt_col='61'
    
class Terra(MODIS):
    """
    Terra MODIS (Moderate Resolution Imaging Spectroradiometer) satellite source.
    """

    def __init__(self, arg):
        super(Terra, self).__init__(arg)

    # instance variables
    id='Terra'
    info='Terra Moderate Resolution Imaging Spectroradiometer (MODIS)'
    info_url='https://www.earthdata.nasa.gov/data/catalog/lpcloud-mod14-061'
    info_url_nrt='https://www.earthdata.nasa.gov/data/catalog/lancemodis-mod14-6.1nrt'
    prefix='MOD'
    geo_prefix='MOD03'
    fire_prefix='MOD14'
    geo_nrt_prefix='MOD03'
    fire_nrt_prefix='MOD14'
    geo_collection_id='C1379767668-LAADS'
    fire_collection_id='C2271754179-LPCLOUD'
    geo_nrt_collection_id='C1426422512-LANCEMODIS'
    fire_nrt_collection_id='C2007630683-LANCEMODIS'
    platform='Terra'

class Aqua(MODIS):
    """
    Aqua MODIS (Moderate Resolution Imaging Spectroradiometer) satellite source.
    """

    def __init__(self, arg):
        super(Aqua, self).__init__(arg)

    # instance variables
    id='Aqua'
    info='Aqua Moderate Resolution Imaging Spectroradiometer (MODIS)'
    info_url='https://www.earthdata.nasa.gov/data/catalog/lpcloud-myd14-061'
    info_url_nrt='https://www.earthdata.nasa.gov/data/catalog/lancemodis-myd14-6.1nrt'
    prefix='MYD'
    geo_prefix='MYD03'
    fire_prefix='MYD14'
    geo_nrt_prefix='MYD03'
    fire_nrt_prefix='MYD14'
    geo_collection_id='C1379841358-LAADS'
    fire_collection_id='C2278858993-LPCLOUD'
    geo_nrt_collection_id='C1426640814-LANCEMODIS'
    fire_nrt_collection_id='C2007630791-LANCEMODIS'
    platform='Aqua'
