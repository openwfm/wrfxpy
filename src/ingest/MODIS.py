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
    fire_nrt_col='6'
    
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
    geo_nrt_prefix='MOD03'
    fire_nrt_prefix='MOD14'
    geo_collection_id='C1379767668-LAADS'
    fire_collection_id='C193529945-LPDAAC_ECS'
    geo_nrt_collection_id='C1426422512-LANCEMODIS'
    fire_nrt_collection_id='C1219250604-LANCEMODIS'
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
    geo_nrt_prefix='MYD03'
    fire_nrt_prefix='MYD14'
    geo_collection_id='C1379841358-LAADS'
    fire_collection_id='C193529462-LPDAAC_ECS'
    geo_nrt_collection_id='C1426640814-LANCEMODIS'
    fire_nrt_collection_id='C1219248602-LANCEMODIS'
    platform='Aqua'
