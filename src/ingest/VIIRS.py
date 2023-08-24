#
# Angel Farguell, CU Denver
#

from __future__ import absolute_import
from ingest.sat_source_cmr import SatSourceCMR

class VIIRS(SatSourceCMR):
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
    geo_prefix='VNP03MOD'
    fire_prefix='VNP14'
    geo_nrt_prefix='VNP03MOD_NRT'
    fire_nrt_prefix='VNP14_NRT'
    geo_collection_id='C1344456864-LAADS'
    fire_collection_id='C1392010612-LPDAAC_ECS'
    geo_nrt_collection_id='C1412085464-LANCEMODIS'
    fire_nrt_collection_id='C1538921885-LANCEMODIS'
    platform='S-NPP'
    geo_col='5110'
    fire_col='5000'
    geo_nrt_col='5001'
    fire_nrt_col='5000'

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
    geo_prefix='VNP03IMG'
    fire_prefix='VNP14IMG'
    geo_nrt_prefix='VNP03IMG_NRT'
    fire_nrt_prefix='VNP14IMG_NRT'
    geo_collection_id='C1344460841-LAADS'
    fire_collection_id=''
    geo_nrt_collection_id='C1412839136-LANCEMODIS'
    fire_nrt_collection_id='C1344295580-LANCEMODIS'
    platform='S-NPP'
    geo_col='5110'
    fire_col='5000'
    geo_nrt_col='5001'
    fire_nrt_col='5001'

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
    prefix='VJ1'
    geo_prefix='VJ103MOD'
    fire_prefix='VJ114'
    geo_nrt_prefix='VJ103MOD_NRT'
    fire_nrt_prefix='VJ114_NRT'
    geo_collection_id='C1624979410-LAADS'
    fire_collection_id=''
    geo_nrt_collection_id='C1604697279-LANCEMODIS'
    fire_nrt_collection_id=''
    platform='NOAA-20'

class NOAA20HR(VIIRS):
    """
    NOAA-20 VIIRS (Visible Infrared Imaging Radiometer Suite) satellite source.
    """

    def __init__(self, arg):
        super(NOAA20, self).__init__(arg)

    # instance variables
    id='NOAA20HR'
    info_url=''
    info='NOAA-20 Visible Infrared Imaging Radiometer Suite (VIIRS)'
    prefix='VJ1HR'
    geo_prefix='VJ103IMG'
    fire_prefix='VJ114IMG'
    geo_nrt_prefix='VJ103IMG_NRT'
    fire_nrt_prefix='VJ114IMG_NRT'
    geo_collection_id='C1626204226-LAADS'
    fire_collection_id=''
    geo_nrt_collection_id='C1604644159-LANCEMODIS'
    fire_nrt_collection_id=''
    platform='NOAA-20'
