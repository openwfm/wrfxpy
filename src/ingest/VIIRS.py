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
    url='https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5200'

class SNPP(VIIRS):
    """
    S-NPP VIIRS (Visible Infrared Imaging Radiometer Suite) satellite source.
    """

    def __init__(self, arg):
        super(SNPP, self).__init__(arg)

    # instance variables
    id='SNPP'
    info='S-NPP Visible Infrared Imaging Radiometer Suite (VIIRS)'
    info_url='https://www.earthdata.nasa.gov/data/catalog/lpcloud-vnp14-002'
    info_url_nrt='https://www.earthdata.nasa.gov/data/catalog/lancemodis-vnp14-nrt-2'
    prefix='VNP'
    geo_prefix='VNP03MOD'
    fire_prefix='VNP14'
    geo_nrt_prefix='VNP03MOD_NRT'
    fire_nrt_prefix='VNP14_NRT'
    geo_collection_id='C2105092427-LAADS'
    fire_collection_id='C2545314536-LPCLOUD'
    geo_nrt_collection_id='C2185511251-LANCEMODIS'
    fire_nrt_collection_id='C2888489803-LANCEMODIS'
    platform='S-NPP'
    geo_col='5200'
    fire_col='5200'
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
    info='High Resolution S-NPP Visible Infrared Imaging Radiometer Suite (VIIRS)'
    info_url='https://www.earthdata.nasa.gov/data/catalog/lpcloud-vnp14img-002'
    info_url_nrt='https://www.earthdata.nasa.gov/data/catalog/lancemodis-vnp14img-nrt-2'
    prefix='VNPHR'
    geo_prefix='VNP03IMG'
    fire_prefix='VNP14IMG'
    geo_nrt_prefix='VNP03IMG_NRT'
    fire_nrt_prefix='VNP14IMG_NRT'
    geo_collection_id='C2105092163-LAADS'
    fire_collection_id='C2734202914-LPCLOUD'
    geo_nrt_collection_id='C2185522599-LANCEMODIS'
    fire_nrt_collection_id='C1886251885-LANCEMODIS'
    platform='S-NPP'
    geo_col='5200'
    fire_col='5200'
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
    info='NOAA-20 Visible Infrared Imaging Radiometer Suite (VIIRS)'
    info_url='https://www.earthdata.nasa.gov/data/catalog/lpcloud-vj114-002'
    info_url_nrt='https://www.earthdata.nasa.gov/data/catalog/lancemodis-vj114-nrt-2'
    prefix='VJ1'
    geo_prefix='VJ103MOD'
    fire_prefix='VJ114'
    geo_nrt_prefix='VJ103MOD_NRT'
    fire_nrt_prefix='VJ114_NRT'
    geo_collection_id='C2105084593-LAADS'
    fire_collection_id='C2545310869-LPCLOUD'
    geo_nrt_collection_id='C2208781576-LANCEMODIS'
    fire_nrt_collection_id='C2888590350-LANCEMODIS'
    platform='NOAA-20'
    geo_col='5201'
    fire_col='5200'
    geo_nrt_col=''
    fire_nrt_col=''

class NOAA20HR(VIIRS):
    """
    NOAA-20 VIIRS (Visible Infrared Imaging Radiometer Suite) satellite source.
    """

    def __init__(self, arg):
        super(NOAA20HR, self).__init__(arg)

    # instance variables
    id='NOAA20HR'
    info='NOAA-20 Visible Infrared Imaging Radiometer Suite (VIIRS)'
    info_url='https://www.earthdata.nasa.gov/data/catalog/lpcloud-vj114img-002'
    info_url_nrt='https://www.earthdata.nasa.gov/data/catalog/lancemodis-vj114img-nrt-2'
    prefix='VJ1HR'
    geo_prefix='VJ103IMG'
    fire_prefix='VJ114IMG'
    geo_nrt_prefix='VJ103IMG_NRT'
    fire_nrt_prefix='VJ114IMG_NRT'
    geo_collection_id='C2105086226-LAADS'
    fire_collection_id='C2734197957-LPCLOUD'
    geo_nrt_collection_id='C2208793489-LANCEMODIS'
    fire_nrt_collection_id='C1907902788-LANCEMODIS'
    platform='NOAA-20'
    geo_col='5201'
    fire_col='5200'
    geo_nrt_col=''
    fire_nrt_col=''
