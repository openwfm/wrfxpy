#
# Angel Farguell, CU Denver
#

from __future__ import absolute_import
from ingest.sat_source_archive import SatSourceArchive

class VIIRSLL(SatSourceArchive):
    """
    VIIRS (Visible Infrared Imaging Radiometer Suite) satellite source.
    """

    def __init__(self, arg):
        super(VIIRSLL, self).__init__(arg)


class SNPPLL(VIIRSLL):
    """
    S-NPP VIIRS (Visible Infrared Imaging Radiometer Suite) Low Latency  satellite source.
    """

    def __init__(self, arg):
        super(SNPPLL, self).__init__(arg)

    # instance variables
    base_url='http://52.200.226.137/archive/allData/5001'
    id='SNPPLL'
    info_url='https://www.nasa.gov/mission_pages/NPP/mission_overview/index.html'
    info='S-NPP Visible Infrared Imaging Radiometer Suite (VIIRS) Low Latency'
    prefix='VNPLL'
    geo_nrt_prefix='VNP03IMG_NRT'
    fire_nrt_prefix='VNP14IMG_NRT'
    platform='S-NPP'

class VJ01LL(VIIRSLL):
    """
    J01 VIIRS (Visible Infrared Imaging Radiometer Suite) Low Latency satellite source.
    """

    def __init__(self, arg):
        super(VJ01LL, self).__init__(arg)

    # instance variables
    base_url='http://52.200.226.137/archive/allData/5200'
    id='VJ01LL'
    info_url='https://www.nasa.gov/mission_pages/NPP/mission_overview/index.html'
    info='J01 Visible Infrared Imaging Radiometer Suite (VIIRS) Low Latency'
    prefix='VJ1LL'
    geo_nrt_prefix='VJ103IMG_NRT'
    fire_nrt_prefix='VJ114IMG_NRT'
    platform='J01'

