#
# Angel Farguell, CU Denver
#

from __future__ import absolute_import
from ingest.sat_source_aws import SatSourceAWS

class GOES16(SatSourceAWS):
    """
    GOES16 satellite source.
    """

    def __init__(self, arg):
        super(GOES16, self).__init__(arg)

    # instance variables
    id='G16'
    info_url=''
    info='GOES16 Geostationary Satellite'
    prefix='G16'
    platform='goes16'
