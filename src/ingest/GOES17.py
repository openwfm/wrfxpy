#
# Angel Farguell, CU Denver
#

from __future__ import absolute_import
from ingest.sat_source_aws import SatSourceAWS

class GOES17(SatSourceAWS):
    """
    GOES16 satellite source.
    """

    def __init__(self, arg):
        super(GOES16, self).__init__(arg)

    # instance variables
    id='G17'
    info_url=''
    info='GOES17 Geostationary Satellite'
    prefix='G17'
    platform='goes17'
