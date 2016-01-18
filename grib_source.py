# Copyright (C) 2013-2016 Martin Vejmelka, UC Denver
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

from utils import ensure_dir, symlink_unless_exists
from datetime import datetime, timedelta
import requests
import os
import os.path as osp


class GribError(Exception):
    """
    Raised when a GribSource cannot retrieve GRIBs.
    """
    pass


class GribSource(object):
    """
    Parent class of all grib sources.
    """

    def __init__(self, ingest_dir):
        """
        Initialize grib source with ingest directory (where GRIB files are stored).

        :param ingest_dir: root of GRIB storage
        """
        self.ingest_dir = osp.abspath(ingest_dir)

    def vtables(self):
        """
        Returns the vtables that must be set for use with this source.
        """
        return {}

    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with this GRIB source.
        """
        return {}

    def retrieve_gribs(self, from_utc, to_utc):
        """
        Attempts to retrieve the GRIB files for the forecast range.

        :param from_utc: forecast start time
        :param to_utc: forecast end time
        :return: a list of paths to local GRIB files
        """
        pass

    def download_grib(self, url_base, rel_path):
        """
        Download a GRIB file from a GRIB service and stream to rel_path in ingest_dir.

        :param url_base: the base URL part of the GRIB service
        :param rel_path: the relative path of the file (w.r.t GRIB base url and w.r.t self.ingest_dir)
        """
        r = requests.get(url_base + '/' + rel_path, stream=True)
        grib_path = osp.join(self.ingest_dir, rel_path)
        with open(ensure_dir(grib_path), 'wb') as f:
            for chunk in r.iter_content(1024 * 1024):
                f.write(chunk)

    def symlink_gribs(self, manifest, wps_dir):
        """
        Make symlinks in the form GRIBFILE.XYZ to all manifest files into wps_dir.

        :param manifest: relative paths (w.r.t. ingest_dir) to GRIB files we want linked
        :param wps_dir: the WPS directory where we want the symlinks to appear
        :return:
        """
        for rel_path, grib_name in zip(manifest, generate_grib_names()):
            symlink_unless_exists(osp.join(self.ingest_dir, rel_path), osp.join(wps_dir, grib_name))


class HRRR(GribSource):
    """
    The HRRR grib source as provided by NOMADS.
    """

    def __init__(self, ingest_dir):
        super(HRRR, self).__init__(ingest_dir)

    def vtables(self):
        """
        Returns the variable tables that must be linked in for use with the HRRR data source.
        :return:
        """
        return {'geogrid_vtable': 'GEOGRID.TBL.HRRR',
                'ungrib_vtable': 'Vtable.HRRR',
                'metgrid_vtable': 'METGRID.TBL.HRRR'}

    def namelist_keys(self):
        """
        Returns the namelist keys that must be modified in namelist.input with HRRR.

        HRRR requires that ''num_metgrid_soil_levels'' is set to 8.
        """
        #return { 'domains' : { 'num_metgrid_soil_levels': 8 }}
        return {}

    def retrieve_gribs(self, from_utc, to_utc):
        """
        Attempts to retrieve the files to satisfy the simulation request from_utc - to_utc.

        Starts with the most recent cycle available an hour ago, then moves further
        into the past.  For each candidate cycle, the filenames are computed, the local cache is
        checked for files that are already there.  The presence of remaining files is checked
        on server, if not available, we try an older cycle, if yes, download is attempted.
        Once all files are downloaded, the manifest is returned.

        :param from_utc: forecast start time
        :param to_utc: forecast end time
        :return: a list of paths to local GRIB files
        """
        # ensure minutes and seconds are zero, simplifies arithmetic later
        from_utc = from_utc.replace(minute=0, second=0)
        to_utc = to_utc.replace(minute=0, second=0)

        # select cycle (at least one hour behind)
        cycle_start = from_utc - timedelta(hours=1)

        # check if the request is even satisfiable
        delta = to_utc - cycle_start
        fc_hours = delta.days * 24 + delta.seconds / 3600

        if fc_hours > 15:
            raise GribError('Unsatisfiable: HRRR only forecasts 15 hours ahead.')

        # computes the relative paths of the desired files (the manifest)
        manifest = self.compute_manifest(from_utc, fc_hours)

        # check what's available locally
        nonlocals = filter(lambda x: not osp.exists(osp.join(self.ingest_dir, x)), manifest)

        # check if GRIBs we don't are available remotely
        url_base = self.remote_url
        unavailables = filter(lambda x: requests.head(url_base + '/' + x).status_code != 200, nonlocals)
        if len(unavailables) > 0:
            raise GribError('Unsatisfiable: GRIBs %s not available.' % repr(unavailables))

        # download all gribs we need
        map(lambda x: self.download_grib(url_base, x), nonlocals)

        # return manifest
        return manifest

    def compute_manifest(self, cycle_start, fc_hours):
        """
        Computes the relative paths of required GRIB files.

        Note, the system is built so that relative paths are the same in local cache
        and in remote system w.r.t. URL base.

        :param cycle_start: UTC time of cycle start
        :param fc_hours: final forecast hour 
        """
        year, mon, day, hour = cycle_start.year, cycle_start.month, cycle_start.day, cycle_start.hour
        return map(lambda x: self.path_tmpl % (year, mon, day, hour, x), range(fc_hours))

    # instance variables
    remote_url = 'http://www.ftp.ncep.noaa.gov/data/nccf/nonoperational/com/hrrr/prod'
    path_tmpl = 'hrrr.%04d%02d%02d/hrrr.t%02dz.wrfprsf%02d.grib2'


def generate_grib_names():
    """
    Keeps generating gribfile names from GRIBFILE.AAA to ZZZ.
    """
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    for c1 in alphabet:
        for c2 in alphabet:
            for c3 in alphabet:
                yield "GRIBFILE." + c1 + c2 + c3
