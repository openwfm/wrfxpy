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

from __future__ import absolute_import
from utils import ensure_dir, make_dir, symlink_unless_exists

import os
import os.path as osp
from six.moves import map


class WRFCloner(object):
    """
    This class is responsible for cloning a (working) WPS/WRF installation
    into a new location that can be used for running jobs without interfering
    with the original installation.
    """

    def __init__(self, args):
        """
        Initialize with WRF/WPS installation directories from which clones will be built.

        #param args: a dictionary of arguments containing the following keys
                sys_install_path: installation directory of the current system
                wrf_install_path: WRF installation directory
                wps_install_path: WPS installation directory
                wrf_serial_install_path: WRF installation directory with real.exe which runs without MPI
        :return:
        """
        self.sys_idir = args['sys_install_path']
        self.wrf_idir = args['wrf_install_path']
        self.wps_idir = args['wps_install_path']
        self.wrf_serial_idir = args.get('wrf_serial_install_path')

    def clone_wps(self, tgt, with_files):
        """
        Clone the WPS installation directory (self.wps_idir) together with the chosen table files vtables
        and additional files with_files.  The WPS clone is created in directory tgt.

        :param tgt: target directory into which WPS is cloned
        :param with_files: a list of files from the WPS source directory that should be symlinked
        :return:
        """
        src = self.wps_idir

        # build a list of all files that are simply symlinked
        symlinks = list(self.wps_exec_files)
        symlinks.extend(with_files)

        # create target directory (and all intermediate subdirs if necessary)
        make_dir(tgt)

        # clone all WPS executables
        list(map(lambda x: symlink_unless_exists(osp.join(src, x), osp.join(tgt, x)), symlinks))

    def clone_wrf(self, tgt, with_files):
        """
        Clone the WRFV3 directory (self.wrf_idir) into tgt together with the additional files with_files.

        :param tgt: target directory into which WRF is cloned
        :param with_files: a list of files from the WPS source directory that should be symlinked
        :return:
        """
        src = osp.join(self.wrf_idir, "run")

        # gather all files to symlink in one place
        symlinks = list(self.wrf_files)
        symlinks.extend(with_files)

        if self.wrf_serial_idir is None:
            symlinks.extend(wrf_serial_files)

        # create target directory (and all intermediate subdirs if necessary)
        make_dir(tgt)

        # symlink all at once
        list(map(lambda x: symlink_unless_exists(osp.join(src, x), osp.join(tgt, x)), symlinks))

        if self.wrf_serial_idir is not None:
            # gather all files to symlink in one place
            symlinks = list(self.wrf_serial_files)
            src = osp.join(self.wrf_serial_idir, "run")
            # symlink all at once
            list(map(lambda x: symlink_unless_exists(osp.join(src, x), osp.join(tgt, x)), symlinks))

    # list of executable file that must be symlinked in WPS directory
    wps_exec_files = ['geogrid.exe', 'metgrid.exe', 'ungrib.exe']

    # list of files that must(?) be symlinked from the WRFV3 directory
    wrf_files = ["CAM_ABS_DATA", "CAM_AEROPT_DATA", "co2_trans", "ETAMPNEW_DATA", "ETAMPNEW_DATA_DBL",
                 "ETAMPNEW_DATA.expanded_rain", "ETAMPNEW_DATA.expanded_rain_DBL", "GENPARM.TBL",
                 "gribmap.txt", "grib2map.tbl", "LANDUSE.TBL", "MPTABLE.TBL",
                 "ozone.formatted", "ozone_lat.formatted", "ozone_plev.formatted",
                 "RRTM_DATA", "RRTM_DATA_DBL", "RRTMG_LW_DATA", "RRTMG_LW_DATA_DBL",
                 "RRTMG_SW_DATA", "RRTMG_SW_DATA_DBL", "SOILPARM.TBL", "tc.exe", "tr49t67", "tr49t85",
                 "tr67t85", "URBPARM.TBL", "URBPARM_UZE.TBL", "VEGPARM.TBL", "wrf.exe"]

    # list of files that must(?) be symlinked from the WRFV3 serial directory
    wrf_serial_files = ["real.exe"]

def test_cloner():
    inst_dir = '/share_home/mvejmelka/Packages/wrf-fire.openwfm.clamping2/'
    w = WRFCloner(osp.join(inst_dir, 'WRFV3'), osp.join(inst_dir, 'WPS'))
    import shutil
    if os.path.exists('test_wps'): shutil.rmtree('test_wps')
    if os.path.exists('test_wrf'): shutil.rmtree('test_wrf')
    w.clone_wps('test_wps', {'geogrid_vtable': 'geogrid/GEOGRID.TBL'}, ['namelist.wps'])
    w.clone_wrf('test_wrf', [])
    return w
