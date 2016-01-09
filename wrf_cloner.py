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

import os
import os.path as pth

class WRFCloner:
  """
  This class is responsible for cloning a (working) WPS/WRF installation
  into a new location that can be used for running jobs without interfering
  with the original installation.
  """

  def __init__(self, wrf_install_dir, wps_install_dir):
    """
    Initialize the cloner with WRF/WPS installation directories from which
    clones will be built.
    """
    self.wrf_idir = wrf_install_dir
    self.wps_idir = wps_install_dir

 
  def clone_wps(self, tgt, vtables, simple_links):
    """
    Clone the WPS installation directory (self.wps_idir) together with the chosen table files vtables
    and additional files with_files.  The WPS clone is created in directory tgt.

    Args:
      tgt: target directory
      vtables: a dictionary with keys from list ['geogrid_vtable', 'ungrib_vtable', 'metgrid_vtable']
      simple_links: a list of files from the WPS source directory that should be symlinked
    """
    src = self.wps_idir
    vtable_locs = self.vtable_locations

    # create target directory (and all intermediate subdirs if necessary)
    os.makedirs(tgt)

    # clone all WPS executables
    map(lambda x: os.symlink(pth.join(src, x), pth.join(tgt, x)), self.wps_exec_files)

    # clone all vtables (build symlink name, ensure directories exist, create the symlink)
    for vt, rel_path in vtables.iteritems():
      tgt_loc = pth.join(tgt, vtable_locs[vt])
      if not pth.exists(tgt_loc):
        os.makedirs(pth.dirname(tgt_loc))
        os.symlink(pth.join(src, rel_path), tgt_loc)

    # symlink in all the simple links
    map(lambda x: os.symlink(pth.join(src, x), pth.join(tgt, x)), simple_links)


  # list of executable file that must be symlinked in WPS directory
  wps_exec_files = [ 'geogrid.exe', 'metgrid.exe', 'ungrib.exe' ]

  # where are the symlink locations for vtable files (name of symlink)
  vtable_locations = { 'geogrid_vtable' : 'geogrid/GEOGRID.TBL',
                       'ungrib_vtable' : 'Vtable',
                       'metgrid_vtable' : 'metgrid/METGRID.TBL' }



def test_cloner():
  inst_dir = '/share_home/mvejmelka/Packages/wrf-fire.openwfm.clamping2/'
  w = WRFCloner(pth.join(inst_dir, 'WRFV3'), pth.join(inst_dir, 'WPS'))
  w.clone_wps('test_wps', { 'geogrid_vtable' : 'geogrid/GEOGRID.TBL.HRRR' }, [ 'namelist.wps' ])
  return w


