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

import netCDF4

class WPSDomainLCC:
    """
    A top-level Lambert Conic Conformal projection domain.
    """
    def __init__(self, dom_id, cfg):
        """
        Initialize a top-level LCC domain.

        :param dom_id: the domain id
        :param cfg: the configuration dictionary
        """
        self.dom_id = dom_id
        if 'precomputed' in cfg:
            self.init_from_precomputed(cfg['precomputed'])
        else:
            self.init_from_dict(cfg)


    def init_from_dict(self, cfg):
        """
        Initialize the domain from a dictionary containing the parameters.

        :param cfg: the parameter dictionary
        """
        self.sr_x = cfg.get('sr_x', 1)
        self.sr_y = cfg.get('sr_y', 1)
        self.geog_res = str(cfg.get('geog_res', '30s'))
        self.parent_id = cfg.get('parent_id', self.dom_id)
        if self.parent_id == self.dom_id:
            # this is a top level domain
            self.sx, self.sy = 1, 1
            self.parent_grid_ratio = 1
            self.dx, self.dy = cfg['dx'], cfg['dy']
            self.ex, self.ey = cfg['nx'], cfg['ny']
            self.ref_lat = cfg['center_lat']
            self.ref_lon = cfg['center_lon']
            self.time_step = cfg.get('time_step', 6*self.dx/1000)
            if 'truelats' in cfg:
                self.truelats = cfg['truelats']
            else:
                self.truelats = [self.ref_lat]*2
            if 'stand_lon' in cfg:
                self.stand_lon = cfg['stand_lon']
            else:
                self.stand_lon = self.ref_lon
        else:
            # this is a child domain, placed normally
            self.parent_grid_ratio = cfg['parent_grid_ratio']
            self.parent_time_step_ratio = cfg['parent_time_ratio']
            self.sx, self.sy = cfg['sx'], cfg['sy']
            self.ex, self.ey = cfg['ex'], cfg['ey']

        # either way, it's not precomputed
        self.precomputed = False


    def init_from_precomputed(self, path):
        """
        Initialize the domain from an existing geo_em.dXX.nc file.

        :param path: path to the geo_em file
        """
        self.precomputed = True
        d = netCDF4.Dataset(path)
        self.sr_x = int(d.getncattr('sr_x'))
        self.sr_y = int(d.getncattr('sr_y'))
        self.geog_res = '30s' # not in geo_em file and does not matter (no processing will be done)
        self.dom_id = int(d.getncattr('grid_id'))
        self.parent_id = int(d.getncattr('parent_id'))
        self.parent_grid_ratio = int(d.getncattr('parent_grid_ratio'))
        self.sx = int(d.getncattr('i_parent_start'))
        self.sy = int(d.getncattr('j_parent_start'))
        self.dx = float(d.getncattr('DX'))
        self.dy = float(d.getncattr('DY'))
        self.ref_lat = float(d.getncattr('CEN_LAT'))
        self.ref_lon = float(d.getncattr('CEN_LON'))
        self.nx = int(d.getncattr('i_parent_end'))
        self.ny = int(d.getncattr('j_parent_end'))
        self.truelats = [float(d.getncattr('TRUELAT1')), float(d.getncattr('TRUELAT2'))] 
        self.stand_lon = float(d.getncattr('STAND_LON'))
        d.close()

    
    def update_wpsnl(self, nml):
        """
        Update the share and geogrid section of the WPS namelist.

        :param nml_geogrid: a dictionary containing the geogrid section of the WPS namelist
        """
        nml_share = nml['share']
        self.update_entry(nml_share, 'subgrid_ratio_x', self.sr_x)
        self.update_entry(nml_share, 'subgrid_ratio_y', self.sr_y)

        # prevent geogrid from re-processing the grid
        self.update_entry(nml_share, 'active_grid', not self.precomputed)

        nml_geogrid = nml['geogrid']
        self.update_entry(nml_geogrid, 'parent_id', self.parent_id)
        self.update_entry(nml_geogrid, 'parent_grid_ratio', self.parent_grid_ratio)
        self.update_entry(nml_geogrid, 'i_parent_start', self.sx)
        self.update_entry(nml_geogrid, 'j_parent_start', self.sy)
        self.update_entry(nml_geogrid, 's_we', 1)
        self.update_entry(nml_geogrid, 's_sn', 1)
        self.update_entry(nml_geogrid, 'e_we', self.ex)
        self.update_entry(nml_geogrid, 'e_sn', self.ey)
        self.update_entry(nml_geogrid, 'geog_data_res', self.geog_res)
        self.update_entry(nml_geogrid, 'dx', self.dx)
        self.update_entry(nml_geogrid, 'dy', self.dy)
        self.update_entry(nml_geogrid, 'map_proj', 'lambert')
        self.update_entry(nml_geogrid, 'ref_lat', self.ref_lat)
        self.update_entry(nml_geogrid, 'ref_lon', self.ref_lon)
        self.update_entry(nml_geogrid, 'truelat1', self.truelats[0])
        self.update_entry(nml_geogrid, 'truelat2', self.truelats[1])
        self.update_entry(nml_geogrid, 'stand_lon', self.stand_lon)


    def update_inputnl(self, nml):
        """
        Update the WRF input namelist according to the domain configuration.

        :param nml: the namelist dictionary
        """
        nml_doms = nml['domains']
        self.update_entry(nml_doms, 'grid_id', self.dom_id)
        self.update_entry(nml_doms, 'parent_id', self.parent_id)
        self.update_entry(nml_doms, 'parent_grid_ratio', self.parent_grid_ratio)
        self.update_entry(nml_doms, 'i_parent_start', self.sx)
        self.update_entry(nml_doms, 'j_parent_start', self.sy)
        self.update_entry(nml_doms, 's_we', 1)
        self.update_entry(nml_doms, 's_sn', 1)
        self.update_entry(nml_doms, 'e_we', self.ex)
        self.update_entry(nml_doms, 'e_sn', self.ey)
        self.update_entry(nml_doms, 'dx', self.dx)
        self.update_entry(nml_doms, 'dy', self.dy)
        self.update_entry(nml_doms, 'sr_x', self.sr_x)
        self.update_entry(nml_doms, 'sr_y', self.sr_y)
        self.update_entry(nml_doms, 'time_step', self.time_step)
    
    
    def update_entry(self, section, key, value):
        """
        Update the corresponding entry (given by dom_id) in the given section.
        Crashes if preceding values are not filled out. i.e. if setting value for domain 2,
        the value for domain 1 must be already filled.

        :param section: the namelist section to process
        :param key: the name of the parameter
        :param value: the value of the parameter
        """
        entries = section[key] if key in section else []
        if type(entries) != list:
            entries = [entries]
        if len(entries) < self.dom_id - 2:
            raise ValueError('Cannot set namelist value for domain %d when previous domains not filled out.' % self.dom_id)
        if len(entries) <= self.dom_id - 1:
            entries.append(value)
        else:
            entries[self.dom_id-1] = value
        section[key] = entries


if __name__ == '__main__':
    import f90nml
    import sys
    d = WPSDomainLCC({'ref_lat' : 39, 'ref_lon' : -105.663, 'dx' : 1000, 'dy' : 1000, 'nx' : 51, 'ny' : 51, 'sg_x' : 50, 'sg_y' : 50})

    wps_nml = f90nml.read(sys.argv[1])

    d.update_wpsnl(wps_nml)
    f90nml.write(wps_nml, sys.argv[2], force=True)

    d = WPSDomainLCC({'precomputed' : 'geo_em.d01.nc'})
    wps_nml = f90nml.read(sys.argv[1])
    d.update_wpsnl(wps_nml)
    f90nml.write(wps_nml, sys.argv[3], force=True)




