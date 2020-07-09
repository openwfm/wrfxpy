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
from __future__ import print_function
import netCDF4
import numpy as np
import pyproj
from six.moves import range


class WPSDomainLCC(object):
    """
    A top-level Lambert Conic Conformal projection domain.
    """
    def __init__(self, dom_id, cfg, parent = None):
        """
        Initialize a top-level LCC domain.

        :param dom_id: the domain id
        :param cfg: the configuration dictionary
        :param parent: the parent domain (if this is a child)
        """
        self.dom_id = dom_id
        self.parent = parent
        self.parent_id = parent.dom_id if parent is not None else self.dom_id
        self.top_level = parent is None
        self.subgrid_ratio = cfg.get('subgrid_ratio', (1, 1))
        self.geog_res = str(cfg.get('geog_res', '30s'))
        self.history_interval = cfg.get('history_interval', 60)
        self.frames_per_outfile = cfg.get('frames_per_outfile', 24)
        if 'precomputed' in cfg:
            self._init_from_precomputed(cfg)
        else:
            self._init_from_dict(cfg)


    def _init_from_dict(self, cfg):
        """
        Initialize the domain from a dictionary containing the parameters.

        :param cfg: the parameter dictionary
        """
        if self.top_level:
            # this is a top level domain
            self.parent_cell_size_ratio = 1
            self.parent_time_step_ratio = 1
            self.cell_size = tuple(cfg['cell_size'])
            self.ref_latlon = tuple(cfg['center_latlon'])
            self.domain_size = tuple(cfg['domain_size'])
            self.time_step = cfg.get('time_step', int(max(6*min(self.domain_size) // 1000, 1)))
            self.parent_start = (1,1)
            if 'truelats' in cfg:
                self.truelats = cfg['truelats']
            else:
                self.truelats = [self.ref_lat]*2
            if 'stand_lon' in cfg:
                self.stand_lon = cfg['stand_lon']
            else:
                self.stand_lon = self.ref_lon
            self._init_projection()
        else:
            self.parent_cell_size_ratio = cfg['parent_cell_size_ratio']
            self.parent_time_step_ratio = cfg['parent_time_step_ratio']
            self.cell_size = [float(x) / self.parent_cell_size_ratio for x in self.parent.cell_size]
            if 'bounding_box' in cfg:
                # this is a child domain that we need to place dynamically given a bounding box
                self.bbox = tuple(cfg['bounding_box'])

                min_lon, min_lat, max_lon, max_lat = self.bbox

                # project the min/max coordinates
                i1, j1 = self.parent.latlon_to_ij(min_lat, min_lon)
                i2, j2 = self.parent.latlon_to_ij(max_lat, max_lon)

                self.parent_start = (int(i1), int(j1))
                self.parent_end = (int(i2+1), int(j2+1))

                if any([x < 0 for x in self.parent_start]) or \
                   self.parent_end[0] >= self.parent.domain_size[0] or \
                   self.parent_end[1] >= self.parent.domain_size[1]:
                    raise ValueError('Cannot place child, out of parents bounds')

                pstart, pend, pgr = self.parent_start, self.parent_end, self.parent_cell_size_ratio
                self.domain_size = (pgr * (pend[0] - pstart[0] + 1) + 1, pgr * (pend[1] - pstart[1] + 1) + 1)
            elif 'parent_start' in cfg and 'parent_end' in cfg:
                # this is a child domain, placed normally
                self.parent_start = tuple(cfg['parent_start'])
                self.parent_end = tuple(cfg['parent_end'])
                pstart, pend, pgr = self.parent_start, self.parent_end, self.parent_cell_size_ratio
                self.domain_size = (pgr * (pend[0] - pstart[0] + 1) + 1, pgr * (pend[1] - pstart[1] + 1) + 1)
            else:
                raise ValueError('Invalid domain specification')

        # either way, it's not precomputed
        self.precomputed = False


    def _init_from_precomputed(self, cfg):
        """
        Initialize the domain from an existing geo_em.dXX.nc file.

        :param path: path to the geo_em file
        """
        path = cfg['precomputed']
        self.precomputed_path = path
        self.history_interval = cfg.get('history_interval', 60)
        self.precomputed = True
        d = netCDF4.Dataset(path)
        self.subgrid_ratio = (int(d.getncattr('sr_x')), int(d.getncattr('sr_y')))
        self.geog_res = '30s' # not in geo_em file and does not matter (no processing will be done)
        self.dom_id = int(d.getncattr('grid_id'))
        self.parent_id = int(d.getncattr('parent_id'))
        self.parent_cell_size_ratio = int(d.getncattr('parent_grid_ratio'))
        if self.top_level:
            # this is a top-level domain
            self.domain_size = (int(d.getncattr('i_parent_end')), int(d.getncattr('j_parent_end')))
            self.cell_size = (float(d.getncattr('DX')), float(d.getncattr('DY')))
            self.ref_latlon = (float(d.getncattr('CEN_LAT')),float(d.getncattr('CEN_LON')))
            self.truelats = (float(d.getncattr('TRUELAT1')), float(d.getncattr('TRUELAT2')))
            self.stand_lon = float(d.getncattr('STAND_LON'))
            self.time_step = cfg.get('time_step', max(1, int(6*min(self.cell_size) // 1000)))
            self.parent_start = (1, 1)
            self.parent_time_step_ratio = 1
            self._init_projection()
        else:
            # this is a child domain
            self.parent_start = (int(d.getncattr('i_parent_start')), int(d.getncattr('j_parent_start')))
            self.parent_end = (int(d.getncattr('i_parent_end')), int(d.getncattr('j_parent_end')))
            # precomputed files don't remember grid size but rather directly parent end/start
            pstart, pend, pgr = self.parent_start, self.parent_end, self.parent_cell_size_ratio
            self.domain_size = (pgr * (pend[0] - pstart[0] + 1) + 1, pgr * (pend[1] - pstart[1] + 1) + 1)
            self.parent_time_step_ratio = cfg['parent_time_step_ratio']
            self.cell_size = tuple(float(x) / self.parent_cell_size_ratio for x in self.parent.cell_size)

        # only used if this is a top-level domain but we read it in anyway
        d.close()


    def _init_projection(self):
        """
        This function is based on code by Pavel Krc <krc@cs.cas.cz>, minor changes applied.
        It initializes the projection for a top-level domain.
        """
        radius = 6370e3
        
        # Spherical latlon used by WRF
        self.latlon_sphere = pyproj.Proj(proj='latlong',
                a=radius, b=radius, towgs84='0,0,0', no_defs=True)

        # Lambert Conformal Conic used by WRF
        self.lambert_grid = pyproj.Proj(proj='lcc',
            lat_1=self.truelats[0],
            lat_2=self.truelats[1],
            lat_0=self.ref_latlon[0],
            lon_0=self.stand_lon,
            a=radius, b=radius, towgs84='0,0,0', no_defs=True)

        grid_size_i = (self.domain_size[0] - 2) * self.cell_size[0]
        grid_size_j = (self.domain_size[1] - 2) * self.cell_size[1]

        grid_center_i, grid_center_j = pyproj.transform(
                self.latlon_sphere, self.lambert_grid,
                self.ref_latlon[1], self.ref_latlon[0])
        
        self.offset_i = grid_center_i - grid_size_i * .5
        self.offset_j = grid_center_j - grid_size_j * .5


    def latlon_to_ij(self, lat, lon):
        """
        Convert latitude and longitude into grid coordinates.

        If this is a child domain, it asks it's parent to do the projection and then
        remaps it into its own coordinate system via parent_start and cell size ratio.

        :param lat: latitude
        :param lon: longitude
        :return: the i, j position in grid coordinates
        """
        if self.top_level:
            proj_i, proj_j = pyproj.transform(self.latlon_sphere, self.lambert_grid,
                    lon, lat)
            return  ((proj_i - self.offset_i) / self.cell_size[0],
                     (proj_j - self.offset_j) / self.cell_size[1])
        else:
            pi, pj = self.parent.latlon_to_ij(lat, lon)
            pcsr, ps = self.parent_cell_size_ratio, self.parent_start
            return ((pi - ps[0] + 1.5) * pcsr - .5,
                    (pj - ps[1] + 1.5) * pcsr - .5)

    
    def ij_to_latlon(self, i, j):
        """
        Convert grid into latitude and longitude coordinates.

        If this is a child domain, it asks it's parent to do the projection and then
        remaps it into its own coordinate system via parent_start and cell size ratio.

        :param i: x grid coordinate
        :param j: y grid coordinate
        :return: the latitude,longitude position in degree coordinates
        """
        if self.top_level:
            lon, lat = pyproj.transform(self.lambert_grid, self.latlon_sphere,
                    i * self.cell_size[0] + self.offset_i,
                    j * self.cell_size[1] + self.offset_j)
            return lat, lon
        else:
            pcsr, ps = self.parent_cell_size_ratio, self.parent_start
            return self.parent.ij_to_latlon((i+.5)/pcsr+ps[0]-1.5, (j+.5)/pcsr+ps[1]-1.5) 

    def bounding_box(self):
        """
        Generate bounding box degree coordinates for the four corners of the domain.

        :return: the four latitude,longitude degree coordinates of the domain corners
        """
        latlon00 = self.ij_to_latlon(-1,-1)
        latlon01 = self.ij_to_latlon(-1,self.domain_size[1]+1)
        latlon11 = self.ij_to_latlon(self.domain_size[0]+1,self.domain_size[1]+1)
        latlon10 = self.ij_to_latlon(self.domain_size[0]+1,-1)
        return (latlon00,latlon01,latlon11,latlon10)
    
    def update_wpsnl(self, nml):
        """
        Update the share and geogrid section of the WPS namelist.

        :param nml_geogrid: a dictionary containing the geogrid section of the WPS namelist
        """
        nml_share = nml['share']
        self._update_entry(nml_share, 'subgrid_ratio_x', self.subgrid_ratio[0])
        self._update_entry(nml_share, 'subgrid_ratio_y', self.subgrid_ratio[1])

        # prevent geogrid from re-processing the grid (HACK: note that all grids must be activated
        # before metgrid runs!)
        self._update_entry(nml_share, 'active_grid', not self.precomputed)

        nml_geogrid = nml['geogrid']
        self._update_entry(nml_geogrid, 'geog_data_res', self.geog_res)
        self._update_entry(nml_geogrid, 'parent_id', self.parent_id)
        self._update_entry(nml_geogrid, 'parent_grid_ratio', self.parent_cell_size_ratio)
        self._update_entry(nml_geogrid, 'i_parent_start', self.parent_start[0])
        self._update_entry(nml_geogrid, 'j_parent_start', self.parent_start[1])
        self._update_entry(nml_geogrid, 's_we', 1)
        self._update_entry(nml_geogrid, 's_sn', 1)
        self._update_entry(nml_geogrid, 'e_we', self.domain_size[0])
        self._update_entry(nml_geogrid, 'e_sn', self.domain_size[1])

        # only for top-level domains
        if self.dom_id == self.parent_id:
            self._update_entry(nml_geogrid, 'dx', self.cell_size[0])
            self._update_entry(nml_geogrid, 'dy', self.cell_size[1])
            self._update_entry(nml_geogrid, 'map_proj', 'lambert')
            self._update_entry(nml_geogrid, 'ref_lat', self.ref_latlon[0])
            self._update_entry(nml_geogrid, 'ref_lon', self.ref_latlon[1])
            self._update_entry(nml_geogrid, 'truelat1', self.truelats[0])
            self._update_entry(nml_geogrid, 'truelat2', self.truelats[1])
            self._update_entry(nml_geogrid, 'stand_lon', self.stand_lon)


    def update_inputnl(self, nml):
        """
        Update the WRF input namelist according to the domain configuration.

        :param nml: the namelist dictionary
        """
        nml_tc = nml['time_control']
        self._update_entry(nml_tc, 'history_interval', self.history_interval)
        self._update_entry(nml_tc, 'frames_per_outfile', self.frames_per_outfile)
        self._update_entry(nml_tc, 'input_from_file', True)

        nml_doms = nml['domains']
        self._update_entry(nml_doms, 'grid_id', self.dom_id)
        self._update_entry(nml_doms, 'parent_id', self.parent_id)
        self._update_entry(nml_doms, 'i_parent_start', self.parent_start[0])
        self._update_entry(nml_doms, 'j_parent_start', self.parent_start[1])
        self._update_entry(nml_doms, 's_we', 1)
        self._update_entry(nml_doms, 's_sn', 1)
        self._update_entry(nml_doms, 'e_we', self.domain_size[0])
        self._update_entry(nml_doms, 'e_sn', self.domain_size[1])
        self._update_entry(nml_doms, 'sr_x', self.subgrid_ratio[0])
        self._update_entry(nml_doms, 'sr_y', self.subgrid_ratio[1])
        self._update_entry(nml_doms, 'dx', self.cell_size[0])
        self._update_entry(nml_doms, 'dy', self.cell_size[1])
        if self.dom_id == self.parent_id:
            # store cell size & time step
            self._update_entry(nml_doms, 'time_step', self.time_step)
        
        # for child domains cell sizes and timesteps are determined relative to parents
        # for top-level domains, these are fixed to 1
        self._update_entry(nml_doms, 'parent_grid_ratio', self.parent_cell_size_ratio)
        self._update_entry(nml_doms, 'parent_time_step_ratio', self.parent_time_step_ratio)

        # ensure boundary conditions are correctly read in (top level is specified, rest is nested)
        nml_bdy = nml['bdy_control']
        self._update_entry(nml_bdy, 'specified', self.top_level)
        self._update_entry(nml_bdy, 'nested', not self.top_level)
    
    
    def _update_entry(self, section, key, value):
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
            raise ValueError('Cannot set namelist value for domain %d, previous domains not filled out.' % self.dom_id)
        if len(entries) <= self.dom_id - 1:
            entries.append(value)
        else:
            entries[self.dom_id-1] = value
        section[key] = entries


class WPSDomainConf(object):
    """
    Represents a domain configuration for WPS, that is one or more domains with
    one top-level domain.
    """

    def __init__(self, cfg):
        """
        Initilize the domain configuration with a config dictionary mapping
        string domain id ("1", "2", ...) to its configuration.

        :param cfg: the configuration dictionary
        """
        self.domains = []

        # process domains in order
        for i in range(1, len(cfg)+1):
            this_dom = cfg[str(i)]
            par_dom = self.domains[this_dom['parent_id']-1] if 'parent_id' in this_dom else None
            self.domains.append(WPSDomainLCC(i, this_dom, par_dom))

    def __len__(self):
        """
        Returns the number of domains.
        """
        return len(self.domains)

    
    def prepare_for_geogrid(self, wps_nml, input_nml = None, wrfxpy_dir = None, wps_dir = None):
        """
        Update the namelists for geogrid processing - write the domain configurations.
        Additionally

        :param wps_nml: the WPS namelist
        :param input_nml: the input namelist, if None that skipped
        :param wrfxpy_dir: the installation directory of wrfxpy (if None, no linking is done)
        :param wps_dir: the WPS working directory
        """
        N = len(self.domains)
        wps_nml['share']['max_dom'] = N

        def extend_entry(nl, key):
            """
            Some keys should be same across all domains.  This function takes the first value
            and copies it max_dom times.
            """
            val = nl[key][0] if type(nl[key]) == list else nl[key]
            nl[key] = [val] * N

        # process the input namelist
        if input_nml is not None:
            doms = input_nml['domains']
            doms['max_dom'] = N
            doms['s_vert'] = [1] * N   
            extend_entry(doms, 'e_vert')

        for dom in self.domains:
            # ensure both namelists is updated with respect to all domains
            dom.update_wpsnl(wps_nml)
            if input_nml is not None:
                dom.update_inputnl(input_nml)

            # if the domains are precomputed, link in the files
            if dom.precomputed and wrfxpy_dir is not None:
                link_tgt = osp.join(wrfxpy_dir, dom.precomputed_path)
                link_loc = osp.join(wps_dir, 'geo_em.d%02d.nc' % dom.dom_id)
                symlink_unless_exists(link_tgt, link_loc) 


    def prepare_for_metgrid(self, wps_nml):
        """
        Set all domains that we use to active.

        :param wps_nml: the WPS namelist
        """
        wps_nml['share']['active_grid'] = [True] * len(self.domains)
        

if __name__ == '__main__':
    import f90nml
    import sys
    import json

    # parse a JSON domain configuration and a namelist.wps file and build the correct wps
    if len(sys.argv) != 3 and len(sys.argv) != 4:
        print(('usage: %s <domains_json_file> <namelist.wps> [namelist.input]' % sys.argv[0]))
        sys.exit(1)

    dcfg = WPSDomainConf(json.load(open(sys.argv[1])))
    
    wps_nml = f90nml.read(sys.argv[2])
    wrf_nml = f90nml.read(sys.argv[3]) if len(sys.argv) == 4 else None
    dcfg.prepare_for_geogrid(wps_nml, wrf_nml)

    f90nml.write(wps_nml, sys.argv[2], force=True)
    if wrf_nml is not None:
        f90nml.write(wps_nml, sys.argv[3], force=True)


