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
            self.init_from_precomputed(cfg)
        else:
            self.init_from_dict(cfg)


    def init_from_dict(self, cfg):
        """
        Initialize the domain from a dictionary containing the parameters.

        :param cfg: the parameter dictionary
        """
        self.parent_id = cfg.get('parent_id', self.dom_id)
        self.subgrid_ratio = cfg.get('subgrid_ratio', (1, 1))
        self.geog_res = str(cfg.get('geog_res', '30s'))
        self.history_interval = cfg.get('history_interval', 60)
        if self.parent_id == self.dom_id:
            # this is a top level domain
            self.parent_cell_size_ratio = 1
            self.parent_time_step_ratio = 1
            self.cell_size = tuple(cfg['cell_size'])
            self.ref_latlon = tuple(cfg['center_latlon'])
            self.domain_size = tuple(cfg['domain_size'])
            self.time_step = cfg.get('time_step', int(max(6*min(self.domain_size)/1000, 1)))
            self.parent_start = (1,1)
            if 'truelats' in cfg:
                self.truelats = cfg['truelats']
            else:
                self.truelats = [self.ref_lat]*2
            if 'stand_lon' in cfg:
                self.stand_lon = cfg['stand_lon']
            else:
                self.stand_lon = self.ref_lon
        elif 'center_latlon' not in cfg:
            # this is a child domain, placed normally
            self.parent_cell_size_ratio = cfg['parent_cell_size_ratio']
            self.parent_time_step_ratio = cfg['parent_time_step_ratio']
            self.parent_start = tuple(cfg['parent_start'])
            self.parent_end = tuple(cfg['parent_end'])
            pstart, pend, pgr = self.parent_start, self.parent_end, self.parent_cell_size_ratio
            self.domain_size = (pgr * (pend[0] - pstart[0] + 1) + 1, pgr * (pend[1] - pstart[1] + 1) + 1)
        else:
            # this is a child domain that we need to place dynamically
            ref_latlon = tuple(cfg['center_latlon'])
            radius = cfg['radius']



        # either way, it's not precomputed
        self.precomputed = False


    def init_from_precomputed(self, cfg):
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
        if self.dom_id == self.parent_id:
            # this is a top-level domain
            self.domain_size = (int(d.getncattr('i_parent_end')), int(d.getncattr('j_parent_end')))
            self.cell_size = (float(d.getncattr('DX')), float(d.getncattr('DY')))
            self.ref_latlon = (float(d.getncattr('CEN_LAT')),float(d.getncattr('CEN_LON')))
            self.truelats = (float(d.getncattr('TRUELAT1')), float(d.getncattr('TRUELAT2')))
            self.stand_lon = float(d.getncattr('STAND_LON'))
            self.time_step = cfg.get('time_step', max(1, int(6*min(self.cell_size)/1000)))
            self.parent_start = (1, 1)
            self.parent_time_step_ratio = 1
        else:
            # this is a child domain
            self.parent_start = (int(d.getncattr('i_parent_start')), int(d.getncattr('j_parent_start')))
            self.parent_end = (int(d.getncattr('i_parent_end')), int(d.getncattr('j_parent_end')))
            # precomputed files don't remember grid size but rather directly parent end/start
            pstart, pend, pgr = self.parent_start, self.parent_end, self.parent_cell_size_ratio
            self.domain_size = (pgr * (pend[0] - pstart[0] + 1) + 1, pgr * (pend[1] - pstart[1] + 1) + 1)
            self.parent_time_step_ratio = cfg['parent_time_step_ratio']

        # only used if this is a top-level domain but we read it in anyway
        d.close()

    
    def update_wpsnl(self, nml):
        """
        Update the share and geogrid section of the WPS namelist.

        :param nml_geogrid: a dictionary containing the geogrid section of the WPS namelist
        """
        nml_share = nml['share']
        self.update_entry(nml_share, 'subgrid_ratio_x', self.subgrid_ratio[0])
        self.update_entry(nml_share, 'subgrid_ratio_y', self.subgrid_ratio[1])

        # prevent geogrid from re-processing the grid (HACK: note that all grids must be activated
        # before metgrid runs!)
        self.update_entry(nml_share, 'active_grid', not self.precomputed)

        nml_geogrid = nml['geogrid']
        self.update_entry(nml_geogrid, 'geog_data_res', self.geog_res)
        self.update_entry(nml_geogrid, 'parent_id', self.parent_id)
        self.update_entry(nml_geogrid, 'parent_grid_ratio', self.parent_cell_size_ratio)
        self.update_entry(nml_geogrid, 'i_parent_start', self.parent_start[0])
        self.update_entry(nml_geogrid, 'j_parent_start', self.parent_start[1])
        self.update_entry(nml_geogrid, 's_we', 1)
        self.update_entry(nml_geogrid, 's_sn', 1)
        self.update_entry(nml_geogrid, 'e_we', self.domain_size[0])
        self.update_entry(nml_geogrid, 'e_sn', self.domain_size[1])

        # only for top-level domains
        if self.dom_id == self.parent_id:
            self.update_entry(nml_geogrid, 'dx', self.cell_size[0])
            self.update_entry(nml_geogrid, 'dy', self.cell_size[1])
            self.update_entry(nml_geogrid, 'map_proj', 'lambert')
            self.update_entry(nml_geogrid, 'ref_lat', self.ref_latlon[0])
            self.update_entry(nml_geogrid, 'ref_lon', self.ref_latlon[1])
            self.update_entry(nml_geogrid, 'truelat1', self.truelats[0])
            self.update_entry(nml_geogrid, 'truelat2', self.truelats[1])
            self.update_entry(nml_geogrid, 'stand_lon', self.stand_lon)


    def update_inputnl(self, nml):
        """
        Update the WRF input namelist according to the domain configuration.

        :param nml: the namelist dictionary
        """
        nml_tc = nml['time_control']
        self.update_entry(nml_tc, 'history_interval', self.history_interval)

        nml_doms = nml['domains']
        self.update_entry(nml_doms, 'grid_id', self.dom_id)
        self.update_entry(nml_doms, 'parent_id', self.parent_id)
        self.update_entry(nml_doms, 'i_parent_start', self.parent_start[0])
        self.update_entry(nml_doms, 'j_parent_start', self.parent_start[1])
        self.update_entry(nml_doms, 's_we', 1)
        self.update_entry(nml_doms, 's_sn', 1)
        self.update_entry(nml_doms, 'e_we', self.domain_size[0])
        self.update_entry(nml_doms, 'e_sn', self.domain_size[1])
        self.update_entry(nml_doms, 'sr_x', self.subgrid_ratio[0])
        self.update_entry(nml_doms, 'sr_y', self.subgrid_ratio[1])
        if self.dom_id == self.parent_id:
            # store cell size & time step
            self.update_entry(nml_doms, 'dx', self.cell_size[0])
            self.update_entry(nml_doms, 'dy', self.cell_size[1])
            self.update_entry(nml_doms, 'time_step', self.time_step)
        
        # for child domains cell sizes and timesteps are determined relative to parents
        # for top-level domains, these are fixed to 1
        self.update_entry(nml_doms, 'parent_grid_ratio', self.parent_cell_size_ratio)
        self.update_entry(nml_doms, 'parent_time_step_ratio', self.parent_time_step_ratio)
    
    
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
    d = WPSDomainLCC({'ref_latlon' : [39, -105.663],
                      'cell_size' : [1000, 1000],
                      'domain_size' : [51, 51],
                      'subgrid_ratio' : [50, 50]})

    wps_nml = f90nml.read(sys.argv[1])

    d.update_wpsnl(wps_nml)
    f90nml.write(wps_nml, sys.argv[2], force=True)

    d = WPSDomainLCC({'precomputed' : 'geo_em.d01.nc', 'time_step' : 6})
    wps_nml = f90nml.read(sys.argv[1])
    d.update_wpsnl(wps_nml)
    f90nml.write(wps_nml, sys.argv[3], force=True)




